function bat_fit(hist, fit_function, prior; nsamples=10^5, nchains=4, params_first=true, norm_expected=true)
    likelihood = let h = hist, f = fit_function
        # Histogram counts for each bin as an array:
        observed_counts = h.weights

        # Histogram binning:
        bin_edges = h.edges[1]
        bin_edges_left = bin_edges[1:end-1]
        bin_edges_right = bin_edges[2:end]
        bin_widths = bin_edges_right - bin_edges_left
        bin_centers = midpoints(bin_edges)

        params -> begin
            if params_first
                expected_counts = bin_widths .* f(params, bin_centers)
            else
                expected_counts = bin_widths .* f(bin_centers, params)
            end
            if norm_expected
                expected_counts ./= sum(expected_counts)
                expected_counts .*= sum(observed_counts)
            end
            ll_value = sum(logpdf.(Poisson.(expected_counts), observed_counts))
            return LogDVal(ll_value)
        end
    end
    parshapes = varshape(prior)
    posterior = PosteriorDensity(likelihood, prior);
    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = nsamples, nchains = nchains)).result;
end

#
# Wrapper for fits that might need several iterations
#
function bat_fit_wrapper(histogram, model, prior, nsamples_list; ch_str="XX", params_first=false, norm_expected=false)
    parameters, errors = [], []
    run = 1
    for ns in nsamples_list
        @info("Ch" * ch_str * " | try " * string(run) * " of " * string(length(nsamples_list)) * ": " * string(ns) * " samples")
        samples = try bat_fit(histogram, model, prior; nsamples=ns, nchains=4, params_first=params_first, norm_expected=norm_expected);
        catch
            IJulia.clear_output(true)
            run += 1
            continue
        end
        IJulia.clear_output(true)
        @info("Success: Ch" * ch_str * " | " * string(ns) * " samples")
        parameters = mode(samples)
        errors = std(samples)
        break
    end
    return [parameters, errors]
end


#
# Peak model with low E tail
#
function peak(E, par)
    n = try par[1] catch; par[1][1] end
    σ = try par[2] catch; par[2][1] end
    µ = try par[3] catch; par[3][1] end
    
    bkg_l = try par[4] catch; par[4][1] end
    bkg_r = try par[5] catch; par[5][1] end
    b = (bkg_r - bkg_l) / (maximum(E) - minimum(E))
    a = bkg_l - b * minimum(E)
    
    c = try par[6] catch; par[6][1] end
    d = try par[7] catch; par[7][1] end
    β = try par[8] catch; par[8][1] end
    
    d = b > 0 ? -d : d
    correction = d < 0 ? -1*minimum(d/2 .* erfc.((E .- µ) ./ (sqrt(2)*σ))) : 0.0
    result = @. n / sqrt(2*pi*σ^2) * exp(-(E - µ)^2 / (2*σ^2)) + # Gaussian
        a + b*E + # Linear bkg
        d/2 * erfc((E - µ) / (sqrt(2)*σ)) + correction +  # Low-E step
        c/(2*β) * exp((E-µ) / β + σ^2 / (2*β^2)) * erfc((E - µ) / (sqrt(2)*σ) + σ / (sqrt(2)*β)) # Low-E tail
end
function peak_twisted(par, E)
    n = try par[1] catch; par[1][1] end
    σ = try par[2] catch; par[2][1] end
    µ = try par[3] catch; par[3][1] end
    
    bkg_l = try par[4] catch; par[4][1] end
    bkg_r = try par[5] catch; par[5][1] end
    b = (bkg_r - bkg_l) / (maximum(E) - minimum(E))
    a = bkg_l - b * minimum(E)
    
    c = try par[6] catch; par[6][1] end
    d = try par[7] catch; par[7][1] end
    β = try par[8] catch; par[8][1] end
    
    d = b > 0 ? -d : d
    correction = d < 0 ? -1*minimum(d/2 .* erfc.((E .- µ) ./ (sqrt(2)*σ))) : 0.0
    result = @. n / sqrt(2*pi*σ^2) * exp(-(E - µ)^2 / (2*σ^2)) + # Gaussian
        a + b*E + # Linear bkg
        d/2 * erfc((E - µ) / (sqrt(2)*σ)) + correction +  # Low-E step
        c/(2*β) * exp((E-µ) / β + σ^2 / (2*β^2)) * erfc((E - µ) / (sqrt(2)*σ) + σ / (sqrt(2)*β)) # Low-E tail
end

function broad_peak(E, par)
    n = try par[1] catch; par[1][1] end
    σ = try par[2] catch; par[2][1] end
    µ = try par[3] catch; par[3][1] end
    
    bkg_l = try par[4] catch; par[4][1] end
    bkg_r = try par[5] catch; par[5][1] end
    b = (bkg_r - bkg_l) / (maximum(E) - minimum(E))
    a = bkg_l - b * minimum(E)
    
    c = try par[6] catch; par[6][1] end
    d = try par[7] catch; par[7][1] end
    β = try par[8] catch; par[8][1] end
    
    
    n2 = try par[9]  catch; par[9][1] end
    σ2 = try par[10] catch; par[10][1] end
    
    d = b > 0 ? -d : d
    correction = d < 0 ? -1*minimum(d/2 .* erfc.((E .- µ) ./ (sqrt(2)*σ))) : 0.0
    result = @. n / sqrt(2*pi*σ^2) * exp(-(E - µ)^2 / (2*σ^2)) + # Gaussian
        n2 / sqrt(2*pi*σ2^2) * exp(-(E - µ)^2 / (2*σ2^2)) + # Gaussian 2
        a + b*E + # Linear bkg
        d/2 * erfc((E - µ) / (sqrt(2)*σ)) .+ correction +  # Low-E step
        c/(2*β) * exp((E-µ) / β + σ^2 / (2*β^2)) * erfc((E - µ) / (sqrt(2)*σ) + σ / (sqrt(2)*β)) # Low-E tail
end

function tail_model(x, par)
    scale = try par[1] catch; par[1][1] end
    σ     = try par[2] catch; par[2][1] end
    μ     = try par[3] catch; par[3][1] end
    C     = try par[4] catch; par[4][1] end
    l     = try par[5] catch; par[5][1] end
    e     = try par[6] catch; par[6][1] end
    t     = try par[7] catch; par[7][1] end
    d     = 0
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + C * (exp(e*(x-l)) + d) / (exp((x-l)*t) +1)
end




function peak_model(x, par)
    scale = par[1]
    σ     = par[2]
    μ     = par[3]
    
    #step = x <= µ ? par[4] : par[5]
    if length(x) > 1
        step = []
        for i in x
            push!(step, i <= µ ? par[4] : par[5])
        end
        #mwa_step = isodd(round(length(x)/4)) ? Int.(round(length(x)/4)) : Int.(round(length(x)/4)-1)
        s = 4 * par[2] / 0.1
        mwa_step = isodd(Int.(round(s))) ? Int.(round(s)) : Int.(round(s)) + 1
        step = mwa(Array{Float64,1}(step),mwa_step)
    else
        step = x <= µ ? par[4] : par[5]
    end
    return scale .* exp.(-0.5 * ((x .- μ).^2) / (σ^2)) ./ (sqrt(2 * π * σ^2)) .+ step
end
function peak_model2(par, x)
    scale = par.scale
    σ     = par.σ
    μ     = par.µ
    cp0   = par.cp0
    #m0    = par.m0
    cp1   = par.cp1
    #m1    = par.m1
#     s     = par.s
    
    if length(x) > 1
        step = []
        for i in x
            push!(step, i <= µ ? cp0 : cp1)
        end
        s = 4 * σ / 0.1
#         mwa_step = isodd(round(length(x)/4)) ? Int.(round(length(x)/4)) : Int.(round(length(x)/4)-1)
        mwa_step = isodd(Int.(round(s))) ? Int.(round(s)) : Int.(round(s)) + 1
        step = mwa(Array{Float64,1}(step),mwa_step)
    else
        step = x <= µ ? cp0 : cp1
    end
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + step
end

function model(x, par)
    scale = par[1]
    σ     = par[2]
    μ     = par[3]
    cp0   = par[4]
    m1    = par[5]
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + m1 * x + cp0
end
function model2(par, x)
    scale = par.scale
    σ     = par.σ
    μ     = par.µ
    cp0   = par.cp0
    m1    = par.m1
    #m2    = par[9]
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + m1 * x + cp0
end
function double_model(x, par)
    scale = par[1]
    σ     = par[2]
    μ     = par[3]
    scale2 = par[4]
    σ2     = par[5]
    μ2     = par[6]
    cp0   = par[7]
    m1    = par[8]
    #m2    = par[9]
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + scale2 * exp(-0.5 * ((x - μ2)^2) / (σ2^2)) / (sqrt(2 * π * σ2^2)) + m1 * x + cp0# + m2 * x
end
function double_model2(par,x)
    scale = par.scale
    σ     = par.σ
    μ     = par.µ
    scale2 = par.scale2
    σ2     = par.σ2
    μ2     = par.µ2
    cp0   = par.cp0
    m1    = par.m1
    #m2    = par[9]
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + scale2 * exp(-0.5 * ((x - μ2)^2) / (σ2^2)) / (sqrt(2 * π * σ2^2)) + m1 * x + cp0# + m2 * x
end
# function tail_model(x, par)
#     scale = par[1]
#     σ     = par[2]
#     μ     = par[3]
#     C     = par[4]
#     l     = par[5]
#     e     = par[6]
#     t     = par[7]
#     d     = 0
#     return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + C * (exp(e*(x-l)) + d) / (exp((x-l)/t) +1)
# end
# function tail_model2(par, x)
#     scale = par.scale
#     σ     = par.σ
#     μ     = par.µ
#     C     = par.C
#     l     = par.l
#     e     = par.e
#     t     = par.t
#     d     = 0
#     return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + C * (exp(e*(x-l)) + d) / (exp((x-l)/t) +1)
# end
@. linmodel(x, p) = p[1]*x + p[2]
@. hypmodel_simple(x, p) = sqrt(p[1]+p[2]/x)
@. hypmodel(x, p) = sqrt(p[1]+p[2]/(x^2))
@. poly3model(x, p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
@. atanmodel(x, p) = p[1] * atan(p[2]*x + p[3]) + p[4]
@. sqrt_fct(x,p) = sqrt( p[1] + p[2] * x + p[3] * x * x )

function fit_bat(hist, fit_function, prior, fit_ranges; nsamples=1*10^5, nchains=4, second_ff=fit_function)

    likelihood = let h = hist, f = fit_function
        # Histogram counts for each bin as an array:
        observed_counts = h.weights

        # Histogram binning:
        bin_edges = h.edges[1]
        bin_edges_left = bin_edges[1:end-1]
        bin_edges_right = bin_edges[2:end]
        bin_widths = bin_edges_right - bin_edges_left
        bin_centers = (bin_edges_right + bin_edges_left) / 2

        params -> begin
            # Log-likelihood for a single bin:
            function bin_log_likelihood(i)
                # Simple mid-point rule integration of fit function `f` over bin:
                expected_counts = bin_widths[i] * f(params, bin_centers[i])
                logpdf(Poisson(expected_counts), observed_counts[i])
            end

            # Sum log-likelihood over bins:
            idxs = eachindex(observed_counts)
            ll_value = bin_log_likelihood(idxs[1])
            for i in idxs[2:end]
                ll_value += bin_log_likelihood(i)
            end

            # Wrap `ll_value` in `LogDVal` so BAT knows it's a log density-value.
            return LogDVal(ll_value)
        end
    end
    parshapes = varshape(prior)
    posterior = PosteriorDensity(likelihood, prior);
    
    pretunesamples = 10000
    max_ncycles = 30
    T = Float64
    tuning = BAT.AdaptiveMetropolisTuning(
        r = 1.0,
        λ = 0.5,
        α = 0.05..0.9,
        β = 1.5,
        c = 1e-4..1e2
    )
    convergence = BAT.BrooksGelmanConvergence(
        threshold = T(5),
        corrected = false
    )
    init = BAT.MCMCInitStrategy(
        init_tries_per_chain = 8..128,
        max_nsamples_init = pretunesamples,
        max_nsteps_init = pretunesamples * 10,
        max_time_init = Inf
    )
    burnin = BAT.MCMCBurninStrategy(
        max_nsamples_per_cycle = pretunesamples,
        max_nsteps_per_cycle = pretunesamples * 10,
        max_time_per_cycle = Inf,
        max_ncycles = max_ncycles
    )
    parshapes = varshape(prior)
    posterior = PosteriorDensity(likelihood, prior);
    nsamples = 1*10^4
    nchains = 4
    
    
    samples = bat_sample(Philox4x((123, 456)), posterior, (nsamples, nchains), MetropolisHastings(),
        max_nsteps = 10 * nsamples,
        max_time = Inf,
        tuning = tuning,
        init = init,
        burnin = burnin,
        convergence = convergence,
    ).result;
    params = mode(samples)[1]
    error = std(samples)[1]
    # println("Truth: $p0")
    println("Mode: $(mode(samples))")
    println("Stddev: $(std(samples))")

    fitf = RadiationSpectra.FitFunction{Float64}( second_ff, 1, length(params));
    # set_parameter_bounds!(fitf, bounds)
    RadiationSpectra.set_fitranges!(fitf, (fit_ranges, ) )
    
    par = []
    err = []
    for x in eachindex(params)
        push!(par, params[x])
        push!(err, error[x])
    end
#     params = [params.scale[1], params.σ[1], params.µ[1], params.scale2[1], params.σ2[1], params.µ2[1], params.cp0[1], params.m1[1]]
#     error = [error.scale[1], error.σ[1], error.µ[1], error.scale2[1], error.σ2[1], error.µ2[1], error.cp0[1], error.m1[1]]
    fitf.fitted_parameters = par
    
    return par, err, fitf
end