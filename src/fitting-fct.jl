
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
    #m2    = par[9]
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
function tail_model(x, par)
    scale = par[1]
    σ     = par[2]
    μ     = par[3]
    C     = par[4]
    l     = par[5]
    e     = par[6]
    t     = par[7]
    d     = 0
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + C * (exp(e*(x-l)) + d) / (exp((x-l)/t) +1)
end
function tail_model2(par, x)
    scale = par.scale
    σ     = par.σ
    μ     = par.µ
    C     = par.C
    l     = par.l
    e     = par.e
    t     = par.t
    d     = 0
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + C * (exp(e*(x-l)) + d) / (exp((x-l)/t) +1)
end
@. linmodel(x, p) = p[1]*x + p[2]
@. hypmodel_simple(x, p) = sqrt(p[1]+p[2]/x)
@. hypmodel(x, p) = sqrt(p[1]+p[2]/(x^2))
@. poly3model(x, p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
@. atanmodel(x, p) = p[1] * atan(p[2]*x + p[3]) + p[4]

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