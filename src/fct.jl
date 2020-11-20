#########################################
# Create mapping dictionary ############
sim_to_channel = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        channels[mapping[ged]["sim"]] = [mapping[ged]["channel"], mapping[ged]["name"]]
    end
    return channels
end;

#########################################
# Create mapping dictionary ############
channel_to_name = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        channels[mapping[ged]["channel"]] = mapping[ged]["name"]
    end
    return channels
end;

#########################################
# Create mapping dictionary ############
parameters = open(parameters_file) do f
    temp = JSON.parse(f)
    temp["GTF45"] = temp["GTF45_2"]
    return temp
end;

function root_to_typedtable(filename; hits_threshold=0.005, number_of_events=0, E_lim=0.3, mult_lim=1) # 0 = all
    file = TFile(filename)

    tree = file["fTree"];

    tt = Table(eventnumber = tree.eventnumber[:],
        hits_iddet   = tree.hits_iddet[:],
        hits_totnum  = tree.hits_totnum[:],
        hits_edep    = tree.hits_edep[:],
        hits_xpos    = tree.hits_xpos[:] .* 10, # mm
        hits_ypos    = tree.hits_ypos[:] .* 10, # mm
        hits_zpos    = tree.hits_zpos[:] .* 10) # mm

    #########################################
    # Filter events without energy deposition
    tt = tt |> @filter(_.hits_totnum != 0) |> Table;
    tt = tt |> @filter(sum(_.hits_edep) >= E_lim) |> Table;

    i = 1
    if number_of_events != 0
        i_max = number_of_events
    else
        i_max = size(tt, 1)
    end
    

    evt_num      = []
    multiplicity = []
    iddet        = []
    totnum       = []
    edep         = []
    position     = []

    p = Progress(i_max, dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=10)
    @info("Number of events = "*string(i_max))
    @info("Separating events by detector id and filter by minimum energy deposition threshold = "*string(hits_threshold)*" MeV")
    #########################################
    # Separate events by detector ###########
    while i <= i_max
        multi = 0
        for det in unique(tt[i].hits_iddet)
            index = findall(x->x==det, tt[i].hits_iddet)
            if sum( tt[i].hits_edep[ index ] ) >= hits_threshold
                multi += 1
                upside_down = 1
                if parameters[sim_to_channel[det][2]]["upside_down"] == true
                    upside_down = -1
                end

                hits_iddet_mapped = []
                for i in index
                    push!(hits_iddet_mapped, sim_to_channel[det][1])
                end
                push!(iddet, Array{Int64,1}(hits_iddet_mapped))
                push!(edep,  tt[i].hits_edep[ index ])
                push!(position, [ SVector{3}(([ tt[i].hits_xpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_x"], 
                                                tt[i].hits_ypos[k] .- parameters[sim_to_channel[det][2]]["detcenter_y"], 
                                                upside_down .* (tt[i].hits_zpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_z"] .+ upside_down * parameters[sim_to_channel[det][2]]["height"]/2) 
                                    ] * u"mm")...) for k in index ])
            end        
        end
        
        j = 1
        while j <= multi
            push!(evt_num, tt[i].eventnumber)
            push!(totnum, tt[i].hits_totnum)
            push!(multiplicity, multi)
            j += 1
        end
        next!(p)
        i += 1
    end
    tt = nothing
    #########################################
    # Create output table ###################
    filter_index = findall(x -> x == mult_lim, multiplicity)

    events = Table(evtno = evt_num[filter_index], 
        event_mult       = multiplicity[filter_index], 
        detno            = iddet[filter_index],
        hits_totnum      = totnum[filter_index],
        edep             = VectorOfArrays(edep[filter_index] .*1000 *u"keV"),
        pos              = position[filter_index]
    )
    return events
end




function read_AE_params(filename::String)
    fit_params = JSON.parsefile(filename)
    tt = Table(detector=[], set=[], energy=[], scale=[], sig=[], mu=[], C=[], l=[], e=[], t=[], scale_err=[], sig_err=[], mu_err=[], C_err=[], l_err=[], e_err=[], t_err=[])
    for det in keys(fit_params) # Detector
        for set in keys(fit_params[det])
            for E in keys(fit_params[det][set]) 
                if "bat_params" in keys(fit_params[det][set][E])
                    append!(tt, Table(detector=[det], set=[set], energy=[parse(Int64, E)], 
                        scale=[fit_params[det][set][E]["bat_params"][1][1]], 
                        sig=[fit_params[det][set][E]["bat_params"][1][2]], 
                        mu=[fit_params[det][set][E]["bat_params"][1][3]], 
                        C=[fit_params[det][set][E]["bat_params"][1][4]], 
                        l=[fit_params[det][set][E]["bat_params"][1][5]], 
                        e=[fit_params[det][set][E]["bat_params"][1][6]], 
                        t=[fit_params[det][set][E]["bat_params"][1][7]], 
                        scale_err=[fit_params[det][set][E]["bat_uncert"][1][1]], 
                        sig_err=[fit_params[det][set][E]["bat_uncert"][1][2]], 
                        mu_err=[fit_params[det][set][E]["bat_uncert"][1][3]], 
                        C_err=[fit_params[det][set][E]["bat_uncert"][1][4]], 
                        l_err=[fit_params[det][set][E]["bat_uncert"][1][5]], 
                        e_err=[fit_params[det][set][E]["bat_uncert"][1][6]], 
                        t_err=[fit_params[det][set][E]["bat_uncert"][1][7]]));
                end
            end
        end
    end
    return tt
end

function get_cut_values(fit_params, ch, set)
    fit_cut_values = []
    fit_survival =[]
    fit_survival_err =[]
    for cut in keys(fit_params[channel_to_name[ch]]["DEP_cuts"][set])
        if !(cut in ["cut_sum", "cut_fit"])
            push!(fit_cut_values, parse(Float64, cut))
            push!(fit_survival, fit_params[channel_to_name[ch]]["DEP_cuts"][set][cut]["survival"])
            push!(fit_survival_err, fit_params[channel_to_name[ch]]["DEP_cuts"][set][cut]["survival_err"])
        end
    end
    sum_cut_values = fit_params[channel_to_name[ch]]["DEP_cuts"][set * "_sum"]["cut_values"]
    sum_survival = fit_params[channel_to_name[ch]]["DEP_cuts"][set * "_sum"]["survival"]
    sum_survival_err = fit_params[channel_to_name[ch]]["DEP_cuts"][set * "_sum"]["survival_err"]
    return fit_cut_values, fit_survival, fit_survival_err, sum_cut_values, sum_survival, sum_survival_err
end
function get_cut(fit_params, ch, set)
    return fit_params[channel_to_name[ch]]["DEP_cuts"][set]["cut_fit"], fit_params[channel_to_name[ch]]["DEP_cuts"][set]["cut_sum"]
end


function find_intersect(samples, threshold::Real, noisefilter = 1)
    if noisefilter < 1 noisefilter = 1 end
    if length(samples) < 1 error("Empty array") end
    i = 1
    x = samples[i]
    counter = 0
    intersect = i
    if x != threshold
        findHigher = x < threshold
        while i < length(samples)
            i += 1
            x = samples[i]
            if (findHigher && (x >= threshold)) || (!findHigher && (x <= threshold))
                if counter == 0
                    intersect = i
                end
                if counter >= noisefilter-1
                    return intersect
                else
                    counter += 1
                end
            else
                counter = 0
            end
        end
        #println("No intersect found")
        return NaN
    end
    return NaN
end

simdir_old      = "/remote/ceph/group/gerda/data/simulation/gerda-mage-sim"
mapping_file    = "$simdir_old/UTILS/det-data/ged-mapping.json"
parameters_file = "$simdir_old/UTILS/det-data/ged-parameters.json"


#>------- Moving Window Average -----------------------<#
#>-----------------------------------------------------<#

function movingaverage(X::Vector,numofele::Int)
    BackDelta = div(numofele,2) 
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(X)
    Y = similar(X)
    for n = 1:len
        lo = max(1,n - BackDelta)
        hi = min(len,n + ForwardDelta)
        Y[n] = mean(X[lo:hi])
    end
    return Y
end

function movingaverage(X::Vector,numofele::Int,iterations::Int)
    i = 1
    while i <= iterations
        X = movingaverage(X,numofele) 
        i += 1
    end
    return X    
end

function mwa(samples::Vector, window::Integer)
#window length is in number of samples and has to be an odd number

    if !isodd(window)
        error("For this method the window size has to be an odd integer")
    end

    newsamples = zeros(eltype(samples), length(samples))
    p = convert(eltype(window), (window - 1) / 2)

	#padding before and after the original vector
    samples = vcat(ones(eltype(samples), p) .* samples[1],
		samples, ones(eltype(samples), p) .* samples[length(samples)])

    newsamples[1] = sum(samples[1:window])
    for i in 2:length(newsamples)
        newsamples[i] = newsamples[i-1] + samples[i+2*p] - samples[i-1]
    end

    newsamples./window #return after normalization

end

function multi_mwa(samples::Vector, window::Integer, reps::Integer)
    for i in 1:1:reps
        samples = mwa(samples, window)
    end
    return samples
end


#>------- Get average maximum -------------------------<#
#>-----------------------------------------------------<#
function get_avg_maximum(pulse, intervall)
    index = findall(x -> x == maximum(pulse), pulse)[1]
    intervall = (index - Int(floor(intervall/2))+1):1:(index + Int(ceil(intervall/2)))
    return mean(pulse[intervall])
end


#>------- Create mapping dictionary -------------------<#
#>-----------------------------------------------------<#
sim_to_channel = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        channels[mapping[ged]["sim"]] = [mapping[ged]["channel"], mapping[ged]["name"]]
    end
    return channels
end;

#>------- Create mapping dictionary -------------------<#
#>-----------------------------------------------------<#
channel_to_name = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        channels[mapping[ged]["channel"]] = mapping[ged]["name"]
    end
    return channels
end;

#>------- Create mapping dictionary -------------------<#
#>-----------------------------------------------------<#
parameters = open(parameters_file) do f
    temp = JSON.parse(f)
    temp["GTF45"] = temp["GTF45_2"]
    return temp
end;


#>------- Create Coax dictionary ----------------------<#
#>-----------------------------------------------------<#
channel_to_coax = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        if mapping[ged]["channel"] in [8,9,10,27,28,29,36,37,38,39]
            channels[mapping[ged]["channel"]] = true
        else
            channels[mapping[ged]["channel"]] = false
        end
    end
    return channels
end;


#>------- Create BEGe dictionary ----------------------<#
#>-----------------------------------------------------<#
channel_to_bege = open(mapping_file) do f
    mapping  = JSON.parse(f)["mapping"]
    channels = Dict()
    for ged in keys(mapping)
        if mapping[ged]["channel"] in [8,9,10,27,28,29,36,37,38,39]
            channels[mapping[ged]["channel"]] = false
        else
            channels[mapping[ged]["channel"]] = true
        end
    end
    return channels
end;




function check_parameters(tt, td, GBP, tau, Cd, Cf)
    cut_events = Table( waveform    = [],
                        evtno       = [],
                        detno       = [],
                        det_ch      = [],
                        edep        = []);
    for i in 1:1:size(tt,1)
        pulse = tt[i].waveform.value
        filtered_pulse = applyElectronics(pulse; Ts = 10e-9, GBP = GBP, tau = tau, Kv = 150e3, Cd = Cd, Cf = Cf, Rf = 500e6)
        inter = find_intersect(filtered_pulse, maximum(filtered_pulse)/2, 5)
        if inter + 399 < length(pulse) && inter - 400 > 0
            cut  = filtered_pulse[(inter-400):1:(inter+400-1)]
            time = (0:10:7990)u"ns"
            cut_pulse = RDWaveform(time, cut)

            append!(cut_events, Table(  waveform    = [cut_pulse],
                                        evtno       = [tt[i].evtno],
                                        detno       = [tt[i].detno],
                                        det_ch      = [tt[i].det_ch],
                                        edep        = [tt[i].edep]))
        end
    end
    if size(cut_events,1) > 0
        sp = zeros(length(cut_events[1].waveform.value))
        for i in 1:1:size(cut_events,1)
            sp .+= cut_events[i].waveform.value
        end
        sp ./= size(cut_events,1)
        sp ./= maximum(sp)
        data_sp = td[1].superpulse.value ./ maximum(td[1].superpulse.value)
        #chi2_a = sum((sp[1:400] .- data_sp[1:400]).^2 ./ data_sp[1:400] )
        #chi2_b = sum((sp[401:end] .- data_sp[401:end]).^2 ./ data_sp[401:end] )
        chi2 = sum(abs.((sp .- data_sp).^2 ./ data_sp) )
        #return abs(chi2_a) + abs(chi2_b), [GBP, tau, Cd, Cf]
        return chi2, [GBP, tau, Cd, Cf]
    else
        return false
    end
end


function check_parameters_AE(tt, td, GBP, tau, Cd, Cf, ch; std_AE=2)
    cut_events = Table( waveform    = [],
                        evtno       = [],
                        detno       = [],
                        det_ch      = [],
                        edep        = []);
    for i in 1:1:size(tt,1)
        pulse = tt[i].waveform.value
        filtered_pulse = applyElectronics(pulse; Ts = 10e-9, GBP = GBP, tau = tau, Kv = 150e3, Cd = Cd, Cf = Cf, Rf = 500e6)
        inter = find_intersect(filtered_pulse, maximum(filtered_pulse)/2, 5)
        if inter + 399 < length(pulse) && inter - 400 > 0
            cut  = filtered_pulse[(inter-399):1:(inter+400)]
            time = (0:10:7990)u"ns"
            cut_pulse = RDWaveform(time, cut)

            append!(cut_events, Table(  waveform    = [cut_pulse],
                                        evtno       = [tt[i].evtno],
                                        detno       = [tt[i].detno],
                                        det_ch      = [tt[i].det_ch],
                                        edep        = [tt[i].edep]))
        end
    end
    if size(cut_events,1) > 0
        sp = zeros(length(cut_events[1].waveform.value))
        if channel_to_bege[ch]
            #> With A/E filter
            A = []
            for i in 1:1:size(cut_events, 1)
                diff_pulse = diff(cut_events.waveform[i].value ./sum(cut_events.edep[i]).val) # normalized to energy
                a = maximum(diff_pulse)
                push!(A, a)
            end
            µ = mean(A)
            sig = std(A)
            good = []
            for i in 1:1:size(cut_events,1)
                if A[i] >= µ - std_AE*sig && A[i] <= µ + std_AE*sig
                    sp  .+= cut_events[i].waveform.value
                    push!(good, i)
                end
            end
            sp ./= length(good)
        else
            #> Normal superposition of all pulses
            for i in 1:1:size(cut_events,1)
                sp .+= cut_events[i].waveform.value
            end
            sp ./= size(cut_events,1)
        end
        sp ./= maximum(sp)
        data_sp = td[1].superpulse.value ./ maximum(td[1].superpulse.value)
        diff_data = diff(data_sp)
        diff_sim  = diff(sp)
        #chi2_a = sum((sp[1:400] .- data_sp[1:400]).^2 ./ data_sp[1:400] )
        #chi2_b = sum((sp[401:end] .- data_sp[401:end]).^2 ./ data_sp[401:end] )
        non_zero = findall(x->x != 0, diff_data)
        chi2 = sum( abs.((diff_sim[non_zero] .- diff_data[non_zero]).^2 ./ diff_data[non_zero]) )
        return chi2, [GBP, tau, Cd, Cf]
    else
        return false
    end
end


function plot_param_set(tt, td, params, ch; std_AE=2) 
    cut_events = Table( waveform    = [],
                        evtno       = [],
                        detno       = [],
                        det_ch      = [],
                        edep        = []);
    for i in 1:1:size(tt,1)
        pulse = tt[i].waveform.value
        filtered_pulse = applyElectronics(pulse; Ts = 10e-9, GBP = params[1], tau = params[2], Kv = 150e3, Cd = params[3], Cf = params[4], Rf = 500e6)
        inter = find_intersect(filtered_pulse, maximum(filtered_pulse)/2, 5)
        cut  = filtered_pulse[(inter-399):1:(inter+400)]
        time = (0:10:7990)u"ns"
        cut_pulse = RDWaveform(time, cut)

        append!(cut_events, Table(  waveform    = [cut_pulse],
                                    evtno       = [tt[i].evtno],
                                    detno       = [tt[i].detno],
                                    det_ch      = [tt[i].det_ch],
                                    edep        = [tt[i].edep]))
    end

    sp = zeros(length(cut_events[1].waveform.value))
    if channel_to_bege[ch]
        #> With A/E filter
        A = []
        for i in 1:1:size(cut_events, 1)
            diff_pulse = diff(cut_events.waveform[i].value ./sum(cut_events.edep[i]).val) # normalized to energy
            a = maximum(diff_pulse)
            push!(A, a)
        end
        µ = mean(A)
        sig = std(A)
        good = []
        for i in 1:1:size(cut_events,1)
            if A[i] >= µ - std_AE*sig && A[i] <= µ + std_AE*sig
                sp  .+= cut_events[i].waveform.value
                push!(good, i)
            end
        end
        sp ./= length(good)
    else
        #> Normal superposition of all pulses
        for i in 1:1:size(cut_events,1)
            sp .+= cut_events[i].waveform.value
        end
        sp ./= size(cut_events,1)
    end
    sp ./= maximum(sp)
    data_sp = td[1].superpulse.value ./ maximum(td[1].superpulse.value);
    diff_data = diff(data_sp);
    diff_sim  = diff(sp);

    chi2 = sum( abs.((diff_sim .- diff_data).^2 ./ diff_data) )
    p1 = plot(sp, label="Simulation: chi² of diff = "*string( round(chi2, digits=2) ));
    p1 = plot!(data_sp, color=:black, label="Superpulse (data)");
    p1 = plot!(title=channel_to_name[ch]*" | get electronics");
    p2 = plot(diff_data, label="Diff(data)");
    p2 = plot!(diff_sim, label="Diff(sim)");
    display(plot(p1,p2, layout=(2,1), xlim=[300, 500]))

end


function applyElectronics(pulse; Ts = 10e-9, GBP = 2750e6, tau = 180e-6, Kv = 150e3, Cd = 50e-12, Cf = 0.35e-12, Rf = 500e6)

    wop = GBP / (2 * pi * Kv)
    Cmod = Cf + Cd
    wmod = 1.0 / (Rf * Cmod)
    alfa = Cmod / (Cf * GBP)

    b0 = 1.0 / alfa
    a2 = 1.0
    a1 = 1.0 / alfa + wop + wmod
    a0 = 1.0 / (tau * alfa) + wmod*wop

    # then the transfer function in the *Laplace* s-domain looks like this:
    #                       b0
    #   T(s) = ----------------------------
    #             a2 * s^2 + a1 * s + a0

    # PolynomialRatio needs z-transform paramters: s- and z-domains can be connected by
    # the bilinear transform:
    #        z - 1
    # s = K -------- , K = Ts/2  , Ts - sampling period
    #        z + 1 
    #        
    # we can then convert T(s) to T(z):
    #              bz2 * z^2 + bz1 * z + bz0
    #   T(z) = -------------------------------
    #              az2 * z^2 + az1 * z + az0
    #

    K = 2/Ts

    az2 = 1.0   # normalized
    az1 = (2*a0 - 2*K^2)/(K^2 + a1*K + a0)
    az0 = (K^2 - a1*K + a0)/(K^2 + a1*K + a0)

    bz2 = b0/(K^2 + a1*K + a0)
    bz1 = 2*b0/(K^2 + a1*K + a0)
    bz0 = b0/(K^2 + a1*K + a0)

    myfilter = PolynomialRatio([bz2, bz1, bz0], [az2, az1, az0])

    filtered = filt(myfilter, vcat([0], diff(pulse)))

end

function add_baseline_and_extend_tail(wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}, n_baseline_samples::Int, total_waveform_length::Int) where {T,U,TV,UV}
    new_signal::Vector{eltype(UV)} = Vector{eltype(UV)}(undef, total_waveform_length)
    new_signal[1:n_baseline_samples] .= zero(eltype(UV))
    if length(wv.value) <= total_waveform_length - n_baseline_samples
        new_signal[n_baseline_samples+1:n_baseline_samples+length(wv.value)] = wv.value
        new_signal[n_baseline_samples+length(wv.value)+1:end] .= wv.value[end]
    else
        new_signal[n_baseline_samples+1:end] = wv.value[1:total_waveform_length - n_baseline_samples]
    end
    new_times = if TV <: AbstractRange
        range( zero(first(wv.time)), step = step(wv.time), length = total_waveform_length )
    else
        error("Not yet definted for timestamps of type `$(TV)`")
    end
    return RDWaveform( new_times, new_signal )
end