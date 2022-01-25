# -*- coding: utf-8 -*-
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")
include("fct-sim.jl")

sim_file        = "/remote/ceph2/group/gerda/data/simulation/gerda-mage-sim/gedet/intrinsic_bege/0nbb/edep/raw-gedet-intrinsic_bege-0nbb-edep-000.root"
config_dir      = "/res/Impurity_Scan/config-dep/"
plots_base_path = "pulses/sim/raw_0vbb/plots/"
base_path_raw   = "pulses/sim/raw_0vbb/"
base_path       = "pulses/sim/0vbb/"
n_sim_events    = 2500
sampling_time   = 1u"ns"
hits_threshold  = 0.005; # MeV
generate_wf    = false;
generate_cl_wf = true;



function simulate_waveforms_safe(events, output_filename, sim)
    wf_p = []
    wf_n = []
    skip = []
    no_skip = []
    @showprogress for evt in events
        edep =  evt.edep
        pos = []
        for sp in evt.pos
            push!(pos, T[sp[1].val, sp[2].val, sp[3].val] ./ 1000)
        end
        pos = CartesianPoint.(pos)
        e = Event(pos, edep) 
        try simulate!(e, sim)
            catch 
            push!(skip, findfirst(x->x == evt, events))
            continue
        end
        push!(no_skip, findfirst(x->x == evt, events))
        push!(wf_p, e.waveforms[1])
        push!(wf_n, e.waveforms[2])
    end
#     return events
    tt = Table(
        evtno = Int32.(vcat(events.evtno[no_skip], events.evtno[no_skip])),
        event_mult = Int32.(vcat(events.event_mult[no_skip], events.event_mult[no_skip])),
        detno = VectorOfArrays(vcat(events.detno[no_skip], events.detno[no_skip])),
        hits_totnum = Int32.(vcat(events.hits_totnum[no_skip], events.hits_totnum[no_skip])),
        edep = VectorOfArrays(vcat(events.edep[no_skip], events.edep[no_skip])),
        pos = VectorOfArrays(vcat(events.pos[no_skip], events.pos[no_skip])),
        waveform = ArrayOfRDWaveforms(Array{typeof(wf_p[1]), 1}(vcat(wf_p, wf_n)))
    );
    if length(skip) > 0
        output_filename = split(output_filename, ".h5")[1] * "_skipped-" * string(length(skip)) * ".h5"
    end
    HDF5.h5open(output_filename, "w") do h5f
        LegendHDF5IO.writedata(h5f, "generated_waveforms", tt)
    end
end


# Search for available channels
channels = []
for c in 6:1:12
    files = glob("pulses/sim/raw_0vbb/raw*/"*string(c)*"-"*channel_to_name[c]*"/");
    if length(files) > 0 && channel_to_bege[c]
        push!(channels, c)
    end
end

function check_charge_deposits(positions, sim)
    good_event = true
    for pos in positions
        point = CartesianPoint{Float32}(pos[1].val, pos[2].val, pos[3].val) ./1000
        if !in(point, sim.detector.semiconductors[1].geometry)
            good_event = false
        end
    end
    return good_event
end

for ch in channels
    Base.run(`clear`)
    @info("Start ch"*string(ch)*" | "*channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
    files = glob(base_path_raw*"*/*-"*channel_to_name[ch]*"/*.h5")
    p = Progress(length(files), dt=0.5,
                    barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                    barlen=10)
    for file in files
        # Check whether the file has already been simulated
        events = h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end
        
        last_file = glob(base_path_raw*"*"*channel_to_name[ch]*"/"*channel_to_name[ch]*"-cl-wf-*"*split(basename(file), ".")[1]*"_*")
        if length(last_file) == 0
            existance_check = false
        else
            last_file       = last_file[end]
            last_event      = parse(Int64,split(split(last_file, "_evts_")[2],"-")[end])
            existance_check = isfile(last_file) && last_event == size(events,1) && stat(last_file).size > 0
        end
        if existance_check == false
            # Read the SSD simulation file of depleted detector
            config_file = glob("res/Impurity_Scan/config-dep/"*channel_to_name[ch]*"*.config")[1]
            sim_output  = glob("res/Impurity_Scan/output/"*channel_to_name[ch]*"*/" )[1]
            # Open *.h5 file containing detector specific events
            sim = @suppress readcreateh5(config_file, sim_output);
            
            #
            ## Check for charge deposit location and filter if necessary
            t = events |> @filter(check_charge_deposits(_.pos, sim)) |> Table
            if size(events,1) != size(t,1)
                HDF5.h5open(file, "w") do h5f
                    LegendHDF5IO.writedata( h5f, "data", Table(evtno = t.evtno, event_mult = t.event_mult, detno = VectorOfArrays(t.detno), hits_totnum = t.hits_totnum, 
                                                                            edep = VectorOfArrays(t.edep), pos = VectorOfArrays(t.pos)))
                end
                events = h5open(file, "r") do h5f
                    LegendHDF5IO.readdata(h5f, "data")
                end
            end

            isdir(base_path_raw*string(ch)*"-"*channel_to_name[ch]) ? "Directory exists" : mkpath(base_path_raw*string(ch)*"-"*channel_to_name[ch])
            file_id = findall(x->x == file, files)[1]
            output_dir       = base_path_raw*string(ch)*"-"*channel_to_name[ch]
            if generate_wf
                output_basename  = channel_to_name[ch]*"-wf-"*split(basename(file), ".")[1]
                i = 1
                j = n_sim_events
                while i <= size(events,1)
                    if j > size(events,1)
                        j = size(events,1)
                    end
                    step_rng = i:1:j
                    str_rng = lpad(i, 6, "0") * "-" * lpad(j, 6, "0")
                    output_filename = joinpath(output_dir, output_basename) * "_evts_" * str_rng * ".h5"
                    @info(str_rng)
                    try @suppress SSD.simulate_waveforms(events[step_rng], sim, output_dir, output_basename * "_evts_" * str_rng, chunk_n_physics_events = 10000, Δt = sampling_time);
                    catch
                        @suppress_err simulate_waveforms_safe(events[step_rng], output_filename)
                    end

                    i += n_sim_events
                    j += n_sim_events
                end
#                 t_with_waveforms = @suppress SSD.simulate_waveforms(events, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
            end

            if generate_cl_wf
                events_clustered = SSD.cluster_detector_hits(events, 0.2u"mm")
                output_basename  = channel_to_name[ch]*"-cl-wf-"*split(basename(file), ".")[1]
                i = 1
                j = n_sim_events
                while i <= size(events,1)
                    if j > size(events,1)
                        j = size(events,1)
                    end
                    step_rng = i:1:j
                    str_rng = lpad(i, 6, "0") * "-" * lpad(j, 6, "0")
                    output_filename = joinpath(output_dir, output_basename) * "_evts_" * str_rng * ".h5"
                    
                    @info(str_rng)
                    try @suppress SSD.simulate_waveforms(events_clustered[step_rng], sim, output_dir, output_basename * "_evts_" * str_rng, chunk_n_physics_events = 10000, Δt = sampling_time);
                    catch
                        @info("Safe simulation")
                        @suppress_err simulate_waveforms_safe(events_clustered[step_rng], output_filename, sim)
                    end

                    i += n_sim_events
                    j += n_sim_events
                end

#                 t_clustered_with_waveforms = try @suppress SSD.simulate_waveforms(events_clustered, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
#                     catch; end
            end 
        end
        next!(p)
    end 
end
