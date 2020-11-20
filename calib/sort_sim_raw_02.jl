cd("/remote/ceph/group/gedet/data/sim/2020/2020-02-06_8380701d_st_ffischer/")
simdir          = "/remote/ceph/group/gerda/data/simulation/gerda-mage-sim"

mapping_file    = simdir*"/UTILS/det-data/ged-mapping.json"
parameters_file = simdir*"/UTILS/det-data/ged-parameters.json"
config_dir      = "/res/Impurity_Scan/config-dep/"
n_sim_events    = 5000

windows = false
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")
include("fct-sim.jl")

plots_base_path = "pulses/sim/raw/plots/"
base_path_raw   = "pulses/sim/raw/"
base_path       = "pulses/sim/"

current_dir = pwd()
cd(simdir)
filenames = glob("calib/*/*/*/*.root")
cd(current_dir)

ch = 0;
isTP = 0;
isBL = 0;
mult = 1;
E = 300;
sampling_time   = 1u"ns"
hits_threshold  = 0.005; # MeV
generate_wf    = true;
generate_cl_wf = true;

# Search for available channels
channels = []
for c in 0:1:36
    files = glob("pulses/sim/raw/raw*/"*string(c)*"-"*channel_to_name[c]*".h5");
    if length(files) > 0
        push!(channels, c)
    end
end

for ch in reverse(channels)
    Base.run(`clear`)
    @info("Start ch"*string(ch)*" | "*channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
    files = glob(base_path_raw*"*/*-"*channel_to_name[ch]*".h5")

    for file in files   
        # Check whether the file has already been simulated
        @time events = h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end;
        last_file = glob(base_path_raw*"*"*channel_to_name[ch]*"/"*channel_to_name[ch]*"-wf-*"*split(file, "/")[end-1]*"*")
        if length(last_file) == 0
            existance_check = false
        else
            last_file       = last_file[end]
            last_event      = parse(Int64,split(split(last_file,"-")[end], ".h5")[1])
            existance_check = isfile(last_file) && last_event == size(events,1) && stat(last_file).size > 0
        end
        if existance_check
            @info("This file has already been simluated")
        else
            # Read the SSD simulation file of depleted detector
            config_file = glob("res/Impurity_Scan/config-dep/"*channel_to_name[ch]*"*.config")[1]
            sim_output  = glob("res/Impurity_Scan/output/"*channel_to_name[ch]*"*/" )[1]
            @info("Load SSD simulation")
            @time sim = @suppress readcreateh5(config_file, sim_output);
            # Open *.h5 file containing detector specific events



            isdir(base_path_raw*string(ch)*"-"*channel_to_name[ch]) ? "Directory exists" : mkpath(base_path_raw*string(ch)*"-"*channel_to_name[ch])
            file_id = findall(x->x == file, files)[1]
            output_dir       = base_path_raw*string(ch)*"-"*channel_to_name[ch]
            if generate_wf
                # IJulia.clear_output(true)
                Base.run(`clear`)
                @info("File $file_id | "*string(length(files)))
                @info("Start simulating waveforms for ch"*string(ch)*" | "*channel_to_name[ch])
                @info(">---------------------------------------------------------------<")
                output_basename  = channel_to_name[ch]*"-wf-"*split(file, "/")[end-1]
                t_with_waveforms = @suppress SSD.simulate_waveforms(events, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
            end

            if generate_cl_wf
                # IJulia.clear_output(true)
                Base.run(`clear`)
                @info("File $file_id | "*string(length(files)))
                @info("Start simulating clustered waveforms for ch"*string(ch)*" | "*channel_to_name[ch])
                @info(">---------------------------------------------------------------<")
                @info("$(sum(length.(events.edep))) hits before clustering")
                events_clustered = SSD.cluster_detector_hits(events, 0.2u"mm")
                @info("$(sum(length.(events_clustered.edep))) hits after clustering")
                output_basename            = channel_to_name[ch]*"-cl-wf-"*split(file, "/")[end-1]
                t_clustered_with_waveforms = @suppress SSD.simulate_waveforms(events_clustered, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
            end 
        end
    end 
end