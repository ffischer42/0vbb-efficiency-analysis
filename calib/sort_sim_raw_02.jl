function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this        = parse(Int64, this)

#
# Important directories
simdir          = "/remote/ceph2/group/gerda/data/mpik/gerda-simulations/gerda-gems-db/pss-dev"
simdir_old      = "/remote/ceph/group/gerda/data/simulation/gerda-mage-sim"
mapping_file    = simdir_old*"/UTILS/det-data/ged-mapping.json"
parameters_file = simdir_old*"/UTILS/det-data/ged-parameters.json"
config_dir      = "../../../2020-02-06_8380701d_st_ffischer/res/Impurity_Scan/config-dep/"

#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
include("../src/worker_fct.jl")

#
# Output paths and filter settings
plots_base_path = "../../waveforms/sim/plots/raw/"
base_path_raw   = "../../waveforms/sim/raw/"
base_path       = "../../waveforms/sim/raw_wf/"
n_sim_events    = 10000;
sampling_time   = 1u"ns";
generate_wf     = false;
generate_cl_wf  = true;
no_new_files    = 0

channels = []
for c in 0:1:36
    files = glob(base_path_raw * "raw*/" * lpad(c, 2, "0") * "-" * channel_to_name[c] * ".h5");
    if length(files) > 0
        push!(channels, c)
    end
end
channels = get_share_for_worker(channels, num_workers, this)

while no_new_files <= 3
    for ch in channels
        @info("Start Ch"*lpad(ch, 2, "0")*" | "*channel_to_name[ch])
        @info(">---------------------------------------------------------------<")
        files = glob(base_path_raw*"*/*-"*channel_to_name[ch]*".h5")

        # Read the SSD simulation file of depleted detector
        config_file = glob(config_dir * channel_to_name[ch] * "*.config")[1]
        sim_output  = glob(config_dir * "../output/"*channel_to_name[ch]*"*/" )[1]
        @info("Load SSD simulation")
        @time sim = @suppress readcreateh5(config_file, sim_output);

        for file in files
            output_dir = base_path * lpad(ch, 2, "0") * "-" * channel_to_name[ch] * "/cl-wf/"
            output_dir *= split(file, "/")[end-1] * "/"
            @time events = h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end;
            last_file = glob(output_dir * "*.h5")
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
                no_new_files = 0
                isdir(output_dir) ? "Directory exists" : output_dir
                file_id = findall(x->x == file, files)[1]
                
                if generate_wf
                    # IJulia.clear_output(true)
                    Base.run(`clear`)
                    @info("File $file_id | " * string(length(files)))
                    @info("Start simulating waveforms for ch" * lpad(ch, 2, "0") * " | " * channel_to_name[ch])
                    @info(">---------------------------------------------------------------<")
                    output_basename  = lpad(ch, 2, "0") * "-" * channel_to_name[ch]*"-wf-"*split(file, "/")[end-1]
                    t_with_waveforms = @suppress SSD.simulate_waveforms(events, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
                end

                if generate_cl_wf
                    # IJulia.clear_output(true)
                    Base.run(`clear`)
                    @info("File $file_id | "*string(length(files)))
                    @info("Start simulating clustered waveforms for ch" * lpad(ch, 2, "0") * " | " * channel_to_name[ch])
                    @info(">---------------------------------------------------------------<")
                    @info("$(sum(length.(events.edep))) hits before clustering")
                    events_clustered = SSD.cluster_detector_hits(events, 0.2u"mm")
                    @info("$(sum(length.(events_clustered.edep))) hits after clustering")
                    output_basename            = lpad(ch, 2, "0") * "-" * channel_to_name[ch] * "-cl-wf-"*split(file, "/")[end-1]
                    t_clustered_with_waveforms = SSD.simulate_waveforms(events_clustered, sim, output_dir, output_basename, chunk_n_physics_events = n_sim_events, Δt = sampling_time);
                end 
            end
        end
    end
    no_new_files += 1
    sleep(300)
end