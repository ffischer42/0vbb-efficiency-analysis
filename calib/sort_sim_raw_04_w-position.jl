# function main(args)
#     return args
# end
# # Load inline arguments for worker unit
# num_workers, this = main(ARGS)
# num_workers = parse(Int64, num_workers)
# this        = parse(Int64, this)

#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
include("../src/worker_fct.jl")

#
# Output paths and filter settings
plots_base_path = "../../waveforms/sim/plots/"
base_path_raw   = "../../waveforms/sim/raw_wf/"
base_path       = "../../waveforms/sim/wf_w_pos/"
plots_bl_path   = "../../waveforms/baselines/plots/"
base_bl_path    = "../../waveforms/baselines/wf/"
# base_bl_path_filter = "../../waveforms/baselines/wf_filter/"
sampling_time   = 1u"ns";
parameters = JSON.parsefile("../dicts/electronics_parameters.json")

channels = []
for c in 0:1:36
    files = glob(base_path_raw * lpad(c, 2, "0") * "-" * ctn[c] * "/w_filter_w_pos/*.h5");
    if length(files) > 0 && ctb[c]
        push!(channels, c)
    end
end
# channels = get_share_for_worker(channels, num_workers, this);
# @info(channels)

for ch in channels
#     IJulia.clear_output(true)
    Base.run(`clear`)
    str_ch = lpad(ch, 2, "0");
    
    
    @info("Start Ch$str_ch | " * ctn[ch])
    @info(">---------------------------------------------------------------<")
    @info("Read in baselines")
    bl_files  = glob(joinpath(base_bl_path, "*/" * str_ch*"-"*ctn[ch]*"/*.h5"));

    bl = Table(run = [], datasetID = [], waveform = [])
    @showprogress 1 "Collecting baselines for Ch$str_ch ..." for file in bl_files
        tmp_data = HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end
        if typeof(tmp_data.run) == Array{Int64,1}
            append!(bl, Table(run = tmp_data.run, datasetID = tmp_data.datasetID, waveform = tmp_data.waveform))
        else
            println(typeof(tmp_data.run))
        end
    end
    bl_dict = JSON.parsefile("../dicts/bl_rms_limits.json");
    bl = bl |> @filter(bl_dict[ctn[ch]]["low_rms_limit"] < rms(_.waveform.value) < bl_dict[ctn[ch]]["high_rms_limit"]) |> Table
    
    
    @info("Number of baselines: " * string(size(bl,1)))
    @info(">---------------------------------------------------------------<")
    files = glob(joinpath(base_path_raw, str_ch * "-" * ctn[ch] * "/w_filter_w_pos/*.h5"));
    number_of_pulses = 0;
    key_str = []
    for file in files
        push!(key_str, split(basename(file), ".root_evts_")[1])
    end
    key_str = unique(key_str)
    @showprogress 1 "Counting simulated waveforms..." for k in key_str
        last_file = basename(glob(joinpath(base_path_raw, str_ch * "-" * ctn[ch] * "/w_filter_w_pos/" * String(k) * "*.h5"))[end])
        number_of_pulses += parse(Int64, split(split(last_file, ".h5")[1], "-")[end])
    end
    
    @info("Found number of simulation pulses: " * string(number_of_pulses))
    @info("Found number of noise baselines: " * string(size(bl,1)))

    # get first_range
    if number_of_pulses <= size(bl,1)
        first_range = [1]
    else
        mult_factor = ceil(number_of_pulses/size(bl,1))
        first = 1
        last = 201
        step = Int(round((last-first)/(mult_factor)))
        first_range = first:step:last
    end
    
    
    @info("This results in "*string( size(bl, 1) *length(first_range) )*" baselines for "*string(number_of_pulses)*" pulses")
    @info(">---------------------------------------------------------------<")
    @info("Indexing the baselines")
    indices = Table(bl=[], first=[])
    for i in 0:1:number_of_pulses
        bl_id = i%size(bl,1) + 1
        first = first_range[Int(ceil((i + 1) / size(bl,1)))]
        append!(indices, Table(bl=[bl_id], first=[first]))
    end
    pulse_id = 1
    @info("Add baselines to pulses")

    
    @showprogress 1 "Adding baselines..." for file in files
        create_file = false
        ch_str = lpad(ch, 2, "0");
        filename = joinpath(
            joinpath(
                joinpath(base_path, 
                    str_ch * "-" * ctn[ch]), 
                split(basename(file), ".root")[1]),
            basename(file)
        )
        !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""

        if isfile(filename)
            if stat(filename).size <= stat(file).size/2
                create_file = true
            end
        else
            create_file = true
        end
        if create_file
            data_w_bl = Table(energy=[], multiplicity=[], timestamp=[], run=[], channel=[], waveform=[], pos=[])

            data = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end
            for event in data
                time = event.waveform.time
                pulse = event.waveform.value ./ get_avg_maximum(event.waveform.value, 10)
                pulse .*= sum(event.energy) * 1000 # Energy is stored in MeV

                bl_id = indices[pulse_id].bl
                first = indices[pulse_id].first

                last = first + length(pulse) -1
                pulse .+= bl[bl_id].waveform.value[first:1:last]

                append!(data_w_bl, Table(energy=[event.energy], 
                                        multiplicity=[event.multiplicity], 
                                        timestamp=[event.timestamp], 
                                        run=[event.run], 
                                        channel=[event.channel], 
                                        waveform=[RDWaveform(time, pulse)],
                                        pos = [event.pos]));
                pulse_id += 1
            end
            HDF5.h5open(filename, "w") do h5f
                LegendHDF5IO.writedata(h5f, "data", Table(energy     = VectorOfArrays(Array{typeof(data_w_bl[1].energy),1}(data_w_bl.energy)),
                                                        multiplicity = Int32.(data_w_bl.multiplicity),
                                                        timestamp    = Int32.(data_w_bl.timestamp),
                                                        run          = Int32.(data_w_bl.run),
                                                        channel      = Int32.(data_w_bl.channel),
                                                        waveform     = ArrayOfRDWaveforms( Array{typeof(data_w_bl[1].waveform), 1}(data_w_bl.waveform) ),
                                                        pos          = VectorOfArrays(Array{typeof(data_w_bl[1].pos),1}(data_w_bl.pos))))
            end
        else
            pulse_id += size(HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end, 1)

        end
    end
end
