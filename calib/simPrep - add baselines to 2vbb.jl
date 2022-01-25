#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
# include("../src/worker_fct.jl")

calib_filepath = "../dicts/calib.json"
calib = JSON.parsefile(calib_filepath);

datasets_str = [
    "0053-0064",
    "0065-0079",
    "0080-0092",
    "0093-0113"
]
det_lib = JSON.parsefile("../dicts/det_lib.json");
bl_dict = JSON.parsefile("../dicts/bl_rms_limits.json");
parameters = JSON.parsefile("../dicts/electronics_parameters.json")

#
# Output paths and filter settings
plots_base_path = "../../waveforms/sim/plots/"
base_path_raw   = "../../../2020-02-06_8380701d_st_ffischer/pulses/sim/raw_2vbb/"
base_path       = "../../waveforms/sim/2vbb_check/"
plots_bl_path   = "../../waveforms/baselines/plots/"
base_bl_path    = "../../waveforms/baselines/wf/"

sampling_time   = 1u"ns";

channels = []
for ch in 0:1:36
    if channel_to_bege[ch] == false || ch in [5,6,7,13]
        continue
    end
    push!(channels, ch)
end
println(channels)
for ch in channels
    Base.run(`clear`)
    ds_str = det_lib[channel_to_name[ch]]["run_str"]

    str_ch = lpad(ch, 2, "0");
    @info("Start Ch$str_ch | " * channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
    @info("Read in baselines")
    bl_files = glob(joinpath(base_bl_path, "*/" * str_ch*"-"*channel_to_name[ch]*"/*.h5"))
    filtered = []
    for file in bl_files
        run = parse(Int64, split(split(basename(file), "-")[2], "run")[end])
        if run in calib[channel_to_name[ch]]["data"][ds_str]["used_runs"]
            push!(filtered, file)
        end
    end
    bl_files = filtered

    bl = Table(run = [], datasetID = [], waveform = [])
    @showprogress 1 "Filtering baselines for Ch$str_ch ..." for file in bl_files
        tmp_data = HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end
        tmp_data = tmp_data |> @filter(bl_dict[channel_to_name[ch]]["low_rms_limit"] < rms(_.waveform.value) < bl_dict[channel_to_name[ch]]["high_rms_limit"]) |> Table
        if size(tmp_data,1) > 0
            append!(bl, Table(run = tmp_data.run, datasetID = tmp_data.datasetID, waveform = tmp_data.waveform))
        end
    end

    bl = bl |> @filter(bl_dict[channel_to_name[ch]]["low_rms_limit"] < rms(_.waveform.value) < bl_dict[channel_to_name[ch]]["high_rms_limit"]) |> Table
    @info("Number of baselines: " * string(size(bl,1)))
    
    @info(">---------------------------------------------------------------<")
    files = glob(joinpath(base_path_raw, string(ch) * "-" * channel_to_name[ch] * "/*w_filter*.h5"));
    number_of_pulses = 0;
    key_str = []
    for file in files
        push!(key_str, split(basename(file), "_evts_")[1])
    end
    key_str = unique(key_str)
    for k in key_str
        last_file = basename(glob(joinpath(base_path_raw, string(ch) * "-" * channel_to_name[ch] * "/" * string(k) * "*.h5"))[end])
        number_of_pulses += parse(Int64, split(split(split(last_file, "_evts_")[end], "-")[end], ".h5")[1])
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
    @info("This results in " * string( size(bl, 1) * length(first_range) ) * " baselines for " * string(number_of_pulses) * " pulses")
    @info(">---------------------------------------------------------------<")
    @info("Indexing the baselines")
#     indices = Table(bl=[], first=[])
#     for i in 0:1:number_of_pulses
#         bl_id = i%size(bl,1) + 1
#         first = first_range[Int(ceil((i + 1) / size(bl,1)))]
#         append!(indices, Table(bl=[bl_id], first=[first]))
#     end
    indices = []
    for i in 0:1:(size(bl, 1) * length(first_range)-1)
        bl_id = i%size(bl,1) + 1
        first = first_range[Int(ceil((i+1) / size(bl,1)))]
        push!(indices, [bl_id, first])
    end
    indices = shuffle(indices)
    pulse_id = 1
    @info("Add baselines to pulses")
    
    @showprogress "Adding baselines..." for file in files
        create_file = false
        ch_str = lpad(ch, 2, "0");
        filename = joinpath(joinpath(base_path, str_ch * "-" * channel_to_name[ch]), basename(file))
        !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""
        if isfile(filename)
            if stat(filename).size <= stat(file).size/2
                create_file = true
            end
        else
            create_file = true
        end
        if create_file
            data_w_bl = Table(energy=[], multiplicity=[], timestamp=[], run=[], channel=[], waveform=[], evtno=[], pos=[])
            file_no_filter = joinpath(dirname(file), split(basename(file), "w_filter")[1] * "w_bl" * split(basename(file), "w_filter")[end])

            data = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end
            data_no_filter = HDF5.h5open(file_no_filter, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end

            for i in eachindex(data)
                event = data[i]
                time = event.waveform.time
                pulse = event.waveform.value ./ get_avg_maximum(event.waveform.value, 10)
                pulse .*= sum(event.energy).val #* 1000 # Energy is stored in MeV

                bl_id = indices[pulse_id][1]
                first = indices[pulse_id][2]

                last = first + length(pulse) - 1
                pulse .+= bl[bl_id].waveform.value[first:1:last]
    #             return plot(RDWaveform(time, pulse))
                append!(data_w_bl, Table(energy=[event.energy], 
                                        multiplicity=[event.multiplicity], 
                                        timestamp=[event.timestamp], 
                                        run=[event.run], 
                                        channel=[event.channel], 
                                        waveform=[RDWaveform(time, pulse)],
                                        evtno=[data_no_filter[i].evtno],
                                        pos=[data_no_filter[i].pos])
                    );
                pulse_id += 1
            end
            HDF5.h5open(filename, "w") do h5f
                LegendHDF5IO.writedata(h5f, "data", Table(energy     = VectorOfArrays(Array{typeof(data_w_bl[1].energy),1}(data_w_bl.energy)),
                                                        multiplicity = Int32.(data_w_bl.multiplicity),
                                                        timestamp    = Int32.(data_w_bl.timestamp),
                                                        run          = Int32.(data_w_bl.run),
                                                        channel      = Int32.(data_w_bl.channel),
                                                        waveform     = ArrayOfRDWaveforms( Array{typeof(data_w_bl[1].waveform), 1}(data_w_bl.waveform) ),
                                                        evtno        = Int32.(data_w_bl.evtno),
                                                        pos          = VectorOfArrays(Array{typeof(data_w_bl[1].pos), 1}(data_w_bl.pos)),
                        ))
            end
        else
            pulse_id += size(HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end, 1)
        end
    end
end
