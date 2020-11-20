function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this = parse(Int64, this)

println("Dataset: " * meta_key_file)
println("Worker " * string(this) * " of " * string(num_workers))
#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
include("../src/worker_fct.jl")
data_type = "cal" # phy
dataset = "v07.01"

if data_type == "phy"
    plots_base_path = "../../waveforms/data/plots/raw_" * dataset * "/"
    base_path_raw   = "../../waveforms/data/raw_" * dataset * "/"
    base_path       = "../../waveforms/data/" * dataset * "/"
elseif data_type == "cal"
    plots_base_path = "../../waveforms/calib/plots/raw_" * dataset * "/"
    base_path_raw   = "../../waveforms/calib/raw_" * dataset * "/"
    base_path       = "../../waveforms/calib/" * dataset * "/"
end
number_of_pulses = 1e4
bl_range = 1:1:200


files = glob(base_path_raw * "*/*.h5");
files = get_share_for_worker(files, num_workers, this)

data_cal = Dict()
for ch in 0:1:36
    data_cal[string(ch)] = Dict()
    data_cal[string(ch)]["data"] = Table(energy = [], run = [], channel = [], AoEvetoed = [], datasetID = [], AoEclassifier = [], A = [], E = [], waveform = [])
    data_cal[string(ch)]["counter"] = 1
end

pro = Progress(length(files), dt=0.5,
                barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                barlen=10);
for file in files
    if stat(file).size > 1000
        data = HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end;
        for ch in unique(data.channel)
            ch_str = lpad(ch, 2, "0");
            temp = data |> @filter(_.channel == ch) |> Table
#             cal_factor = cal[string(ch)]["value"]
            waveforms = []
            A = []
            E = []
            for wf in temp.waveform
                pulse   = -1 .*float.(wf.value)
                pulse .-= sum(pulse[bl_range])/length(bl_range)
#                 pulse .*= cal_factor
                push!(A, maximum(multi_mwa(diff(pulse),5,3)))
                push!(E, maximum(multi_mwa(diff(pulse),201,6)))
                push!(waveforms, RDWaveform(wf.time, pulse))
            end
            
            append!(data_cal[string(ch)]["data"], Table(energy   = temp.energy,
                                                        run      = temp.run,
                                                        channel  = temp.channel,
                                                        AoEvetoed= temp.AEvetoed,
                                                        datasetID= temp.datasetID,
                                                        AoEclassifier = temp.AEclassifier,
                                                        A        = A,
                                                        E        = E,
                                                        waveform = waveforms))
            if size(data_cal[string(ch)]["data"],1) >= number_of_pulses
                filename = base_path * channel_to_name[ch] * "-" * ch_str * "/"
                !isdir(filename) ? mkpath(filename) : "path exists"
                filename *= lpad(data_cal[string(ch)]["counter"], 4, "0") * ".h5"
                data_cal[string(ch)]["counter"] += 1
                HDF5.h5open(filename, "w") do h5f
                    LegendHDF5IO.writedata( h5f, "data", Table( energy       = float.(data_cal[string(ch)]["data"].energy),
                                                                run          = Int.(data_cal[string(ch)]["data"].run),
                                                                channel      = Int.(data_cal[string(ch)]["data"].channel),
                                                                AoEvetoed    = Int.(data_cal[string(ch)]["data"].AoEvetoed),
                                                                datasetID    = Int.(data_cal[string(ch)]["data"].datasetID),
                                                                AoEclassifier= float.(data_cal[string(ch)]["data"].AoEclassifier),
                                                                A            = float.(data_cal[string(ch)]["data"].A),
                                                                E            = float.(data_cal[string(ch)]["data"].E),
                                                                waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(data_cal[string(ch)]["data"].waveform))))
                end
                data_cal[string(ch)]["data"] = Table(energy = [], run = [], channel = [], AoEvetoed = [], datasetID = [], AoEclassifier = [], A = [], E = [], waveform = [])
            end
        end
    end
    next!(pro)
end
for ch in 0:1:36
    ch_str = lpad(ch, 2, "0");
    if size(data_cal[string(ch)]["data"],1) > 0
        filename = base_path * channel_to_name[ch] * "-" * ch_str * "/"
        !isdir(filename) ? mkpath(filename) : "path exists"
        filename *= lpad(data_cal[string(ch)]["counter"], 4, "0") * ".h5"
        data_cal[string(ch)]["counter"] += 1
        HDF5.h5open(filename, "w") do h5f
            LegendHDF5IO.writedata( h5f, "data", Table( energy       = float.(data_cal[string(ch)]["data"].energy),
                                                        run          = Int.(data_cal[string(ch)]["data"].run),
                                                        channel      = Int.(data_cal[string(ch)]["data"].channel),
                                                        AoEvetoed     = Int.(data_cal[string(ch)]["data"].AoEvetoed),
                                                        datasetID    = float.(data_cal[string(ch)]["data"].datasetID),
                                                        AoEclassifier= float.(data_cal[string(ch)]["data"].AoEclassifier),
                                                        A            = float.(data_cal[string(ch)]["data"].A),
                                                        E            = float.(data_cal[string(ch)]["data"].E),
                                                        waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(data_cal[string(ch)]["data"].waveform))))
        end
    end
end