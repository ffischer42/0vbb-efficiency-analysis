function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this = parse(Int64, this)
# num_workers, this = 10,1
println("Worker " * string(this) * " of " * string(num_workers))

include("../src/init.jl")
include("../src/fct.jl")
include("../src/fitting-fct.jl")
include("../src/worker_fct.jl")
data_type = "cal"
dataset = "v07.01"
datasets = [
    [53,64],
    [65,79],
    [80,92],
    [93,113]
]

channels = []
for ch in 0:1:36
    if channel_to_bege[ch] && ch != 6
        push!(channels, ch)
    end
end

channels = get_share_for_worker(channels, num_workers, this)
# channels = [0]
@info("")

for ch in channels
    ch_str = lpad(ch, 2, "0");
    plots_path = "../plots/data/"
    base_path  = "../../waveforms/calib/v07.01/"
    base_path_AE  = "../../waveforms/data/calib_AE/"

    calib_filepath = "../dicts/calib.json"
    calib = isfile(calib_filepath) ? JSON.parsefile(calib_filepath) : Dict()
    det_lib = JSON.parsefile("../dicts/det_lib.json")

    ch_str = lpad(ch, 2, "0");
    ds = det_lib[ctn[ch]]["ds"]

    high_cut = 4
    set = "data";
    if !haskey(calib, ctn[ch])
        calib[ctn[ch]] = Dict()
        calib[ctn[ch]]["data"] = Dict()
        calib[ctn[ch]]["sim"] = Dict()
    end
    
    
    filepath = joinpath(base_path_AE, ch_str * "-" * ctn[ch] * "-20MWA.h5")
    !isdir(dirname(filepath)) ? mkpath(dirname(filepath)) : ""
    if !isfile(filepath)

        files = glob(joinpath(base_path, ch_str * "-" * ctn[ch] * "/*.h5"));

        numofele = 5
        BackDelta5 = div(numofele,2)
        ForwardDelta5 = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
        numofele = 201
        BackDelta201 = div(numofele,2)
        ForwardDelta201 = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1

        all_data = Table(energy = [], run = [], channel = [], AoEvetoed = [], datasetID = [], AoEclassifier = [], A = [], E = [])
        @showprogress 1 "Read data for Ch" * ch_str * "..." for file in files
            tmp = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end
            tmp = tmp |> @filter(datasets[ds][1] <= _.run <= datasets[ds][2]) |> Table
            if channel_to_bege[ch] && ch != 6
                if 0 in unique(tmp.datasetID)
                    tmp = tmp |> @filter(_.datasetID == 0) |> Table
                    A = []
                    E = []
                    for wf in tmp.waveform
                        push!(A, maximum(movingaverage(diff(wf.value),5,BackDelta5,ForwardDelta5,3)))
                        push!(E, maximum(movingaverage(diff(wf.value),201,BackDelta201,ForwardDelta201,20)))
                    end
                    if size(tmp,1) > 1
                        append!(all_data, Table(energy = tmp.energy, run = tmp.run, channel = tmp.channel, AoEvetoed = tmp.AoEvetoed, datasetID = tmp.datasetID, AoEclassifier = tmp.AoEclassifier, A = A, E = E))
                    end
                end
            else
                A = []
                E = []
                for wf in tmp.waveform
                    push!(A, maximum(movingaverage(diff(wf.value),5,BackDelta5,ForwardDelta5,3)))
                    push!(E, maximum(movingaverage(diff(wf.value),201,BackDelta201,ForwardDelta201,20)))
                end
                if size(tmp,1) > 1
                    append!(all_data, Table(energy = tmp.energy, run = tmp.run, channel = tmp.channel, AoEvetoed = tmp.AoEvetoed, datasetID = tmp.datasetID, AoEclassifier = tmp.AoEclassifier, A = A, E = E))
                end
            end    
        end

        HDF5.h5open(filepath, "w") do h5f
            LegendHDF5IO.writedata(h5f, "data", Table(
                A = float.(all_data.A),
                E = float.(all_data.E),
                run = Array{Int64, 1}(all_data.run),
                AoEvetoed = Array{Int64, 1}(all_data.AoEvetoed),
                datasetID = Array{Int64, 1}(all_data.datasetID),
                AoEclassifier = float.(all_data.AoEclassifier),
                energy = float.(all_data.energy)
            ))
        end
    end
end