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
data_type = "cal";
set = "sim";
high_cut = 4;

plots_path = "../plots/sim/"
base_path  = joinpath("../../waveforms/sim/", "wf")
base_path_AE  = joinpath("../../waveforms/sim/", "calib_AE")

calib_filepath = "../dicts/calib.json"
calib = JSON.parsefile(calib_filepath)
cut_lib_filepath = "../dicts/cut_lib.json"
cut_lib = JSON.parsefile(cut_lib_filepath)
AE_cal_filepath = "../dicts/AE_cal.json"
AE_cal = JSON.parsefile(AE_cal_filepath)
det_lib_filepath = "../dicts/det_lib.json"
det_lib = JSON.parsefile(det_lib_filepath)

channels = []
for ch in 0:1:36
    if ctb[ch] && ch != 6
        push!(channels, ch)
    end
end
for ch in channels
    ch_str = lpad(ch, 2, "0");
    if !channel_to_bege[ch]
        println("STOP here! This is a coax detector.")
    end
    if !haskey(calib, channel_to_name[ch])
        calib[channel_to_name[ch]] = Dict()
        calib[channel_to_name[ch]]["data"] = Dict()
        calib[channel_to_name[ch]]["sim"] = Dict()
    end

    filepath = joinpath(base_path_AE, ch_str * "-" * ctn[ch] * "-" * string(this) * ".h5")
    !isdir(dirname(filepath)) ? mkpath(dirname(filepath)) : ""
    if !isfile(filepath)
        numofele_A = 5
        BackDelta5 = div(numofele_A,2)
        ForwardDelta5 = isodd(numofele_A) ? div(numofele_A,2) : div(numofele_A,2) - 1
        numofele_E = 201
        BackDelta201 = div(numofele_E,2)
        ForwardDelta201 = isodd(numofele_E) ? div(numofele_E,2) : div(numofele_E,2) - 1

        files = glob(joinpath(base_path, ch_str * "-" * channel_to_name[ch] * "/*/*.h5"));
        files = get_share_for_worker(files, num_workers, this)
        data = Table(energy = [], A = [], E = [])
        @showprogress 1 "Read data and calculate A & E for Ch$ch_str" for file in files
            tmp = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end
            A = []
            E = []
            for wf in tmp.waveform
                push!(A, maximum(movingaverage(diff(wf.value),numofele_A,BackDelta5,ForwardDelta5,3)))
                push!(E, maximum(movingaverage(diff(wf.value),numofele_E,BackDelta201,ForwardDelta201,13)))
            end
            if size(tmp,1) > 1
                append!(data, Table(energy = sum.(tmp.energy)*1000, A = A, E = E))
            end
        end
        HDF5.h5open(filepath, "w") do h5f
            LegendHDF5IO.writedata(h5f, "data", Table(
                A = float.(data.A),
                E = float.(data.E),
                energy = float.(data.energy)
            ))
        end
        data = nothing
    end
end