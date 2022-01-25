#
## Needed for this script:
# - Energy calibration done
# - A/E calibration done

function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this = parse(Int64, this)
# num_workers, this = 10,1
println("Worker " * string(this) * " of " * string(num_workers))

#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
include("../src/worker_fct.jl")

#
## Choose "2vbb" or "0vbb"
bb = "0vbb"

if bb == "2vbb"
    base_path = "../../waveforms/sim/2vbb_check/"
    output_path = "../../waveforms/sim/2vbb_AoE/"
else
    base_path = "../../waveforms/sim/0vbb_check/"
    output_path = "../../waveforms/sim/0vbb_AoE/"
end
include("../src/fct.jl")
include("../src/fitting-fct.jl")
set = "sim"

calib_filepath = "../dicts/calib.json"
calib = JSON.parsefile(calib_filepath)
cut_lib_filepath = "../dicts/cut_lib.json"
cut_lib = JSON.parsefile(cut_lib_filepath)
AE_cal_filepath = "../dicts/AE_cal.json"
AE_cal = JSON.parsefile(AE_cal_filepath)
sf_lib_filepath = "../dicts/sf_lib.json"
sf_lib = JSON.parsefile(sf_lib_filepath);

channels = eval_ch

channels = get_share_for_worker(channels, num_workers, this)

for ch in channels
    ch_str = lpad(ch, 2, "0")

    Base.run(`clear`)
    @info("This worker has: " * string(channels))
    @info("---")
    @info("Start with Ch$ch_str")


    files = glob(joinpath(base_path, ch_str * "-" * ctn[ch] * "/*.h5"));
    
    numofele = 5
    BackDelta5 = div(numofele,2)
    ForwardDelta5 = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    numofele = 201
    BackDelta201 = div(numofele,2)
    ForwardDelta201 = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    
    data = Table(A=[], E=[], E_unsmeared=[], energy=[], evtno=[], pos=[])

    @showprogress "Calculating A & E " for file in files
        tmp = HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end

        A = []
        E = []
        for wf in tmp.waveform
            push!(A, maximum(movingaverage(diff(wf.value),5,BackDelta5,ForwardDelta5,3)))
            push!(E, maximum(movingaverage(diff(wf.value),201,BackDelta201,ForwardDelta201,13)))
        end
        E_cal = linmodel(E, calib[ctn[ch]][set]["lin_cal"][1])
        energy = []
        for e in sum.(tmp.energy)
            push!(energy, e.val)
        end
        append!(data, Table(A=A, E=E, E_unsmeared=E_cal, energy=energy, evtno=tmp.evtno, pos=tmp.pos))
    end

    smearing_dict = JSON.parsefile("../dicts/smearing.json")["Gauss"]
    p0 = [smearing_dict[string(ch)]["params"]["0"], smearing_dict[string(ch)]["params"]["1"], smearing_dict[string(ch)]["params"]["2"]]
    E = []
    @showprogress "Smearing energy " for e in data.E_unsmeared
        d = Normal(0, sqrt_fct(e,p0))
        push!(E, e + rand(d))
    end
    AoE = data.A ./ E
    AoE ./= calib[ctn[ch]][set]["AE_norm"]
    AoE ./= linmodel(E, AE_cal[ctn[ch]][set]["lin_fit"])
    AoE ./= AE_cal[ctn[ch]][set]["DEP_norm"]
    AoE .-= 1
    AoE ./= hypmodel(E, AE_cal[ctn[ch]][set]["sig_fit"])
    
    filename = joinpath(output_path, ch_str * "-" * ctn[ch] * "-AE_smeared.h5")
    !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""
    HDF5.h5open(filename, "w") do h5f
        LegendHDF5IO.writedata(h5f, "data", Table(A = Array{Float64,1}(data.A), 
                                                E = Array{Float64,1}(E),
                                                AoE = Array{Float64,1}(AoE),
                                                E_unsmeared = Array{Float64,1}(data.E_unsmeared),
                                                energy = Array{Float64,1}(data.energy),
                                                evtno = Array{Int64}(data.evtno),
                                                pos = VectorOfArrays(Array{typeof(data.pos[1])}(data.pos))
            )
        )
    end
end
