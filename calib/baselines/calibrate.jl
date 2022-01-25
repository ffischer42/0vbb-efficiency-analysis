function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this        = parse(Int64, this)

#
# Load packages and functions
include("../../src/init.jl")
include("../../src/fct.jl")
include("../../src/worker_fct.jl")
# using Distributed
#
# Output paths and filter settings
plots_base_path = "../../../waveforms/baselines/plots/raw/"
base_path_raw   = "../../../waveforms/baselines/raw_wf/"
base_path       = "../../../waveforms/baselines/wf/"
# addprocs(3)
# nprocs()
    
    
    
cal = JSON.parsefile("../../dicts/A_cal.json");
channels = 0:1:36
channels = get_share_for_worker(channels, num_workers, this)
for ch in channels
    files = glob(base_path_raw * "*/" * lpad(ch,2,"0") * "-" * channel_to_name[ch] * "/*.h5");
#     IJulia.clear_output(true)
    Base.run(`clear`)
    @info("Channel Ch" * lpad(ch, 2, "0") * " in progress")
    @showprogress for file in files
        filename = base_path * split(split(file, "/raw_wf/")[2], "-bl-raw.h5")[1] * ".h5"
        run_str = split(basename(file), "-")[2]
        run = parse(Int64, split(run_str, "run")[end])
        if haskey(cal, channel_to_name[ch])
            if haskey(cal[channel_to_name[ch]], run_str) && !isfile(filename)
                data = HDF5.h5open(file, "r") do h5f
                    LegendHDF5IO.readdata(h5f, "data")
                end
                cal_pulses = []
                map(x->push!(cal_pulses, RDWaveform(x.time, (float.(x.value) .- mean(x.value)) ./ cal[channel_to_name[ch]][run_str][1] )), data.waveform);
                
                !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""
                HDF5.h5open(filename, "w") do h5f
                    LegendHDF5IO.writedata(h5f, "data", Table(
                            energy = data.energy,
                            multiplicity = data.multiplicity,
                            timestamp = data.timestamp,
                            run = data.run,
                            datasetID = data.datasetID,
                            channel = data.channel,
                            waveform = StructArray{RDWaveform}(Array{RDWaveform,1}(cal_pulses))
                        )
                    )
                end
            end
        end
    end
end