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

#
# Output paths and filter settings
plots_base_path = "../../../waveforms/baselines/plots/raw/"
base_path_raw   = "../../../waveforms/baselines/raw_wf/"
base_path       = "../../../waveforms/baselines/raw_wf/"

run_strs = glob("../../../waveforms/calib/raw_v07.01/*")
runs = []
for run in run_strs
    isdir(run) ? push!(runs, basename(run)) : ""
end

run_keys  = []
for run in runs
    current_dir = pwd()
    ddir = "/remote/ceph/group/gerda/data/phase2/blind/v07.01/gen"
    cd(ddir*"/tier4/all/phy/$run/")
    files = glob("*.root")
    for meta_key in files
        push!(run_keys, split(meta_key, "-phy-all-tier")[1])
    end 
    cd(current_dir)
end
@info("Run keys are gathered! Number of files: " * string(length(run_keys)))
run_keys = get_share_for_worker(run_keys, num_workers, this)
@info("This worker takes: " * string(length(run_keys)))

for id in run_keys
    Base.run(`clear`)
    # IJulia.clear_output(true)
    int_run = parse(Int64, split(split(id, "-")[2], "run")[end])
    str_run = lpad(int_run, 4, "0");
    ddir = "/remote/ceph/group/gerda/data/phase2/blind/v07.01/gen"
    channels = []
    for ch in 0:1:36
        filename = base_path_raw * "run" * str_run * "/" * lpad(ch, 2, "0") * "-" * channel_to_name[ch] * "/"
        !isdir(filename) ? mkpath(filename) : ""
        filename *= id * "-phy-bl-raw.h5"
        !isfile(filename) ? push!(channels, ch) : ""
    end
    println(length(channels))
    if length(channels) > 0
        filename1 = ddir*"/tier1/ged/phy/run"*str_run*"/"*id*"-phy-ged-tier1.root";
        filename4 = ddir*"/tier4/all/phy/run"*str_run*"/"*id*"-phy-all-tier4.root";
        tree = TFile(filename4)["tier4"];
        current_file_number = findall(x->x == id, run_keys)[1]
        @info("File "*string(current_file_number)*" of "*string(length(run_keys)))
        @info(">---------------------------------------------------------------------<")
        @info("Load waveforms")
        @time events = TypedTables.Table(raw2mgtevent.(TFile(filename1)["MGTree"].event[:]))
        @info("Start selecting & storing")
        waveforms    = events.fAuxWaveforms;
        for ch in channels
            index = findall(x->x[ch + 1] == 0, tree.energy[1:end])
            if length(index) > 0
                energy       = sum.(tree.energy[index])
                multiplicity = tree.multiplicity[index]
                timestamp    = tree.timestamp[index]
                channel      = Int.(zeros(length(index)) .+ ch)
                
                datasetID    = []
                map(x->push!(datasetID, x[ch+1]), tree.datasetID[index])
                
                waveform  = [] 
                map(x->push!(waveform, x[ch+1].wf), waveforms[index])

                filtered_dict = Table( 
                    energy       = energy,
                    multiplicity = multiplicity, 
                    timestamp    = timestamp,
                    run          = Int.(zeros(length(index)) .+ int_run),
                    datasetID    = datasetID,
                    channel      = channel, 
                    waveform     = waveform
                )
                events   = nothing;
                
                filename = base_path_raw * "run" * str_run * "/" * lpad(ch, 2, "0") * "-" * channel_to_name[ch] * "/"
                filename *= id * "-phy-bl-raw.h5"
                HDF5.h5open(filename, "w") do h5f
                    LegendHDF5IO.writedata( h5f, "data", Table( 
                            energy       = float.(filtered_dict.energy),
                            multiplicity = float.(filtered_dict.multiplicity), 
                            timestamp    = float.(filtered_dict.timestamp),
                            run          = Int.(filtered_dict.run),
                            datasetID    = Int.(filtered_dict.datasetID),
                            channel      = Int.(filtered_dict.channel), 
                            waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(filtered_dict.waveform))
                        )
                    )
                end

                @info("Ch" * lpad(ch, 2, "0") * "-" * channel_to_name[ch] * " | " * string(size(filtered_dict, 1)) * " baselines stored!")
            end
        end
    end
end