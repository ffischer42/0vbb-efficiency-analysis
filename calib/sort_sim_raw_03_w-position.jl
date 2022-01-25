# -*- coding: utf-8 -*-
function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this        = parse(Int64, this)


#
##
# This script was only executed for Ch01 for demonstrations
##
#

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
base_path_raw   = "../../waveforms/sim/raw_wf/"
base_path       = "../../waveforms/sim/raw_wf_w-position/"
sampling_time   = 1u"ns";
no_new_files    = 0
parameters = JSON.parsefile("../dicts/electronics_parameters.json")

channels = []
for c in [1]#0:1:36
    files = glob(base_path_raw * lpad(c, 2, "0") * "-" * ctn[c] * "/cl-wf/raw*/*.h5");
    if length(files) > 0
        push!(channels, c)
    end
end
# channels = get_share_for_worker(channels, num_workers, this)

for ch in channels
    Base.run(`clear`)
    @info("Start Ch" * lpad(ch, 2, "0") * " | " * ctn[ch])
    @info(">---------------------------------------------------------------<")

    files = glob(base_path_raw * lpad(ch, 2, "0") * "-" * ctn[ch] * "/cl-wf/raw*/*.h5")
    files = get_share_for_worker(files, num_workers, this)
    
    pro = Progress(length(files), dt=0.5,
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=10)
    for file in files
        file_nr = findall(x->x==file, files)[1]
        filename  = dirname(file)
        filename *= ctn[ch] * "_w_filter" * split(basename(file), ctn[ch])[2]
        filename = split(filename, "/cl-wf/")[1] * "/w_filter_w_pos/" * split(filename, "/cl-wf/")[2]
        if !isfile(filename)
            data_raw = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "generated_waveforms")
            end
            first_half = 1:Int(size(data_raw,1) / 2) 
            waveform = add_baseline_and_extend_tail.(data_raw.waveform[first_half], 6000, 12000)

            #
            # Reduce time resolution to fit data
            waveform_w_bl = []
            map(x-> push!(waveform_w_bl, RDWaveform(x.time[1:10:end], x.value[1:10:end])), waveform)
            GBP = parameters[ctn[ch]]["par"]["GBP"]*1e6
            tau = parameters[ctn[ch]]["par"]["tau"]*1e-6
            Cd  = parameters[ctn[ch]]["par"]["Cd"]*1e-12
            Cf  = parameters[ctn[ch]]["par"]["Cf"]*1e-12

            data = Table( 
                energy       = [],
                multiplicity = [],
                timestamp    = [],
                run          = [],
                channel      = [],
                waveform     = [],
                pos          = []
            )
            for i in eachindex(waveform_w_bl)
                pulse = waveform_w_bl[i].value
                if sum(pulse) > 0
                    filtered_pulse = applyElectronics(pulse; Ts = 10e-9, GBP = GBP, tau = tau, Kv = 150e3, Cd = Cd, Cf = Cf, Rf = 500e6)
                    inter = find_intersect(filtered_pulse, maximum(filtered_pulse)/2, 5)
                    cut   = filtered_pulse[(inter-399):1:(inter+400)]
                    time  = (0:10:7990)u"ns"
                    cut_pulse = RDWaveform(time, cut)
                    append!(data, Table(
                            energy       = [data_raw[i].edep],
                            multiplicity = [data_raw[i].event_mult],
                            timestamp    = [data_raw[i].evtno],
                            run          = [0],
                            channel      = [unique(data_raw[i].detno)[1]],
                            waveform     = [cut_pulse],
                            pos          = [data_raw[i].pos]
                        )
                    )
                end


            end
            events = Table(
                energy       = VectorOfArrays(Array{typeof(data[1].energy),1}(data.energy)),
                multiplicity = Int32.(data.multiplicity),
                timestamp    = Int32.(data.timestamp),
                run          = Int32.(data.run),
                channel      = Int32.(data.channel),
                waveform     = ArrayOfRDWaveforms( Array{typeof(data[1].waveform), 1}(data.waveform) ),
                pos          =  VectorOfArrays(Array{typeof(data[1].pos),1}(data.pos)),
            )

            !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""
            HDF5.h5open(filename, "w") do h5f
                LegendHDF5IO.writedata(h5f, "data", events)
            end
        end
        next!(pro)
    end
end
