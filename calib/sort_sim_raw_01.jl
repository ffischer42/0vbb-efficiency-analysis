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
config_dir      = "/res/Impurity_Scan/config-dep/"

#
# Load packages and functions
include("../src/init.jl")
include("../src/fct.jl")
include("../src/worker_fct.jl")

#
# Output paths and filter settings
plots_base_path = "../../waveforms/sim/plots/raw/"
base_path_raw   = "../../waveforms/sim/raw/"
base_path       = "../../waveforms/sim/processed/"
n_sim_events    = 5000
isTP = 0;
isBL = 0;
mult = 1;
sampling_time   = 1u"ns"
hits_threshold  = 0.005; # MeV
E = 400;


#
# Load workers log file
log_file = base_path_raw * "log-" * string(this) * ".json"
if !isdir(base_path_raw)
    mkpath(base_path_raw)
end
if !isfile(log_file)
    open(log_file, "w") do f
        JSON.print(f, Dict(), 4)
    end
    global log = Dict()
else
    global log = JSON.parsefile(log_file)
end;

#
# Collect simulation filepaths
current_dir = pwd()
cd(simdir)
filenames = glob("calib/*/*/*/*.root")
cd(current_dir)

filenames = get_share_for_worker(filenames, num_workers, this)

for filename in filenames
    E_lim = E
    output_path = base_path_raw*basename(filename)*"/";
    if !haskey(log, basename(filename))
        Base.run(`clear`)
        @info(string(findfirst(x->x == filename, filenames)) * " of " * string(length(filenames)) * " in progress!")
        file = TFile(joinpath(simdir, filename))
        tree = file["fTree"];

        tt = Table(eventnumber = tree.eventnumber[:],
            hits_iddet   = tree.hits_iddet[:],
            hits_edep    = tree.hits_edep[:],
            hits_xpos    = tree.hits_xpos[:] .* 10, # mm
            hits_ypos    = tree.hits_ypos[:] .* 10, # mm
            hits_zpos    = tree.hits_zpos[:] .* 10) # mm

        tt = tt |> @filter(length(_.hits_edep) != 0 && sum(_.hits_edep) >= E_lim/1000 && length(unique(_.hits_iddet)) == 1) |> Table;

        data = Table(
            evtno        = [],
            multiplicity = [],
            detno        = [],
            hits_totnum  = [],
            edep         = [],
            pos          = []
        )
        prog = Progress(size(tt,1), dt=0.5,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
        for i in eachindex(tt)
            multi = 1
            det = tt[i].hits_iddet[1]
            parameters[sim_to_channel[det][2]]["upside_down"] == true ? upside_down = -1 : upside_down = 1
            append!(data, Table(
            evtno        = [tt[i].eventnumber],
            multiplicity = [multi],
            detno        = [Array{Int64,1}(zeros(length(tt[i].hits_iddet)) .+ sim_to_channel[det][1])],
            hits_totnum  = [length(tt[i].hits_edep)],
            edep         = [tt[i].hits_edep],
            pos          = [ SVector{3}(([ tt[i].hits_xpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_x"], 
                                        tt[i].hits_ypos[k] .- parameters[sim_to_channel[det][2]]["detcenter_y"], 
                                        upside_down .* (tt[i].hits_zpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_z"] .+ upside_down * parameters[sim_to_channel[det][2]]["height"]/2) 
                            ] * u"mm")...) for k in eachindex(tt[i].hits_xpos) ]
                )
            )
            next!(prog)
        end
        tt = nothing
        dets = []
        for detno in data.detno
            push!(dets, unique(detno)[1])
        end
        dets = unique(dets)
        prog = Progress(length(dets), dt=0.5,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
        @info("Creating files for each detector")
        for detno in dets
            det = channel_to_name[detno]
            t = data |> @filter(unique(_.detno)[1] == detno) |> Table;
            !isdir(output_path) ? mkpath(output_path) : ""
            HDF5.h5open(output_path * lpad(detno, 2, "0") * "-" * det * ".h5", "w") do h5f
                LegendHDF5IO.writedata( h5f, "data", Table(
                    evtno       = t.evtno, 
                    event_mult  = t.multiplicity, 
                    detno       = VectorOfArrays(t.detno), 
                    hits_totnum = t.hits_totnum, 
                    edep        = VectorOfArrays(t.edep), 
                    pos         = VectorOfArrays(t.pos)))
            end
            next!(prog)
        end
        #
        # Log the progress
        log[basename(filename)] = "Done!"
        open(log_file, "w") do f
            JSON.print(f, log, 4)
        end
    end
end