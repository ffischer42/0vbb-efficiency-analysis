# -*- coding: utf-8 -*-
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")

sim_file        = "/remote/ceph2/group/gerda/data/simulation/gerda-mage-sim/gedet/intrinsic_bege/0nbb/edep/raw-gedet-intrinsic_bege-0nbb-edep-000.root"
config_dir      = "/res/Impurity_Scan/config-dep/"
plots_base_path = "pulses/sim/raw_0vbb/plots/"
base_path_raw   = "pulses/sim/raw_0vbb/"
base_path       = "pulses/sim/0vbb/"
n_sim_events    = 5000
sampling_time   = 1u"ns"
hits_threshold  = 0.005; # MeV
generate_wf    = true;
generate_cl_wf = true;

output_path     = base_path_raw*basename(sim_file)*"/";

if !isdir(output_path) 
    mkpath(output_path)
end

file = TFile(sim_file)
tree = file["fTree"];

step_lim = 1000000;
steps = []
i = [1]
while i[1] <= length(tree)
    if i[1] + step_lim -1 <= length(tree)
        push!(steps, i[1]:i[1] + step_lim - 1)
    else
        push!(steps, i[1]:length(tree))
    end
    i[1] += step_lim
end

for s in steps

#     IJulia.clear_output(true)
    Base.run(`clear`)
    @info(">------------------------------------------------------<")
    @info("File to be processed: "*basename(sim_file))
    @info("File size: "*string(round(stat(sim_file).size/1024/1024, digits=2))*" mb")
    @info(">------------------------------------------------------<")
    step_count = findall(x->x == s, steps)[1]
    step_count_str = lpad(step_count, 3, "0");
    
    tree = file["fTree"][s];


    tt = Table(eventnumber = tree.eventnumber[:],
        hits_iddet   = tree.hits_iddet[:],
        hits_totnum  = tree.hits_totnum[:],
        hits_edep    = tree.hits_edep[:],
        hits_xpos    = tree.hits_xpos[:] .* 10, # mm
        hits_ypos    = tree.hits_ypos[:] .* 10, # mm
        hits_zpos    = tree.hits_zpos[:] .* 10) # mm

    #########################################
    # Filter events without energy deposition
    tt = tt |> @filter(_.hits_totnum != 0) |> Table;
    tt = tt |> @filter(length(unique(_.hits_iddet)) == 1) |> Table;
#     tt = tt |> @filter(sum(_.hits_edep) > 0.700 && sum(_.hits_edep) < 1.300) |> Table

    i = 1
    i_max = size(tt, 1)


    evt_num      = []
    multiplicity = []
    iddet        = []
    totnum       = []
    edep         = []
    position     = []

    p = Progress(i_max, dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=10)
    @info("File " * step_count_str * " of " * string(length(steps)))
    @info("Number of events = "*string(i_max))
    @info("Separating events by detector id and filter by minimum energy deposition threshold = "*string(hits_threshold)*" MeV")
    #########################################
    # Separate events by detector ###########
    while i <= i_max
        multi = 0
        for det in unique(tt[i].hits_iddet)
            index = findall(x->x==det, tt[i].hits_iddet)
            if sum( tt[i].hits_edep[ index ] ) >= hits_threshold
                multi += 1
                upside_down = 1
                if parameters[sim_to_channel[det][2]]["upside_down"] == true
                    upside_down = -1
                end

                hits_iddet_mapped = []
                for in in index
                    push!(hits_iddet_mapped, sim_to_channel[det][1])
                end
                push!(iddet, Array{Int64,1}(hits_iddet_mapped))
                push!(edep,  tt[i].hits_edep[ index ])
                push!(position, [ SVector{3}(([ tt[i].hits_xpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_x"], 
                                                tt[i].hits_ypos[k] .- parameters[sim_to_channel[det][2]]["detcenter_y"], 
                                                upside_down .* (tt[i].hits_zpos[k] .- parameters[sim_to_channel[det][2]]["detcenter_z"] .+ upside_down * parameters[sim_to_channel[det][2]]["height"]/2) 
                                    ] * u"mm")...) for k in index ])
            end        
        end
        j = 1
        while j <= multi
            push!(evt_num, tt[i].eventnumber)
            push!(totnum, tt[i].hits_totnum)
            push!(multiplicity, multi)
            j += 1
        end
        next!(p)
        i += 1
    end
    tt = nothing
    #########################################
    # Create output table ###################
    filter_index = findall(x -> x == 1, multiplicity)

    data = Table(evtno = evt_num[filter_index], 
        event_mult       = multiplicity[filter_index], 
        detno            = iddet[filter_index],
        hits_totnum      = totnum[filter_index],
        edep             = VectorOfArrays(edep[filter_index] .*1000 *u"keV"),
        pos              = position[filter_index]
    )

    dets = []
    for detno in data.detno
        push!(dets, unique(detno)[1])
    end
    dets = unique(dets)
    p = Progress(length(dets), dt=0.5,
                    barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                    barlen=10)
    @info("Creating files for each detector")
    for detno in dets
        det = channel_to_name[detno]
        t = data |> @select(:evtno, :event_mult, :detno, :hits_totnum, :edep, :pos) |> @filter(unique(_.detno)[1] == detno) |> Table;
        file_path = output_path * string(detno) * "-" * det * "/"
        if !isdir(file_path)
            mkpath(file_path)
        end
        file_path *= step_count_str * ".h5"
        HDF5.h5open(file_path, "w") do h5f
            LegendHDF5IO.writedata( h5f, "data", Table(evtno = t.evtno, event_mult = t.event_mult, detno = VectorOfArrays(t.detno), hits_totnum = t.hits_totnum, 
                                                                    edep = VectorOfArrays(t.edep), pos = VectorOfArrays(t.pos)))
        end
        next!(p)
    end
end
