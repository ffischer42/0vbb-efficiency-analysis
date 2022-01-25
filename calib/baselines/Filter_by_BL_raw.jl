windows = false
cd("/remote/ceph/group/gedet/data/sim/2020/2020-02-06_8380701d_st_ffischer")
include("gerda_data_processing/init.jl");
include("gerda_data_processing/fct.jl");

# All
ch_table = Table(ch = [0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36], 
                run = [86,86,86,86,86,78,78,67,66,65,86,86,86,86,86,87,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,83,86,86,86,86]);
#= # Just outliers
ch_table = Table(ch = [5, 6, 7, 7, 7, 13,32], 
                run = [78,78,67,66,65,87,83]);
# Just run0086
ch_table = Table(ch = [0, 1, 2, 3, 4, 8, 9, 10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,33,34,35,36], 
                run = [86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86]); =#
#= ch_table = Table(ch = [0], 
                run = [86]); =#
current_dir = pwd()
cal = JSON.parsefile("data_pulses/FEP_pulses/calibration/cal.json");

pre_cut = 50        # keV
multiplicity = 1; 
std_intersect = 2;


for s in unique(ch_table.run)
    
    set = findall(x->x==s, ch_table.run)[1]
    int_run = ch_table[set].run;
    str_run = lpad(int_run, 4, "0");

    ddir = "/remote/ceph/group/gerda/data/phase2/blind/v04.00/gen"
    cd(ddir*"/tier4/all/phy/run"*str_run*"/")
    files = glob("*.root")
    run_keys  = []
    for meta_key in files
        push!(run_keys, split(meta_key, "-phy-all-tier")[1])
    end 
    cd(current_dir)
    @info("Run keys for "* str_run *" are gathered! Number of files: "*string(length(run_keys)))
    #= pro = Progress(length(run_keys), dt=0.5,
                        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                        barlen=10) =#
    for id in run_keys
        filtered_dict = Dict();
        for ch in unique(ch_table.ch)
            filtered_dict[ch] = Table(  energy       = [],
                                        multiplicity = [], 
                                        timestamp    = [], 
                                        #det          = [],
                                        channel      = [], 
                                        waveform     = []);
        end
        channel_index   = findall(x->x==ch_table[set].run, ch_table.run);
        unique_channels = []
        new_id = true
        for eval_channel in unique(ch_table.ch[channel_index])
            ch_str = lpad(eval_channel, 2, "0");
            filename = "data_pulses/baselines_raw/run"*str_run*"/"*ch_str*"-"*channel_to_name[eval_channel]*"/"
            filename *= id*"-phy-bl-raw.h5"
            if !isfile(filename)
                push!(unique_channels, eval_channel)
            else
                @info("File already exist: "*filename)
            end
        end
        unique_channels = unique(unique_channels)



        filename1 = ddir*"/tier1/ged/phy/run"*str_run*"/"*id*"-phy-ged-tier1.root";
        filename4 = ddir*"/tier4/all/phy/run"*str_run*"/"*id*"-phy-all-tier4.root"

        #>------- Tier4 ---------------------------------------<#
        #>-----------------------------------------------------<#
        tree = TFile(filename4)["tier4"];

        Base.run(`clear`)
        current_file_number = findall(x->x == id, run_keys)[1]
        @info("File "*string(current_file_number)*" of "*string(length(run_keys)))
        @info(">---------------------------------------------------------------------<")

        

        for eval_channel in unique_channels
            index = [];
            #push!(index, findall(y->sum(y) >= E - ΔE && sum(y) <= E + ΔE, tree.energy[1:end]));     # Select energy range
            push!(index, findall(x->x[eval_channel + 1] == 0, tree.energy[1:end]));                  # Select channel
            #push!(index, findall(x->x <= 1, tree.multiplicity[1:end]));                             # Select multiplicity
            #push!(index, findall(x->x==0, tree.isBL[1:end]));                                       # Select "is baseline"
            #push!(index, findall(x->x==0, tree.isTP[1:end]));                                       # Select "no test pulse"
            # Add additional filters here
            # push!(index, findall(..., ...))

            filtered_index = index[1]
            #=
            for j in 2:1:length(index)-1
                for r in 1:1:length(filtered_index)
                    if !(filtered_index[r] in index[j])
                        filtered_index[r] = 0
                    end
                end
            end
            =#
            #filtered_index = filtered_index[findall(x->x != 0, filtered_index)];

            if length(filtered_index) > 0
                energy       = sum.(tree.energy[filtered_index])
                multiplicity = tree.multiplicity[filtered_index]
                timestamp    = tree.timestamp[filtered_index]
                eventNumber  = tree.eventNumber[filtered_index]
                #ch           = findall.(x-> x > 5, tree.energy[filtered_index]);
                eventChannelNumber = tree.eventChannelNumber[filtered_index]
                channel  = zeros(length(filtered_index)) .+ eval_channel
                #detector = []
                #for c in ch
                #    push!(channel, c[1]-1)
                #    push!(detector, channel_to_name[c[1]-1])
                #end
                events    = TypedTables.Table(raw2mgtevent.(TFile(filename1)["MGTree"].event[filtered_index]))
                waveforms = events.fAuxWaveforms;
                waveform  = [] 
                for w in waveforms
                    wf_time  = w[eval_channel+1].wf.time
                    #bl       = sum(w[eval_channel+1].wf.value) / length(w[eval_channel+1].wf.value)
                    wf_value = w[eval_channel+1].wf.value #.- bl
                    #wf_value .*= cal[string(eval_channel)]["value"]
                    push!(waveform, RDWaveform(wf_time, wf_value))
                end 
                events   = nothing;
                append!(filtered_dict[eval_channel], Table( energy       = energy,
                                                            multiplicity = multiplicity, 
                                                            timestamp    = timestamp, 
                                                            #det          = detector,
                                                            channel      = channel, 
                                                            waveform     = waveform));
            end
            @info("Ch"*string(eval_channel)*"-"*channel_to_name[eval_channel]*" | Events: "*string(length(filtered_dict[eval_channel].energy)))
            #= next!(pro) =#
        end

        

        @info(">---------------------------------------------------------------------<")
        @info("Start storing the raw baselines")

        for eval_channel in unique_channels
            if length(filtered_dict[eval_channel].energy) > 0
                
                #>------- Store channel -------------------------------<#
                #>-----------------------------------------------------<#
                ch_str = lpad(eval_channel, 2, "0");
                filename = "data_pulses/baselines_raw/run"*str_run*"/"*ch_str*"-"*channel_to_name[eval_channel]*"/"
                if !isdir(filename)
                    mkpath(filename)
                end

                #= str_start = lpad(1, 7, "0")
                str_end   = lpad(length(bl_index), 7, "0")
                filename *= str_start*"-"*str_end*".h5" =#
                bl_index = 1:1:length(filtered_dict[eval_channel].energy)
                filename *= id*"-phy-bl-raw.h5"
                HDF5.h5open(filename, "w") do h5f
                    LegendHDF5IO.writedata( h5f, "data", Table( energy       = float.(filtered_dict[eval_channel].energy[bl_index]),
                                                                multiplicity = float.(filtered_dict[eval_channel].multiplicity[bl_index]), 
                                                                timestamp    = float.(filtered_dict[eval_channel].timestamp[bl_index]),
                                                                #det          = string.(channel_data.det),
                                                                channel      = float.(filtered_dict[eval_channel].channel[bl_index]), 
                                                                waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(filtered_dict[eval_channel].waveform[bl_index])),
                                                                #wf_cut       = StructArray{RDWaveform}(Array{RDWaveform,1}(channel_data.wf_cut))
                                                                ))
                end
            
                @info("Ch"*string(eval_channel)*"-"*channel_to_name[eval_channel]*" | "*string(length(bl_index))*" baselines stored!")
            else
                @info("Ch"*string(eval_channel)*"-"*channel_to_name[eval_channel]*" | No baselines found!")
            end
        end
    end

end

include("Filter_by_BL_converter.jl")