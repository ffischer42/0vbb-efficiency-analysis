include("src/init.jl")
include("src/fct.jl")
data_type = "cal"
dataset = "v07.01"

meta_key_file = "datasets/run0083-run0092-cal-analysis.txt";
meta_keys = CSV.read(meta_key_file);
channels = 0:1:36
event_step = Int(1e4)

current_dir = pwd()
ddir = "/remote/ceph/group/gerda/data/phase2/blind/" * dataset * "/gen/"
cd(ddir)
filenames1 = []
filenames4 = []
for meta_key in meta_keys[1]
    filename1 = ddir * glob("tier1/ged/" * data_type * "/" * split(meta_key, "-")[2] * "/" * meta_key * "*.root")[1]
    filename4 = ddir * glob("tier4/all/" * data_type * "/" * split(meta_key, "-")[2] * "/" * meta_key * "*.root")[1]
    push!(filenames1, filename1)
    push!(filenames4, filename4)
end
cd(current_dir)

if data_type == "phy"
    base_path = "pulses/data/raw_" * dataset * "/"
elseif data_type == "cal"
    base_path = "pulses/calib/raw_" * dataset * "/"
end
log_file = base_path * "log.json"
if !isdir(base_path)
    mkpath(base_path)
end
if !isfile(log_file)
    open(log_file, "w") do f
        JSON.print(f, Dict(), 4)
    end
    global log = Dict()
else
    global log = JSON.parsefile(log_file)
end;

for i in eachindex(filenames4)
    start_t = now()
    run = parse(Int64, split(split(filenames4[i], "gerda-run")[2], "-")[1])
    filename1 = filenames1[i]
    filename4 = filenames4[i]
    
    @info("Run " * string(run) * " | file " * string(i) * " of total " * string(length(filenames4)))
    log = JSON.parsefile(log_file)
    if length(keys(log)) >= 1
        estimation = 0
        for k in keys(log)
            estimation += log[k]
        end
        estimation /= length(keys(log))
        fileinfo = stat(filename4)
        est_t = round(fileinfo.size * estimation / 60, digits=1)
        @info("Estimated time for this file: " * string(est_t) * " min")
        total_filesize = 0
        for f in i:length(filenames4)
            total_filesize += stat(filenames4[f]).size
        end
        est_t = total_filesize * estimation / 60
        unit = " s"
        if est_t/60/60/24 > 1
            est_t /= 60*60*24
            unit = " d"
        elseif est_t/60/60 > 1
            est_t /= 60*60
            unit = " h"
        elseif est_t/60 > 1
            est_t /= 60
            unit = " min"
        end
        @info("Estimated time till completion: " * string(round(est_t, digits=1)) * unit)
    end
    @info("------------------------------")
    

    if !(filename4 in keys(JSON.parsefile(log_file)))
        @info("Load Tier4")
        @time tier4 = Table(TFile(filename4)["tier4"]);
        temp_E = sum.(tier4.energy[:])

        index_1 = findall(x->x == 1, tier4.multiplicity[:])
        index_2 = findall(y->y > 300, temp_E[index_1])
        index_3 = findall(x->x == 0, tier4.isBL[index_1[index_2]])
        index_4 = findall(x->x == 0, tier4.isTP[index_1[index_2[index_3]]])
        filtered_index = index_1[index_2[index_3[index_4]]]
        if length(filtered_index) > 0
            @info("Number of events: " * string(length(temp_E)))
            @info("Number of events after filter: " * string(length(filtered_index)))
            steps = []
            s = 1
            while s <= length(filtered_index)
                a = s
                b = s + event_step - 1 <= length(filtered_index) ? s+event_step-1 : length(filtered_index)
#                 println(string(filtered_index[a]) * " - " * string(filtered_index[b]))
                push!(steps, [a,b])
                s += event_step
            end
#             return
            for s in steps
                result = Table( energy       = [],
                        run          = [],
                        channel      = [],
                        AEvetoed     = [],
                        datasetID    = [],
                        AEclassifier = [],
                        waveform     = [])
#                 IJulia.clear_output(true)
                Base.run(`clear`)
                run_str = split(basename(filename4), "-")[2]
                @info("Run " * run_str * " | file " * string(i) * " of total " * string(length(filenames4)))
                @info("Step " * string(findfirst(x->x == s, steps)) * " of " * string(length(steps)))
                file = base_path * run_str * "/"
                if !isdir(file)
                    mkpath(file)
                end
                file *= basename(filename4) * lpad(findfirst(x->x == s, steps), 4, "0") * ".h5"
                if !isfile(file)
                    
                    tmp = filtered_index[s[1]:s[2]]
                    @info("Event " * string(tmp[1]) * " until " * string(tmp[end]))
                    tier4 = Table(TFile(filename4)["tier4"])[tmp];
                    @info("Load Tier1")
                    @time treeTier1 = TFile(filename1)["MGTree"].event;
                    @info("Apply first filter to Tier1 & Tier4")
                    @time waveforms = TypedTables.Table(raw2mgtevent.(treeTier1[tmp])).fAuxWaveforms
        #             waveforms = waveforms[tmp]
                    if length(tmp) > 0
                        run_str = split(basename(filename4), "-")[2]
                        run = parse(Int64, split(run_str, "run")[2])

                        pro = Progress(length(channels), dt=0.5,
                            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                            barlen=10);
                        for ch in channels

                            tmp = findall(y->y[ch + 1] != 0, tier4.energy)
                            if length(tmp) > 0
                                waveform = []
                                for i in tmp
                                    if data_type == "phy"
                                        push!(waveform, waveforms[i][ch + 1].wf)
                                    else
                                        push!(waveform, waveforms[i][1].wf)
                                    end
                                end

                                AEvetoed     = []
                                datasetID    = []
                                AEclassifier = []
                                for i in tmp
                                    push!(AEvetoed,     tier4[i].isAoEvetoed[ch+1])
                                    push!(datasetID,    tier4[i].datasetID[ch+1])
                                    push!(AEclassifier, tier4[i].AoEclassifier[ch+1])
                                end

                                append!(result, Table(  energy       = sum.(tier4.energy[tmp]),
                                                        run          = Int.(zeros(length(tmp)) .+ run),
                                                        channel      = Int.(zeros(length(tmp)) .+ ch),
                                                        AEvetoed     = AEvetoed,
                                                        datasetID    = datasetID,
                                                        AEclassifier = AEclassifier,
                                                        waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(waveform))))
                            end
                            next!(pro)
                        end
                    end
                    if size(result,1) > 0
                        file = base_path * run_str * "/"
                        file *= basename(filename4) * lpad(findfirst(x->x == s, steps), 4, "0") * ".h5"
                        HDF5.h5open(file, "w") do h5f
                            LegendHDF5IO.writedata( h5f, "data", Table( energy       = float.(result.energy),
                                                                        run          = Int.(result.run),
                                                                        channel      = Int.(result.channel),
                                                                        AEvetoed     = Int.(result.AEvetoed),
                                                                        datasetID    = Int.(result.datasetID),
                                                                        AEclassifier = float.(result.AEclassifier),
                                                                        waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(result.waveform))))
                        end
                    end
                end
            end
            log = JSON.parsefile(log_file)
            dt = now() - start_t
            fileinfo = stat(filenames4[1])
            fileinfo.size
            log[filename4] = (dt.value / 1000) / fileinfo.size
            open(log_file, "w") do f
                JSON.print(f, log, 4)
            end
        end
    end
#     IJulia.clear_output(true)
    Base.run(`clear`)
    log = JSON.parsefile(log_file)
    @info(string(length(log)) * " of " * string(length(filenames4)) * " done!")
end