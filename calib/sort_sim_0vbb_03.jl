# -*- coding: utf-8 -*-
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")
include("fct-sim.jl")

plots_base_path = "pulses/sim/raw_0vbb/plots/"
base_path_raw   = "pulses/sim/raw_0vbb/"
base_path       = "pulses/sim/0vbb/"

parameters = JSON.parsefile("res/parameters.json")
bool_AE = false
std_intersect = 2;

# Search for available channels
channels = []
for c in 0:1:36 # change to 0:1:36 for ALL channels
    files = glob(base_path_raw*string(c)*"-"*channel_to_name[c]*"/*-wf*.h5");
    if length(files) > 0
        push!(channels, c)
    end
end

for ch in channels
    Base.run(`clear`)
#     IJulia.clear_output(true)
    str_ch = lpad(ch, 2, "0");
    @info("Processing Ch" * str_ch * " | "*channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
#     files = glob(base_path_raw*string(ch)*"-"*channel_to_name[ch]*"/"*channel_to_name[ch]*"-wf*.h5");
    files = glob(base_path_raw*string(ch)*"-"*channel_to_name[ch]*"/"*channel_to_name[ch]*"-cl-wf*.h5");
#     append!(files, files2)
#     return length(files)
    # prog = Progress(length(files), dt=0.5,
    #                         barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
    #                         barlen=10)

    for file in files
        file_nr = findall(x->x==file, files)[1]
        filename  = split(file, basename(file))[1]
        filename *= channel_to_name[ch] * "_w_filter" * split(basename(file), channel_to_name[ch])[2]
        if !isfile(filename)
            @info("Processing file "*string(file_nr)*"/"*string(length(files)) * " for Ch " * string(ch) * " | " * channel_to_name[ch])
            @info(">---------------------------------------------------------------<")
            data_raw = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "generated_waveforms")
            end

            t  = Table( waveform    = [],
                        evtno       = [],
                        event_mult  = [],
                        detno       = [],
                        det_ch      = [],
                        hits_totnum = [],
                        edep        = [],
                        pos         = []);
            temp = data_raw |> @filter(sum(_.waveform.value) > 0) |> Table
            
            if length(unique(data_raw.evtno)) != length(temp.evtno)
                evtno = unique(data_raw.evtno)
                @info("Sorting positive and negative waveforms")
                pro = Progress(length(evtno), dt=0.5,
                                barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                                barlen=10)
                for i in 1:1:length(evtno)
                    temp = data_raw |> @filter(_.evtno == evtno[i]) |> Table
                    if sum(temp[1].waveform.value) > 0
                        a = 1
                        b = 2
                    else 
                        a = 2
                        b = 1
                    end
                    wf     = add_baseline_and_extend_tail.(temp.waveform, 6000, 12000)
                    time   = (0:10:11990)u"ns"
                    value  = wf[a].value[1:10:end]
                    wf_cut = [RDWaveform(time, value)]

                    append!(t, Table(waveform   = wf_cut,
                                    evtno       = [evtno[i]],
                                    event_mult  = [temp[1].event_mult],
                                    detno       = [temp[1].detno],
                                    det_ch      = [ch],
                                    hits_totnum = [temp[1].hits_totnum],
                                    edep        = [temp[1].edep],
                                    pos         = [temp[1].pos]))
                    next!(pro)
                end
            else
                waveforms_cut = []
                wf     = add_baseline_and_extend_tail.(temp.waveform, 6000, 12000)
                time   = (0:10:11990)u"ns"
                for i in 1:1:size(temp,1)                    
                    value  = wf[i].value[1:10:end]
                    push!(waveforms_cut, RDWaveform(time, value))
                end
                t  = Table(waveform = waveforms_cut,
                        evtno       = temp.evtno,
                        event_mult  = temp.event_mult,
                        detno       = temp.detno,
                        det_ch      = zeros(length(temp.evtno)) .+ ch,
                        hits_totnum = temp.hits_totnum,
                        edep        = temp.edep,
                        pos         = temp.pos);
            end
                
            events = Table( waveform    = ArrayOfRDWaveforms( Array{typeof(t.waveform[1]), 1}(t.waveform) ),
                            evtno       = Int32.(t.evtno),
                            event_mult  = Int32.(t.event_mult),
                            detno       = VectorOfArrays(Array{Array{Int32,1},1}(t.detno)),
                            det_ch      = Int32.(t.det_ch),
                            hits_totnum = Int32.(t.hits_totnum),
                            edep        = VectorOfArrays(Array{typeof(t[1].edep),1}(t.edep)),
                            pos         = VectorOfArrays(Array{typeof(t[1].pos),1}(t.pos)));


            @info("Store intermediate file")
            filename  = split(file, basename(file))[1]
            filename *= channel_to_name[ch] * "_w_bl" * split(basename(file), channel_to_name[ch])[2]
            HDF5.h5open(filename, "w") do h5f
                LegendHDF5IO.writedata(h5f, "data", events)
            end 
            @info("Apply filter and time alignment")
            GBP = parameters[string(ch)]["filter"]["parameters"]["GBP"]*1e6
            tau = parameters[string(ch)]["filter"]["parameters"]["tau"]*1e-6
            Cd  = parameters[string(ch)]["filter"]["parameters"]["Cd"]*1e-12
            Cf  = parameters[string(ch)]["filter"]["parameters"]["Cf"]*1e-12
            
            if tau != 0 && GBP != 0
                data = Table( energy     = [],
                            multiplicity = [],
                            timestamp    = [],
                            run          = [],
                            channel      = [],
                            waveform     = [])
                filename  = split(file, basename(file))[1]
                filename *= channel_to_name[ch] * "_w_bl" * split(basename(file), channel_to_name[ch])[2]
                tt = HDF5.h5open(filename, "r") do h5f
                    LegendHDF5IO.readdata(h5f, "data")
                end
                pro = Progress(size(tt,1), dt=0.5,
                            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                            barlen=10)
                for i in 1:1:size(tt,1)
                    pulse = tt[i].waveform.value
                    if sum(pulse) > 0
                        filtered_pulse = applyElectronics(pulse; Ts = 10e-9, GBP = GBP, tau = tau, Kv = 150e3, Cd = Cd, Cf = Cf, Rf = 500e6)
                        inter = find_intersect(filtered_pulse, maximum(filtered_pulse)/2, 5)
                        cut  = filtered_pulse[(inter-399):1:(inter+400)]
                        time = (0:10:7990)u"ns"
                        cut_pulse = RDWaveform(time, cut)

                        append!(data, Table(energy       = [tt[i].edep],
                                            multiplicity = [tt[i].event_mult],
                                            timestamp    = [tt[i].evtno],
                                            run          = [0],
                                            channel      = [tt[i].det_ch],
                                            waveform     = [cut_pulse]))
                        next!(pro)
                    end
                end
                @info(length(data.waveform))
                @info("Storing waveforms...")
                events = Table(energy        = VectorOfArrays(Array{typeof(data[1].energy),1}(data.energy)),
                                multiplicity = Int32.(data.multiplicity),
                                timestamp    = Int32.(data.timestamp),
                                run          = Int32.(data.run),
                                channel      = Int32.(data.channel),
                                waveform     = ArrayOfRDWaveforms( Array{typeof(data[1].waveform), 1}(data.waveform) ))


                filename  = split(file, basename(file))[1]
                filename *= channel_to_name[ch] * "_w_filter" * split(basename(file), channel_to_name[ch])[2]
                HDF5.h5open(filename, "w") do h5f
                    LegendHDF5IO.writedata(h5f, "data", events)
                end
            end
        else
            @info("File "*string(file_nr)*"/"*string(length(files))*" has already been processed!")
        end
        Base.run(`clear`)
#         IJulia.clear_output(true)
        
        #next!(prog)
    end 
end
