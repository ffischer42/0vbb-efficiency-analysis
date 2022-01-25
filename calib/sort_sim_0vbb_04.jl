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
for c in 0:1:36
    files = glob(base_path_raw*string(c)*"-"*channel_to_name[c]*"/*w_filter-cl*.h5");
    if length(files) > 0
        push!(channels, c)
    end
end


for ch in channels
#     IJulia.clear_output(true)
    Base.run(`clear`)
    str_ch = lpad(ch, 2, "0");
    @info("Processing Ch" * str_ch * " | " * channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
    @info("Read in baselines")
    files  = glob("data_pulses/baselines_filtered/"*str_ch*"-"*channel_to_name[ch]*"/*.h5");

    bl = Table(energy=[], multiplicity=[], timestamp=[], channel=[], waveform=[])

    @showprogress for file in files
        append!(bl, HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end)
    end

    # get number of pulses
    
    files = glob(base_path_raw*string(ch)*"-"*channel_to_name[ch]*"/*w_filter-cl*.h5");

    number_of_pulses = 0;
    key_str = []
    for file in files
        push!(key_str, split(basename(file), "_evts_")[1])
    end
    key_str = unique(key_str)
    for key in key_str
        key_files = glob(base_path_raw*string(ch)*"-"*channel_to_name[ch]*"/"*key*"*.h5")
        last_file = split(basename(key_files[end]), ".h5")[1]
        skipped = 0
        for single_file in key_files
            if length(split(basename(single_file), "_skipped-")) > 1
                skipped += parse(Int64, split(split(basename(single_file), "_skipped-")[end], ".h5")[1])
#                 println(basename(single_file))
#                 return parse(Int64, split(split(split(last_file, "_evts_")[2], "_")[1], "-")[end])
            end
        end
        number_of_pulses += parse(Int64, split(split(split(last_file, "_evts_")[2], "_")[1], "-")[end]) #- skipped
    end
    #
    ## work-around for skipped events.. just increasing the number of indices
    number_of_pulses += 10000
    @info("Find number of simulation pulses: " * string(number_of_pulses))
    @info("Find number of noise baselines: " * string(size(bl,1)))
    # get first_range
    if number_of_pulses <= size(bl,1)
        first_range = [1]
    else
        mult_factor = ceil(number_of_pulses/size(bl,1))
        first = 1
        last = 201
        step = Int(floor((last-first)/(mult_factor)))
        first_range = first:step:last
    end
    
    
    @info("This results in "*string( size(bl, 1) *length(first_range) )*" baselines for "*string(number_of_pulses)*" pulses")
    @info(">---------------------------------------------------------------<")
    @info("Indexing the baselines")
    indices = Table(bl=[], first=[])
    for i in 0:1:number_of_pulses
        bl_id = i%size(bl,1) + 1
        first = first_range[Int(ceil((i + 1) / size(bl,1)))]
        append!(indices, Table(bl=[bl_id], first=[first]))
    end
    
    pulse_id = 1
    @info("Add baselines to pulses")

    prog = Progress(length(files), dt=0.5,
                    barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                    barlen=10)
    for file in files
        create_file = false
        ch_str = lpad(ch, 2, "0");
        filename  = base_path*ch_str*"-"*channel_to_name[ch]*"/"
        if !isdir(filename)
            mkpath(filename)
        end
        filename *= channel_to_name[ch] * "_w_noise" * split(basename(file), channel_to_name[ch])[2]
        if isfile(filename)
            if stat(filename).size <= stat(file).size/2
                create_file = true
            end
        else
            create_file = true
        end
        if create_file
            data_w_bl = Table(energy=[], multiplicity=[], timestamp=[], run=[], channel=[], waveform=[])

            data = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end

            for event in data
                time = event.waveform.time
                pulse = event.waveform.value ./ get_avg_maximum(event.waveform.value, 10)
                pulse .*= sum(event.energy).val

                bl_id = indices[pulse_id].bl
                first = indices[pulse_id].first

                last = first + length(pulse) -1
                pulse .+= bl[bl_id].waveform.value[first:1:last]


                append!(data_w_bl, Table(energy=[event.energy], 
                                        multiplicity=[event.multiplicity], 
                                        timestamp=[event.timestamp], 
                                        run=[event.run], 
                                        channel=[event.channel], 
                                        waveform=[RDWaveform(time, pulse)]));
                pulse_id += 1
            end

            HDF5.h5open(filename, "w") do h5f
                LegendHDF5IO.writedata(h5f, "data", Table(energy        = VectorOfArrays(Array{typeof(data_w_bl[1].energy),1}(data_w_bl.energy)),
                                                        multiplicity = Int32.(data_w_bl.multiplicity),
                                                        timestamp    = Int32.(data_w_bl.timestamp),
                                                        run          = Int32.(data_w_bl.run),
                                                        channel      = Int32.(data_w_bl.channel),
                                                        waveform     = ArrayOfRDWaveforms( Array{typeof(data_w_bl[1].waveform), 1}(data_w_bl.waveform) )))
            end
        else
            pulse_id += size(HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end, 1)

        end
#         IJulia.clear_output(true)
        Base.run(`clear`)
        @info(string(pulse_id)*" pulses done")
        next!(prog)
    end 
end
