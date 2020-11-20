cd("/remote/ceph/group/gedet/data/sim/2020/2020-02-06_8380701d_st_ffischer/")
simdir          = "/remote/ceph/group/gerda/data/simulation/gerda-mage-sim"

mapping_file    = simdir*"/UTILS/det-data/ged-mapping.json"
parameters_file = simdir*"/UTILS/det-data/ged-parameters.json"
config_dir      = "/res/Impurity_Scan/config-dep/"
n_sim_events    = 5000

windows = false
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")
include("fct-sim.jl")

plots_base_path = "pulses/sim/raw/plots/"
base_path_raw   = "pulses/sim/raw/"
base_path       = "pulses/sim/"

current_dir = pwd()
cd(simdir)
filenames = glob("calib/*/*/*/*.root")
cd(current_dir)

ch = 0;
isTP = 0;
isBL = 0;
mult = 1;
E = 300;
sampling_time   = 1u"ns"
hits_threshold  = 0.005; # MeV
generate_wf    = true;
generate_cl_wf = true;



# Search for available channels
channels = []
for c in 5:1:36 # change to 0:1:36 for ALL channels
    files = glob("pulses/sim/raw/"*string(c)*"-"*channel_to_name[c]*"/*w_filter-*.h5");
    if length(files) > 0
        push!(channels, c)
    end
end

for ch in channels
    Base.run(`clear`)
    @info("Start ch"*string(ch)*" | "*channel_to_name[ch])
    @info(">---------------------------------------------------------------<")
    @info("Read in baselines")
    str_ch = lpad(ch, 2, "0");
    files  = glob("data_pulses/baselines_filtered/"*str_ch*"-"*channel_to_name[ch]*"/*.h5");

    bl = Table(energy=[], multiplicity=[], timestamp=[], channel=[], waveform=[])
    prog = Progress(length(files), dt=0.5,
                    barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                    barlen=10)
    for file in files
        append!(bl, HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end)
        next!(prog)
    end

    # get number of pulses
    
    files = glob("pulses/sim/raw/"*string(ch)*"-"*channel_to_name[ch]*"/*w_filter-*.h5");
    number_of_pulses = 0;
    key_str = []
    for file in files
        push!(key_str, split(basename(file), ".root_evts_")[1])
    end
    key_str = unique(key_str)
    for key in key_str
        last_file = basename(glob("pulses/sim/raw/"*string(ch)*"-"*channel_to_name[ch]*"/"*key*"*.h5")[end])
        number_of_pulses += parse(Int64, split(split(last_file, ".h5")[1], "-")[end])
    end
    @info("Find number of simulation pulses: "*string(number_of_pulses))
    @info("Find number of noise baselines: "*string(size(bl,1)))
    # get first_range
    if number_of_pulses <= size(bl,1)
        first_range = [1]
    else
        mult_factor = ceil(number_of_pulses/size(bl,1))
        first = 1
        last = 201
        step = Int(round((last-first)/(mult_factor)))
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
        filename  = "pulses/sim/"*ch_str*"-"*channel_to_name[ch]*"/"
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
        Base.run(`clear`)
        @info(string(pulse_id)*" pulses done")
        next!(prog)
    end 
end