# -*- coding: utf-8 -*-
include("gerda_data_processing/init.jl")
include("gerda_data_processing/fct.jl")
include("fct-sim.jl")

plots_base_path = "pulses/sim/raw_0vbb/plots/"
base_path_raw   = "pulses/sim/raw_0vbb/"
base_path       = "pulses/sim/0vbb/"
folders = glob(base_path*"*-*")
# cal = JSON.parsefile("data_pulses/FEP_pulses/calibration/cal.json");
fit_params = JSON.parsefile("pulses/fit_params.json")
Ecal = JSON.parsefile("pulses/cal_dict.json");

channels = []
# Find all detectors in data/raw/
for folder in folders
    ch = parse(Int64, split(split(folder, "/")[end], "-")[1])
    ch in channels ? "Channel is already in the list" : push!(channels, ch)
end


for ch in channels
    ch_str = lpad(ch, 2, "0");
#     IJulia.clear_output(true)
    Base.run(`clear`)
    @info("Start calculating A and E for Ch" * string(ch) * " | " * channel_to_name[ch])
    parameters = Table(A=[], E=[], energy=[])

    files = glob(base_path * ch_str * "-" * channel_to_name[ch] * "/*.h5");

    pro = Progress(length(files), dt=0.5,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)

    for file in files
        if stat(file).size > 1000
            data = HDF5.h5open(file, "r") do h5f
                LegendHDF5IO.readdata(h5f, "data")
            end
#             pulses = diff.(data.waveform.value);
#             A = maximum.(movingaverage.(pulses, 5, 3));
#             E = maximum.(movingaverage.(pulses, 201, 6));
            A = []
            E = []
            for wf in data.waveform
                push!(A, maximum(diff(multi_mwa(float.(wf.value),5,3))))
                push!(E, maximum(diff(multi_mwa(float.(wf.value),201,6))))
            end
            energy = []
            for e in sum.(data.energy)
                push!(energy, e.val)
            end
            append!(parameters, Table(A=A, E=E, energy=energy))
        end
        next!(pro)
    end

#     for e in parameters.energy
#         push!(energy, e.val)
#     end
    @info("Store file")
    filename = base_path * ch_str * "-" * channel_to_name[ch] * "-AE_Euncal.h5"
    HDF5.h5open(filename, "w") do h5f
        LegendHDF5IO.writedata(h5f, "data", Table(A    = Array{Float64,1}(parameters.A), 
                                                E_uncal= Array{Float64,1}(parameters.E), 
                                                energy = Array{Float64,1}(parameters.energy)))
    end
end
