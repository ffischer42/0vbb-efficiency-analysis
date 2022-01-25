@info("Load packages...")
windows = false
cd("/remote/ceph/group/gedet/data/sim/2020/2020-02-06_8380701d_st_ffischer")
include("gerda_data_processing/init.jl");
include("gerda_data_processing/fct.jl");


cal = JSON.parsefile("data_pulses/FEP_pulses/calibration/cal.json");

plots_base_path = "data_pulses/baselines_filtered/plots/RMS/"

pre_cut = 50        # keV
multiplicity = 1; 
std_intersect = 2;
number_of_pulses = 1e4


# Find all detectors in baselines_raw
folders = glob("data_pulses/baselines_raw/*/*")
channels = []
for folder in folders
    ch = parse(Int64, split(split(folder, "/")[end], "-")[1])
    ch in channels ? "Channel is already in the list" : push!(channels, ch)
end
# Search in baselines_raw for all files for detector X
for eval_channel in channels
    Base.run(`clear`)
    ch_str = lpad(eval_channel, 2, "0");
    raw_files = glob("data_pulses/baselines_raw/*/" * ch_str * "-" * channel_to_name[eval_channel] * "/*")
    @info("Starting with Ch"*ch_str*" | "*channel_to_name[eval_channel]*" including "*string(length(raw_files))*" files:")
    data = Table( energy     = [],
                multiplicity = [],
                timestamp    = [],
                channel      = [],
                waveform     = [])

    # Read in all files -> data
    pro = Progress(length(raw_files), dt=0.5,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
    for file in raw_files
        append!(data, HDF5.h5open(file, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end)
        next!(pro)
    end
    @info(string(size(data, 1))*" raw baselines found")
    @info("Start X-Talk filtering")
 
    # Energy calibration
    waveforms = []
    for wf in data.waveform
        pulse = wf.value .- sum(wf.value)/length(wf.value)
        pulse .*= cal[string(eval_channel)]["value"]
        push!(waveforms, RDWaveform(wf.time, pulse))
    end
    # Filter for X-Talk
    rms = []
    for wf in waveforms
        pulse = wf.value
        if maximum(abs.(pulse)) < pre_cut
            temp = sqrt(sum(pulse.^2) / length(pulse))
            push!(rms, temp)
        end
    end

    f = StatsBase.fit(Histogram, rms, closed=:left, minimum(rms):0.03:maximum(rms));
    x_hist = float.(f.edges[1][1:end-1]);
    y_hist = float.(f.weights);
    fit = curve_fit(model, x_hist, y_hist, [mean(rms), 1, maximum(y_hist)], lower=[minimum(x_hist), 0.01, maximum(y_hist)/2])
    x_fit = minimum(rms):0.01:maximum(rms)
    y_fit = model(x_fit, fit.param)
    max_x = findall(x -> x == maximum(y_fit), y_fit)[1]
    @info("Creating histogram and fit")
    p = plot(x_hist, y_hist, label="Baseline RMS");
    p = plot!(x_fit, y_fit, label="Fit");
    p = scatter!([x_fit[max_x]], [y_fit[max_x]], label="Maximum at "*string(round(x_fit[max_x], digits=2)))
    p = vline!([x_fit[max_x] - 2 * fit.param[2], x_fit[max_x] + 2 * fit.param[2]], label="2 sigma", legend = :outertopright)
    p = plot!(xlim=[x_fit[max_x] - 5 * fit.param[2], x_fit[max_x] + 5 * fit.param[2]])
    if !isdir(plots_base_path)
        mkpath(plots_base_path)
    end
    ch_str = lpad(eval_channel, 2, "0");
    @suppress savefig(p, plots_base_path*ch_str*"-"*channel_to_name[eval_channel]*".png");
    @suppress savefig(p, plots_base_path*ch_str*"-"*channel_to_name[eval_channel]*".pdf");


    bl_index = []
    for bl in 1:1:length(waveforms)
        wf = waveforms[bl]
        temp = sqrt(sum(wf.value.^2) / length(wf.value))
        if temp >= fit.param[1] - 2*fit.param[2] && temp <= fit.param[1] + 2*fit.param[2]
            push!(bl_index, bl)
        end
    end
    
    ch_str = lpad(eval_channel, 2, "0");
    filename = "data_pulses/baselines_filtered/"*ch_str*"-"*channel_to_name[eval_channel]*"/"
    if !isdir(filename)
        mkpath(filename)
    end
    bl_parts = []
    if length(bl_index) > number_of_pulses
        for i in 1:1:Int(ceil(length(bl_index)/number_of_pulses))
            if i*number_of_pulses <= length(bl_index)
                part = ((i-1)*number_of_pulses + 1):1:i*number_of_pulses
            else
                part = ((i-1)*number_of_pulses + 1):1:length(bl_index)
            end
            push!(bl_parts, bl_index[Int.(part)])
        end
    else
        push!(bl_parts, bl_index)
    end
    @info("Storing data")
    pro = Progress(length(bl_parts), dt=0.5,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
    for part in bl_parts
        part_num  = findall(x->x == part, bl_parts)[1]
        if part_num != length(bl_parts)
            init  = (part_num-1)*number_of_pulses +1
            final = part_num*number_of_pulses
        else
            init = (part_num-1)*number_of_pulses +1
            final = length(bl_index)
        end
        str_start = lpad( Int(init) , 7, "0")
        str_end   = lpad( Int(final), 7, "0")
        filename = "data_pulses/baselines_filtered/"*ch_str*"-"*channel_to_name[eval_channel]*"/"
        filename *= str_start*"-"*str_end*".h5"
        HDF5.h5open(filename, "w") do h5f
            LegendHDF5IO.writedata( h5f, "data", Table( energy       = float.(data.energy[part]),
                                                        multiplicity = float.(data.multiplicity[part]), 
                                                        timestamp    = float.(data.timestamp[part]),
                                                        channel      = float.(data.channel[part]), 
                                                        waveform     = StructArray{RDWaveform}(Array{RDWaveform,1}(waveforms[part]))
                                                        ))
        end
        next!(pro)
    end

end