ENV["GKSwstype"]="100"

function main(args)
    return args
end
# Load inline arguments for worker unit
num_workers, this = main(ARGS)
num_workers = parse(Int64, num_workers)
this = parse(Int64, this)
# num_workers, this = 10,1
println("Worker " * string(this) * " of " * string(num_workers))



include("../src/init.jl")
include("../src/fct.jl")
include("../src/fitting-fct.jl")
include("../src/worker_fct.jl")
data_type = "cal";
set = "sim";
high_cut = 4;

plots_path = "../plots/sim/"
base_path  = joinpath("../../waveforms/sim/", "wf")
base_path_AE  = joinpath("../../waveforms/sim/", "calib_AE")

calib_filepath = "../dicts/calib.json"
AE_cal_filepath = "../dicts/AE_cal.json"
cut_lib_filepath = "../dicts/cut_lib.json"

channels = []
for ch in 0:1:36
    if ctb[ch] && !(ch in [5,6,7])
        push!(channels, ch)
    end
end

channels = get_share_for_worker(channels, num_workers, this)
error_filepath = "error_" * string(this) * ".json"
error_log = !isfile(error_filepath) ? Dict() : JSON.parsefile(error_filepath)
open(error_filepath, "w") do f
    JSON.print(f, error_log, 4)
end

for ch in channels
    calib = JSON.parsefile(calib_filepath)
    cut_lib = JSON.parsefile(cut_lib_filepath)
    AE_cal = JSON.parsefile(AE_cal_filepath)
    ch_str = lpad(ch, 2, "0");
    if !ctb[ch]
        println("STOP here! This is a coax detector.")
    else
        if !haskey(calib, ctn[ch])
            calib[ctn[ch]] = Dict()
            calib[ctn[ch]]["data"] = Dict()
            calib[ctn[ch]]["sim"] = Dict()
        end

        calib = JSON.parsefile(calib_filepath)
        filepath = joinpath(base_path_AE, ch_str * "-" * ctn[ch] * "_AE_calibrated_smeared.h5")
        data = HDF5.h5open(filepath, "r") do h5f
            LegendHDF5IO.readdata(h5f, "data")
        end;
        plots_base_path = joinpath(plots_path, ch_str * "-" * ctn[ch] * "/" *set * "/")
        slice_lib_filepath = "../dicts/slice_lib.json"
        slice_filepath = "../dicts/slices_sim/"
        slice_lib = JSON.parsefile(slice_lib_filepath)
        slices = [600, 630, 660, 800, 895, 925, 955, 985, 1015, 1045,
                1130, 1160, 1190, 1220, 1250, 1280, 1310, 1340, 1370, 1400, 1430, 1460, 1545,
                1665, 1695, 1725, 1755, 1785, 1815, 1845, 1875, 1905, 1935, 
                1965, 1995, 2025, 2055,
                2145, 2175, 2205, 2235, 2265, 2295, 2325];
        overwrite = true
        for slice in eachindex(slices)
#         for slice in eachindex(slices)[1:1]
            slice_lib_filepath = joinpath(slice_filepath, ch_str * "/" * lpad(slices[slice], 4, "0") * ".json")
            if !isfile(slice_lib_filepath) || overwrite

                slice_lib = Dict()
                !haskey(slice_lib, ctn[ch]) ? slice_lib[ctn[ch]] = Dict() : ""
                !haskey(slice_lib[ctn[ch]], set) ? slice_lib[ctn[ch]][set] = Dict() : ""
                !haskey(slice_lib[ctn[ch]][set], "AE_slices") ? slice_lib[ctn[ch]][set]["AE_slices"] = Dict() : ""

                !isdir(dirname(slice_lib_filepath)) ? mkpath(dirname(slice_lib_filepath)) : ""
                open(slice_lib_filepath, "w") do f
                    JSON.print(f, slice_lib, 4)
                end
                slice_lib = JSON.parsefile(slice_lib_filepath)
                calib = JSON.parsefile(calib_filepath)
                @info("Start fitting slice " * string(slice) * " of " * string(length(slices)) * " | Detector: Ch" * ch_str * " - " * ctn[ch] * " (" * set * ")")

                E = deepcopy(data.E)
                A = deepcopy(data.A)
                AoE = A ./ E
                AoE ./= calib[ctn[ch]][set]["AE_norm"]
                index = findall(x -> x >= slices[slice] && x < slices[slice] + 30, E);
                rng_start = 0.7
                rng_end = 1.05
                rng_step = 0.001
                h = fit(Histogram, AoE[index], rng_start-0.1:rng_step:rng_end+0.1);

                if slice > 1
                    last_file = joinpath(slice_filepath, ch_str * "/" * lpad(slices[slice-1], 4, "0") * ".json")
                    tmp = JSON.parsefile(last_file)
                    use_past = haskey(tmp[ctn[ch]], set)
                    if use_past
                        use_past = haskey(tmp[ctn[ch]][set]["AE_slices"], string(slices[slice - 1]))
                    end
                else
                    use_past = false
                    last_file = joinpath(slice_filepath, ch_str * "/" * lpad(slices[slice], 4, "0") * ".json")
                end
                if slice == 1 || !use_past
                    prior = NamedTupleDist(
                        n = 20000..150000,
                        σ = 0.01..0.075,
                        μ = 0.98..1.02,
                        bkg_l = 0.1..3.0,
                        bkg_r = 10..500.0,
                        c = 50000..200000,
                        d = 100..10000.0,
                        β = 0.1..0.9
                    )
                else
                    slice_lib_last = JSON.parsefile(last_file)
                    par_dict = slice_lib_last[ctn[ch]][set]["AE_slices"][string(slices[slice - 1])]["peak"][1]
                    p0 = (
                        n = par_dict["n"],
                        σ = par_dict["σ"],
                        μ = par_dict["μ"],
                        bkg_l = par_dict["bkg_l"],
                        bkg_r = par_dict["bkg_r"],
                        c = par_dict["c"],
                        d = par_dict["d"],
                        β = par_dict["β"]
                    )
                    vary = 2
                    println(p0)
                    prior = NamedTupleDist(
                        n = p0.n/(2*vary)..p0.n*(2*vary),
                        σ = p0.σ/vary..p0.σ*vary,
                        μ = 0.98..1.02,
                        bkg_l = p0.bkg_l/vary..p0.bkg_l*vary,
                        bkg_r = p0.bkg_r/vary..p0.bkg_r*vary,
                        c = p0.c/(2*vary)..p0.c*vary,
                        d = p0.d/vary..p0.d*vary,
                        β = p0.β/vary..p0.β*vary
                    )
                end
                fail = false
                @info("Start fitting slice " * string(slice) * " of " * string(length(slices)) * " | Detector: Ch" * ch_str * " - " * ctn[ch] * " (" * set * ")")
                samples = try bat_fit(h, peak, prior; nsamples=5*10^4, nchains=4, params_first=false, norm_expected=false);
                catch
                    Base.run(`clear`)
                    @info("Start fitting slice " * string(slice) * " of " * string(length(slices)) * " | Detector: Ch" * ch_str * " - " * ctn[ch] * " (" * set * ")")
                    @info("Second try")
                    try bat_fit(h, peak, prior; nsamples=5*10^5, nchains=4, params_first=false, norm_expected=false);
                    catch
                        Base.run(`clear`)
                        @info("Start fitting slice " * string(slice) * " of " * string(length(slices)) * " | Detector: Ch" * ch_str * " - " * ctn[ch] * " (" * set * ")")
                        @info("Third try")
                        try bat_fit(h, peak, prior; nsamples=2*10^6, nchains=4, params_first=false, norm_expected=false);
                        catch
                            fail = true
                        end
                    end
                end
                p2 = plot();
                if !fail
                    par, err = mode(samples), std(samples)
                    x_fit = float.(midpoints(h.edges[1]))
                    y_fit = peak(x_fit, par)
                    y_fit ./= sum(y_fit)
                    y_fit .*= sum(h.weights)

                    par_tail = (n = 0, σ = par.σ, μ = par.µ, bkg_l = par.bkg_l, bkg_r = par.bkg_r, c = par.c, d = par.d, β = par.β)
                    y_tail = peak(x_fit, par_tail)
                    y_tail ./= sum(y_tail)
                    y_tail .*= (sum(peak(x_fit, par_tail)) / sum(peak(x_fit, par)))
                    y_tail .*= sum(h.weights)

                    par_gauss = (n = par.n, σ = par.σ, μ = par.µ, bkg_l = par.bkg_l, bkg_r = par.bkg_r, c = 0, d = par.d, β = par.β)
                    y_gauss = peak(x_fit, par_gauss)
                    y_gauss ./= sum(y_gauss)
                    y_gauss .*= (sum(peak(x_fit, par_gauss)) / sum(peak(x_fit, par)))
                    y_gauss .*= sum(h.weights)

                    !haskey(slice_lib[ctn[ch]][set]["AE_slices"], string(slices[slice])) ? slice_lib[ctn[ch]][set]["AE_slices"][string(slices[slice])] = Dict() : ""
                    slice_lib[ctn[ch]][set]["AE_slices"][string(slices[slice])]["peak"] = [par, err]
                    slice_filepath
                    open(slice_lib_filepath, "w") do f
                        JSON.print(f, slice_lib, 4)
                    end

                    E_range = string(slices[slice]) * "-" * string(slices[slice] + 30)
                    p2 = plot!(h, st=:step, label="Data");
                    p2 = plot!(x_fit, y_fit, label="Fit", lw=2);
                    p2 = plot!(x_fit, y_tail, label="");
                    p2 = plot!(x_fit, y_gauss, label="");
                    p2 = plot!(xlabel="Normalized A/E [a.u.]", ylabel="Samples", title=E_range * " keV | Ch$ch_str", framestyle=:box, legend=:topleft);
                    filename = joinpath(plots_base_path, "AE_slices/" * lpad(slices[slice], 4, "0") * "-peak_model.png")
                    !isdir(dirname(filename)) ? mkpath(dirname(filename)) : ""
                    savefig(p2, filename);
                    Base.run(`clear`)
                else
                    !haskey(error_log, ch_str) ? error_log[ch_str] = [] : ""
                    push!(error_log[ch_str], slice)
                    open(error_filepath, "w") do f
                        JSON.print(f, error_log, 4)
                    end
                end
            end
        end
    end
end