using Pkg
try using ArraysOfArrays 
catch 
    Pkg.add("ArraysOfArrays")
    using ArraysOfArrays 
end
try using CSV
catch
    Pkg.add("CSV")
    using CSV
end
try using Dates
catch
    Pkg.add("Dates")
    using Dates
end
try using Distributions
catch
    Pkg.add("Distributions")
    using Distributions
end
try using DSP
catch
    Pkg.add("DSP")
    using DSP
end
# using GERDAMetadata
try using Glob
catch
    Pkg.add("Glob")
    using Glob
end
try using HDF5
catch
    Pkg.add("HDF5")
    using HDF5
end
try using Interpolations
catch
    Pkg.add("Interpolations")
    using Interpolations
end
try using JSON
catch
    Pkg.add("JSON")
    using JSON
end
try using LaTeXStrings
catch
    Pkg.add("LaTeXStrings")
    using JSON
end
try using LegendDataTypes
catch
    Pkg.add(url="https://github.com/legend-exp/LegendDataTypes.jl")
    using LegendDataTypes
end
try using LegendHDF5IO
catch
    Pkg.add(url="https://github.com/legend-exp/LegendHDF5IO.jl")
    using LegendHDF5IO
end
using LegendHDF5IO: readdata, writedata
# try using UpROOT
# catch
#     Pkg.add("UpROOT")
#     using UpROOT
# end
# try using LegendUpROOTIO
# catch
#     Pkg.add(url="https://github.com/legend-exp/LegendUpROOTIO.jl")
#     using LegendUpROOTIO
# end
# using LegendUpROOTIO: read_mgtevent, raw2mgtevent, mgdo2legend
try using LsqFit
catch
    Pkg.add("LsqFit")
    using LsqFit
end
try using Measurements
catch
    Pkg.add("Measurements")
    using Measurements
end
try using Plots
catch
    Pkg.add("Plots")
    using Plots
end
try using ProgressMeter
catch
    Pkg.add("ProgressMeter")
    using ProgressMeter
end
try using Query
catch
    Pkg.add("Query")
    using Query
end
try using RadiationDetectorSignals
catch
    Pkg.add("RadiationDetectorSignals")
    using RadiationDetectorSignals
end
try using RadiationSpectra
catch
    Pkg.add("RadiationSpectra")
    using RadiationSpectra
end
try using SpecialFunctions
catch
    Pkg.add("SpecialFunctions")
    using SpecialFunctions
end
try using SolidStateDetectors
catch
    Pkg.add("SolidStateDetectors")
    using SolidStateDetectors
end
try using StaticArrays
catch
    Pkg.add("StaticArrays")
    using StaticArrays
end
try using StatsBase
catch
    Pkg.add("StatsBase")
    using StatsBase
end
try using StatsPlots
catch
    Pkg.add("StatsPlots")
    using StatsPlots
end
try using StructArrays
catch
    Pkg.add("StructArrays")
    using StructArrays
end
try using Suppressor
catch
    Pkg.add("Suppressor")
    using Suppressor
end
try using TypedTables
catch
    Pkg.add("TypedTables")
    using TypedTables
end
try using Unitful
catch
    Pkg.add("Unitful")
    using Unitful
end
try using Random
catch
    Pkg.add("Random")
    using Random
end
try using LinearAlgebra
catch
    Pkg.add("LinearAlgebra")
    using LinearAlgebra
end
try using Statistics
catch
    Pkg.add("Statistics")
    using Statistics
end
try using BAT
catch
    Pkg.add(url="https://github.com/bat/BAT.jl")
    using BAT
end
try using IntervalSets
catch
    Pkg.add("IntervalSets")
    using IntervalSets
end
try using ValueShapes
catch
    Pkg.add("ValueShapes")
    using ValueShapes
end
try using IJuliaBell
catch
    Pkg.add("IJuliaBell")
    using IJuliaBell
end

SSD = SolidStateDetectors
T = Float64
eval_ch = [0, 1, 2, 3, 4, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 30, 31, 32, 33, 34, 35]
