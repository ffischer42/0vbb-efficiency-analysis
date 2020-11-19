using ArraysOfArrays
using CSV
using Dates
using Distributions
using DSP
using GERDAMetadata
using Glob
using HDF5
using Interpolations
using JSON
using LegendDataTypes
using LegendHDF5IO
using LegendHDF5IO: readdata, writedata
using UpROOT, LegendUpROOTIO
using LegendUpROOTIO: read_mgtevent, raw2mgtevent, mgdo2legend
using LsqFit
using Measurements
using Plots
using ProgressMeter
using Query
using RadiationDetectorSignals
using RadiationSpectra
using SolidStateDetectors
using StaticArrays
using StatsBase
using StatsPlots
using StructArrays
using Suppressor
using TypedTables
using Unitful

using Random, LinearAlgebra, Statistics
using BAT, IntervalSets
using ValueShapes
using Random123

SSD = SolidStateDetectors
T = Float32
@. gauss(x,p) = p[3]*(sqrt(2*pi*p[2]))^(-1) .* exp.(-(x .- p[1]).^2 ./ (2*p[2]^2));

@. model(x,p) = p[3]*(sqrt(2*pi*p[2]))^(-1) .* exp.(-(x .- p[1]).^2 ./ (2*p[2]^2))


function applyElectronics(pulse; Ts = 10e-9, GBP = 2750e6, tau = 180e-6, Kv = 150e3, Cd = 50e-12, Cf = 0.35e-12, Rf = 500e6)

    wop = GBP / (2 * pi * Kv)
    Cmod = Cf + Cd
    wmod = 1.0 / (Rf * Cmod)
    alfa = Cmod / (Cf * GBP)

    b0 = 1.0 / alfa
    a2 = 1.0
    a1 = 1.0 / alfa + wop + wmod
    a0 = 1.0 / (tau * alfa) + wmod*wop

    # then the transfer function in the *Laplace* s-domain looks like this:
    #                       b0
    #   T(s) = ----------------------------
    #             a2 * s^2 + a1 * s + a0

    # PolynomialRatio needs z-transform paramters: s- and z-domains can be connected by
    # the bilinear transform:
    #        z - 1
    # s = K -------- , K = Ts/2  , Ts - sampling period
    #        z + 1 
    #        
    # we can then convert T(s) to T(z):
    #              bz2 * z^2 + bz1 * z + bz0
    #   T(z) = -------------------------------
    #              az2 * z^2 + az1 * z + az0
    #

    K = 2/Ts

    az2 = 1.0   # normalized
    az1 = (2*a0 - 2*K^2)/(K^2 + a1*K + a0)
    az0 = (K^2 - a1*K + a0)/(K^2 + a1*K + a0)

    bz2 = b0/(K^2 + a1*K + a0)
    bz1 = 2*b0/(K^2 + a1*K + a0)
    bz0 = b0/(K^2 + a1*K + a0)

    myfilter = PolynomialRatio([bz2, bz1, bz0], [az2, az1, az0])

    filtered = filt(myfilter, vcat([0], diff(pulse)))

end

function add_baseline_and_extend_tail(wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}, n_baseline_samples::Int, total_waveform_length::Int) where {T,U,TV,UV}
    new_signal::Vector{eltype(UV)} = Vector{eltype(UV)}(undef, total_waveform_length)
    new_signal[1:n_baseline_samples] .= zero(eltype(UV))
    if length(wv.value) <= total_waveform_length - n_baseline_samples
        new_signal[n_baseline_samples+1:n_baseline_samples+length(wv.value)] = wv.value
        new_signal[n_baseline_samples+length(wv.value)+1:end] .= wv.value[end]
    else
        new_signal[n_baseline_samples+1:end] = wv.value[1:total_waveform_length - n_baseline_samples]
    end
    new_times = if TV <: AbstractRange
        range( zero(first(wv.time)), step = step(wv.time), length = total_waveform_length )
    else
        error("Not yet definted for timestamps of type `$(TV)`")
    end
    return RDWaveform( new_times, new_signal )
end