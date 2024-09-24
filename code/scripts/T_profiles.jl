######################################################
###        VERTICAL PROFILES FROM SCREEN DATA      ###
###            author: Michi Haugeneder            ###
######################################################
#=
Currently only implemented for preprocessed HEFEX2-profiles
=#
using PyCall, Statistics, LaTeXStrings, Dates
#StatsBase,ImageFiltering,
import PyPlot
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
#animation = pyimport("matplotlib.animation")

importdir = joinpath(@__DIR__, "..")
pathtofile = "/media/michi/MHaugeneder_SSD_1/Documents/data/hefex2/profiles/10min_avg/"

include(joinpath(importdir, "src", "ir_evaluation.jl"))
include(joinpath(importdir, "src", "ir_rawdata_processing.jl"))
#include(joinpath(importdir, "src", "general.jl"))
import .irev
import .irraw
#import .gen
PyPlot.pygui(true)

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################
#which sequence should be loaded (string)
ncfile = "230823_1100_to_1110_ver1.nc"
targetfile = "230823_1100_to_1110_Tavg.nc"
readavg = true #true, if median and quartiles have already been written to file
avgfile = "230823_1100_to_1110_Tavg.nc"
#for plotting
heightlowestpoint = 0.2 #m
cmperpxlz = 0.586
#framespersec = 29.9
######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###                   IMPORT STUFF                 ###
######################################################
if !readavg
    profile = irraw.readprofilefromnetcdf(joinpath(pathtofile, ncfile))

    #temporal averaging profile
    prof_avg = zeros(Float64, size(profile, 1))
    prof_low_quart = zeros(Float64, size(profile, 1))
    prof_up_quart = zeros(Float64, size(profile, 1))

    for i in 1:size(profile, 1)
        (prof_low_quart[i], prof_avg[i], prof_up_quart[i]) = quantile(profile[i, :], [0.25, 0.5, 0.75])
    end

    #generate vector containing heights [m]
    z = collect(size(profile, 1)-1:-1:0) .* (cmperpxlz/100) .+ heightlowestpoint

    #save averaged profiles
    irraw.saveavgprofiletonetcdf(joinpath(pathtofile, targetfile), z, prof_avg, prof_low_quart, prof_up_quart)
else
    (z, prof_avg, prof_low_quart, prof_up_quart) = irraw.readavgprofilefromnetcdf(joinpath(pathtofile, avgfile))
end

######################################################
###                    PLOT DATA                   ###
######################################################
#plot temperature profiles
cmap = PyPlot.get_cmap("tab10");
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("IRscreen 23.08.2023 1100-1110 UTC, median and 25%-75%")
ax.plot(prof_avg, z, c=cmap(0))
ax.fill_betweenx(z, prof_low_quart, prof_up_quart, color=cmap(0), ec="face", alpha=0.5)
ax.set_xlabel(L"T_{air}~\mathrm{[^\circ C]}")
ax.set_ylabel(L"z~\mathrm{[m]}")
ax.grid()