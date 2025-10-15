######################################################
###                AGGREGATE FRAMES                ###
###            author: Michi Haugeneder            ###
######################################################
#=
Aggregate frames so that the resulting .nc-file is not too large
=#
using NCDatasets,PyCall, Statistics, LaTeXStrings, Dates #DelimitedFiles, Dates, Random, FFTW, Dates
#StatsBase,ImageFiltering,
import PyPlot
#GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
#animation = pyimport("matplotlib.animation")

importdir = joinpath(@__DIR__, "..")
pathtofile = "/home/haugened/Documents/data/sos/250917_for_Eli/"

include(joinpath(importdir, "src", "ir_evaluation.jl"))
#include(joinpath(importdir, "src", "general.jl"))
import .irev
#import .gen
PyPlot.pygui(true)

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################

ncfile = "230517_120045_frames14000_to_17999.nc"
outfile = joinpath(pathtofile, "230517_120045_frames14000_to_17999_avg300.nc")

######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###                   IMPORT STUFF                 ###
######################################################
(nrows, ncols, nframes) = irev.getdims(joinpath(pathtofile, ncfile), "irdata")

rowrange = collect(75:539)
colrange = collect(1:ncols)
framestoaverage = 30*10

nnewframes = div(nframes, framestoaverage)
agg_array = fill(NaN, length(rowrange), length(colrange), nnewframes)

for i in 1:nnewframes
    frameindices = (i - 1) * framestoaverage + 1:i * framestoaverage
    println("Processing frames ", frameindices[1], " to ", frameindices[end])
    IRdata = irev.loadexcerptfromNetCDF4(joinpath(pathtofile, ncfile), "irdata", rowrange, colrange, frameindices)
    agg_array[:, :, i] = mean(IRdata, dims=3)[:, :, 1]
end

irev.saveirasnetcdf(agg_array, "irdata", outfile)
