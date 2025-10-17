######################################################
###                AGGREGATE FRAMES                ###
###            author: Michi Haugeneder            ###
######################################################
#=
Aggregate frames so that the resulting .nc-file is not too large
=#
using Statistics
#import PyPlot

importdir = joinpath(@__DIR__, "..")
pathtofile = "/home/haugened/Documents/data/sos/250917_for_Eli/230518_132940_frames0_to_54999"

include(joinpath(importdir, "src", "ir_evaluation.jl"))
#include(joinpath(importdir, "src", "general.jl"))
import .irev
#import .gen

rowrange = collect(121:599)
colrange = collect(1:1024)
framestoaverage = 30*10

#######################################################
#save aggregate for single file

#ncfile = "230517_120045_frames14000_to_17999.nc"
#outfile = joinpath(pathtofile, "230517_120045_frames14000_to_17999_avg300.nc")

(nrows, ncols, nframes) = irev.getdims(joinpath(pathtofile, ncfile), "irdata")

nnewframes = div(nframes, framestoaverage)
agg_array = fill(NaN, length(rowrange), length(colrange), nnewframes)

for i in 1:nnewframes
    frameindices = (i - 1) * framestoaverage + 1:i * framestoaverage
    println("Processing frames ", frameindices[1], " to ", frameindices[end])
    IRdata = irev.loadexcerptfromNetCDF4(joinpath(pathtofile, ncfile), "irdata", rowrange, colrange, frameindices)
    agg_array[:, :, i] = mean(IRdata, dims=3)[:, :, 1]
end

irev.saveirasnetcdf(agg_array, "irdata", outfile)

#######################################################
#aggregate for multiple nc files in a folder (not over file boundaries!)
ncfolder = pathtofile #joinpath(pathtofile, "nc")
ncs = readdir(ncfolder)
ncs = sort(ncs)
nfiles = length(ncs)
nnewframes_total = 0
for i in 1:nfiles
    (r, c, f) = irev.getdims(joinpath(ncfolder, ncs[i]), "irdata")
    nnewframes_total += div(f, framestoaverage)
end
agg_array_total = fill(NaN, length(rowrange), length(colrange), nnewframes_total)
currentframe = 1
for i in 1:nfiles
    (r, c, f) = irev.getdims(joinpath(ncfolder, ncs[i]), "irdata")
    nnewframes = div(f, framestoaverage)
    for j in 1:nnewframes
        frameindices = (j - 1) * framestoaverage + 1:j * framestoaverage
        println("Processing file ", ncs[i], ", frames ", frameindices[1], " to ", frameindices[end])
        IRdata = irev.loadexcerptfromNetCDF4(joinpath(ncfolder, ncs[i]), "irdata", rowrange, colrange, frameindices)
        agg_array_total[:, :, currentframe] = mean(IRdata, dims=3)[:, :, 1]
        currentframe += 1
    end
end
irev.saveirasnetcdf(agg_array_total, "irdata", joinpath(ncfolder, "all_sequences_agg.nc"))