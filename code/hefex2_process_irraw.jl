######################################################
###      EXTRACTING PROFILES FROM IR-RAWDATA       ###
###        author: Michi Haugeneder                ###
######################################################
#=
TODO before executing this script:
 - Convert the raw .irb-files to single .txt-files using IRB2ASCII-software
 - Check for completeness of the extracted data (all frames there)

What this script does (execute line by line, because variables need to be adapted!!):
 - Extract vertical and horizontal profiles averaged for a given column (or row)
=#
using ProgressMeter, Statistics
using ReadWriteDlm2, Dates, NCDatasets #DelimitedFiles, 
import PyPlot

if gethostname() == "Michi-T450s" || gethostname() == "x1carbon5"
    importdir = "/home/michi/Documents/slf/ir/code/"
    pathtofile = "/media/michi/hefex2_haugeneder2/converted/"#"/home/michi/Documents/slf/data/"
    targetfile_stam = "/media/michi/MHaugeneder_SSD_1/Documents/data/hefex2/profiles/10min_avg/"#"/home/michi/Documents/slf/hefex2/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ir/code/"
    #pathtofile = "/home/haugened/Documents/data21_unordered/ir_todo"
    pathtofile = "/media/haugened/hefex2_haugeneder2/converted/230822_185513/0055_to_0155/"
    targetfile_stam =  "/home/haugened/Documents/data/hefex2/"
elseif gethostname() == "SLFW27953"
    importdir = "C:/Users/haugened/Documents/ibl_patch_snow/code/"
    pathtofile = "E:/converted/230822_095028/"
    targetfile_stam = "C:/Users/haugened/Documents/hefex2_converted/"
end
include(joinpath(importdir, "src", "ir_rawdata_processing.jl"))
include(joinpath(importdir, "src", "ir_evaluation.jl"))
import .irraw
import .irev

#include(joinpath(importdir, "src", "ir_evaluation.jl"))
#import .irev

"""
    sortnew(lengthin::Int64)::Vector{Int64}

Reorder already created files due to initially wrong sorting. Input length in frames,
output new sorting
"""
function sortnew(lengthin::Int64)::Vector{Int64}
    filesperseq = 1200
    nrtotalseq = div(lengthin, filesperseq)
    lastseq = mod(lengthin, filesperseq)
    intframesperseq = collect(1:filesperseq)
    stringframesperseq = Vector{String}(undef, length(intframesperseq))
    for i in 1:length(intframesperseq)
        if intframesperseq[i]<1000
            stringframesperseq[i] = lpad(intframesperseq[i], 3, "0")
        else
            stringframesperseq[i] = lpad(intframesperseq[i], 4, "0")
        end
    end
    stringsort_old = sort(stringframesperseq)
    sortperm_new = sortperm(parse.(Int64, stringsort_old))
    sortperm_out = zeros(Int64, lengthin)
    #print(stringsort_old[sortperm_new])
    for k in 1:nrtotalseq
        sortperm_out[(1:1200).+(k-1)*filesperseq] = sortperm_new .+ (k-1)*filesperseq
    end
    
    #treatment for last (not complete) sequence, if exists
    if lastseq != 0
        intframesperseq_last = collect(1:lastseq)
        stringframesperseq_last = Vector{String}(undef, length(intframesperseq_last))
        for i in 1:length(intframesperseq_last)
            if intframesperseq_last[i]<1000
                stringframesperseq_last[i] = lpad(intframesperseq_last[i], 3, "0")
            else
                stringframesperseq_last[i] = lpad(intframesperseq_last[i], 4, "0")
            end
        end
        stringsort_old_last = sort(stringframesperseq_last)
        #print(stringsort_old[sortperm_new])
        sortperm_new_last = sortperm(parse.(Int64, stringsort_old_last))
        sortperm_out[(1:lastseq).+(nrtotalseq)*filesperseq] = sortperm_new_last .+ (nrtotalseq)*filesperseq
    end
    return sortperm_out
end

#targetfile_hor1 =  joinpath(targetfile_stam, "230822_185513_0055_hor1.nc")
#targetfile_hor2 =  joinpath(targetfile_stam, "230822_185513_0055_hor2.nc")
targetfile_vert1 = joinpath(targetfile_stam, "230822_1900_to_1910_avg_verT.nc")

#tmpfile_hor1 = string(targetfile_hor1[1:end-3], "_tmp.nc")
#tmpfile_hor2 = string(targetfile_hor2[1:end-3], "_tmp.nc")
tmpfile_vert1 = string(targetfile_vert1[1:end-3], "_tmp.nc")

evaluationfolder = joinpath(pathtofile, "230822_1900_to_1910")#pathtofile
println()
println("-----------S-T-A-R-T-------------")

println("IRB2ASCII: Decimal seperator ., between values tab!!")

allfiles_tmp = readdir(evaluationfolder)
allfiles_tmp = allfiles_tmp[occursin.(r"irdata\d{2,3}_\d{3,4}\.txt", allfiles_tmp)]
numberoffiles = length(allfiles_tmp)
sort!(allfiles_tmp)
new_sort = sortnew(numberoffiles)
allfiles = allfiles_tmp[new_sort]

#visualize the hundred-fourtysecond frame of the first sequence
irtmp = irraw.fromIRloadIRpictureIRB2ASCII(joinpath(evaluationfolder, allfiles[142]), 0)
irev.showheatmaptrad(irtmp, 2, 15)

#set the coordinates of the profiles to extract (carful with x/y and row/column)
#hor_prof_left_row = [475, 660]
#hor_prof_left_col = [65, 65]
#hor_prof_right_row = hor_prof_left_row
#hor_prof_right_col = [510, 510]
#hor_prof_halfwidth = 25
vert_prof_top_row = [130]
vert_prof_col = [800]
vert_prof_bottom_row = [700]
vert_prof_halfwidth = 40

#calculate column ranges to extract
#hor_prof1_row_range = collect(-hor_prof_halfwidth:hor_prof_halfwidth).+hor_prof_left_row[1]
#hor_prof2_row_range = collect(-hor_prof_halfwidth:hor_prof_halfwidth).+hor_prof_left_row[2]
#hor_prof1_col_range = collect(hor_prof_left_col[1]:hor_prof_right_col[1])
#hor_prof2_col_range = collect(hor_prof_left_col[2]:hor_prof_right_col[2])
vert_prof_row_range = collect(vert_prof_top_row[1]:vert_prof_bottom_row[1])
vert_prof_col_range = collect(-vert_prof_halfwidth:vert_prof_halfwidth).+vert_prof_col[1]

#hor_prof1 = fill(NaN, size(hor_prof1_col_range, 1), numberoffiles)
#hor_prof2 = fill(NaN, size(hor_prof2_col_range, 1), numberoffiles)
vert_prof = fill(NaN, size(vert_prof_row_range, 1), numberoffiles)

#hor1_tmp_med = fill(NaN, size(hor_prof1_col_range, 1))
#hor2_tmp_med = fill(NaN, size(hor_prof2_col_range, 1))
vert_tmp_med = fill(NaN, size(vert_prof_row_range, 1))

temp_save_idcs = zeros(Int64, 9)
for ix in 1:length(temp_save_idcs)
    temp_save_idcs[ix] = round(Int, ix/10*numberoffiles)
end

@showprogress "extracting profiles" for i in 1:numberoffiles
    a = readdlm2(joinpath(evaluationfolder, allfiles[i]), '\t'; decimal='.', skipstart=1, use_mmap=true)
    arr = Array{Float64}(a[:, 1:end-1])

    #extract profiles
    #hor1_tmp = arr[hor_prof1_row_range, hor_prof1_col_range]
    #hor2_tmp = arr[hor_prof2_row_range, hor_prof2_col_range]
    vert_tmp = arr[vert_prof_row_range, vert_prof_col_range]

    #take medians
    #=for j in 1:size(hor1_tmp, 2)
        hor1_tmp_med[j] = median(hor1_tmp[:,j])
        hor2_tmp_med[j] = median(hor2_tmp[:,j])
    end=#
    for k in 1:size(vert_tmp, 1)
        vert_tmp_med[k] = median(vert_tmp[k,:])
    end

    #hor_prof1[:, i] = hor1_tmp_med
    #hor_prof2[:, i] = hor2_tmp_med
    vert_prof[:, i] = vert_tmp_med

    #store temporary file every 10% in case script crashes
    if i in temp_save_idcs
    #    saveprofiletonetcdf(tmpfile_hor1, hor_prof1)
    #    saveprofiletonetcdf(tmpfile_hor2, hor_prof2)
        irraw.saveprofiletonetcdf(tmpfile_vert1, vert_prof)        
    end
end

#saveprofiletonetcdf(targetfile_hor1, hor_prof1)
#saveprofiletonetcdf(targetfile_hor2, hor_prof2)
saveprofiletonetcdf(targetfile_vert1, vert_prof)

#removing temporary files
#rm(tmpfile_hor1, force=true)
#rm(tmpfile_hor2, force=true)
rm(tmpfile_vert1, force=true)

println("-----------D-O-N-E---------------")
##########################################################
#=
hor_prof1_old = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/old_wrong_sorting/230822_185513_2325_hor1.nc")
hor_prof2_old = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/old_wrong_sorting/230822_185513_2325_hor2.nc")
vert_prof_old = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/old_wrong_sorting/230822_185513_2325_ver1.nc")

hor_prof1_sortperm = sortnew(size(hor_prof1_old, 2))
#hor_prof2_sortperm = sortnew(size(hor_prof2_old, 2))
#vert_prof_sortperm = sortnew(size(vert_prof_old, 2))

hor_prof1_new = hor_prof1_old[:, hor_prof1_sortperm]
hor_prof2_new = hor_prof2_old[:, hor_prof1_sortperm]
vert_prof_new = vert_prof_old[:, hor_prof1_sortperm]

saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/new_sorting/a230822_185513_2325_hor1.nc", hor_prof1_new)
saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/new_sorting/a230822_185513_2325_hor2.nc", hor_prof2_new)
saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/new_sorting/a230822_185513_2325_ver1.nc", vert_prof_new)
=#

#=
#concentating two profiles
hor_prof1_1 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_to_463_hor1.nc")
hor_prof2_1 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_to_463_hor2.nc")
vert_prof_1 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_to_463_ver1.nc")

hor_prof1_2 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_464_to_end_hor1.nc")
hor_prof2_2 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_464_to_end_hor2.nc")
vert_prof_2 = readprofilefromnetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_464_to_end_ver1.nc")

hor_prof1_conc = cat(hor_prof1_1, hor_prof1_2, dims=2)
hor_prof2_conc = cat(hor_prof2_1, hor_prof2_2, dims=2)
vert_prof_conc = cat(vert_prof_1, vert_prof_2, dims=2)

saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_hor1.nc", hor_prof1_conc)
saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_hor2.nc", hor_prof2_conc)
saveprofiletonetcdf("/home/haugened/Documents/data/hefex2/converted_for_Rebecca/230822_185513_2325_ver1.nc", vert_prof_conc)
=#
