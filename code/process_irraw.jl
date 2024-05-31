######################################################
###      CONVERTING IR-RAWDATA TO MAT-ARRAYS       ###
###        author: Michi Haugeneder                ###
######################################################
#=
notes:
 scans complete directory for IRBIS-exported .txt-files and creates .mat arrays
=#

importdir = joinpath(@__DIR__)
pathtofile = "/home/haugened/Documents/data/hefex2/IR_video/230822_1537_01/bea"
include(joinpath(importdir, "src", "ir_rawdata_processing.jl"))
import .irraw

##########################################################
###            VARIABLES TO BE CHANGED                 ###
##########################################################
#evaluationfolder = joinpath(pathtofile,"ir_data_spring21/210428_111646/")

evaluationlist = [pathtofile]
println()
println("-----------S-T-A-R-T-------------")

println("IRB2ASCII: Decimal seperator ., between values tab!!")

##########################################################
###     CONVERT SINGLE .TXT-FILES TO .MAT ARRAYS       ###
##########################################################
for evaluationfolder in evaluationlist
    #Sort the exported .txt-files into folders and check for completeness
    irraw.sortandcheckforcompleteness(evaluationfolder)

    #create vector with subfolder names
    stringseq = irraw.sequencesin(evaluationfolder)
    for i in stringseq
        #irraw.creatematarrayIRB2ASCII(i, "a")
        targetstring = joinpath(evaluationfolder, string(i, ".nc"))
        irraw.createnetcdfIRB2ASCII(i, targetstring)
    end
    #=
   #only SLFLINUX24: running the MATLAB-script to compress the output .mat-files
    println("running MATLAB script to compress exported .mat-files...")
    command = string("-nodisplay -r \"cd('/home/haugened/Documents/MATLAB/'); resave('", evaluationfolder, "'); quit\"");
    matcommand = `matlab $command`;
    run(matcommand);=#
    println("Done")
end

println("-----------D-O-N-E---------------")
##########################################################
