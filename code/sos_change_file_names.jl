######################################################
###     CHANGE FILE NAMES FOR DANNY'S IR DATA      ###
###        author: Michi Haugeneder                ###
######################################################
#=
Change the file names of the files that I got from Danny Hogan (5.11.24)
so that they appear in the correct order in IRBIS (adding a leading "0")
=#

location = "/home/michi/Documents/slf/tmp/20230124_ir_surface/"

#read all files in folder
files = readdir(location)

#extract file number
filenumber = zeros(Int, length(files))
for i in 1:length(filenumber)
    filenumber[i] = parse(Int, files[i][14:end-4])
end

#add leading zero
newfilenumber = lpad.(string.(filenumber), 4, '0')

#concentate with beginning to obtain new file names
newfilenames = string.(files[1][1:13], newfilenumber)
totalnewfilenames = string.(newfilenames, ".irb")

#find out which file names to update
to_update = .!(files .== totalnewfilenames)

#rename
for i in 1:length(files)
    if to_update[i]
        mv(joinpath(location, files[i]), string(joinpath(location, newfilenames[i]), ".irb"), force=true)
    end
end