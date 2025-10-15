
#directory containing subfolders with 1000 files in each
upperdir = "/home/haugened/Documents/data/sos/250917_for_Eli/txt"

#folders = readdir(upperdir)
#######################################################
#generate subfolders each containing defined number of frames
framesperfolder = 3000
files = readdir(upperdir)
nfiles = length(files)
nfolders = ceil(Int, nfiles / framesperfolder)
println("Number of files: ", nfiles)
println("Number of folders to be created: ", nfolders)
for i in 1:nfolders
    foldername = lpad(string(i), 3, '0')
    mkdir(joinpath(upperdir, foldername))
end
files = sort(files)
for (k, i) in enumerate(files)
    folderindex = div(k - 1, framesperfolder) + 1
    foldername = lpad(string(folderindex), 3, '0')
    mv(joinpath(upperdir, i), joinpath(upperdir, foldername, i); force=true)
end
println("Done")
######################################################
#=
#move files from subfolders to upper directory
for i in folders
    println("Processing folder: ", i)
    files = readdir(joinpath(upperdir, i))
    for j in files
        mv(joinpath(upperdir, i, j), joinpath(upperdir, j))
    end
end

#remove empty subfolders
for i in folders
    rm(joinpath(upperdir, i); force=true, recursive=true)
end
=#
#=
#thin out files, keep only every keepevery-th file
keepevery = 30
files = readdir(upperdir)
nfiles = length(files)
println("Number of files before thinning: ", nfiles)
for (k, i) in enumerate(files)
    if mod(k, keepevery) != 1
        rm(joinpath(upperdir, i); force=true)
    end
end
files = readdir(upperdir)
nfiles = length(files)
println("Number of files after thinning: ", nfiles)
println("Done")
=#
######################################################
#=
#compare to converted files
dir_conv = "/home/haugened/Documents/data/sos/250917_for_Eli/converted/230518/"
files_conv = readdir(dir_conv)
nfiles_conv = length(files_conv)
println("Number of converted files: ", nfiles_conv)
println("Difference in number of files: ", nfiles - nfiles_conv)
println("Done")

#find and print missing files (discard last 4 characters of files and files_conv, which are .txt and .irb)
missing = String[]
for i in files
    short_i = i[1:end-4]
    found = false
    for j in files_conv
        short_j = j[1:end-4]
        if short_i == short_j
            found = true
            break
        end
    end
    if !found
        push!(missing, i)
    end
end
println("Missing files: ")
for i in missing
    println(i)
end
println("Number of missing files: ", length(missing))
println("Done")
=#
######################################################
#move irdataxxxxx.txt to xxxxx.txt
for i in readdir(upperdir)
    println("folder: ", i)
    files = readdir(joinpath(upperdir,i))
    #rename files to consecutive numbering starting from 00001.txt
    nfiles = length(files)
    sortedfiles = sort(files)
    for (k, j) in enumerate(sortedfiles)
        newname = lpad(string(k), 7, '0') * ".txt"
        mv(joinpath(upperdir, i, j), joinpath(upperdir, i, newname); force=true)
    end
    files = readdir(upperdir)
end