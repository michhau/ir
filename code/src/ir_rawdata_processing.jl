######################################################
###    MODULE HANDLING IR-RAWDATA: CONVERTING TO   ###
###             NetCDF4-ARRAYS                     ###
###        author: Michi Haugeneder                ###
######################################################
module irraw

using ReadWriteDlm2, DelimitedFiles, Dates, NCDatasets
#using MAT

export sortandcheckforcompleteness, sequencesin, fromdirloadpointdata,
creatematarray, creatematarrayIRB2ASCII, createnetcdfIRB2ASCII,
saveprofiletonetcdf, saveavgprofiletonetcdf,
readprofilefromnetcdf, readavgprofilefromnetcdf

"""
    sorttxtfilestofolders(pathtofolderall::String)

Sort the unsorted .txt files that were exported from IRBIS 3 plus to
'pathtotopfolderall' into the suited ordered folders (referring to sequence). Not exposed.
"""
function sorttxtfilestofolders(pathtofolderall::String)
    println("Sorting files to corresponding subfolders")
    #regex for the 'xx_points_slice.txt'-files
    sliceregex = string(".*", "_points_slice.txt")
    elementsrawdir = readdir(pathtofolderall)
    println(elementsrawdir)
    for i in elementsrawdir
        if isfile(joinpath(pathtofolderall, i))
            #seperate file stam from numeration (skip file ending)
            elemnumber = 0
            elemstam = ""
            for j in 4:12
                if tryparse(Int, string(i[end-j])) == nothing
                    elemnumber = tryparse(Int, i[end-j+1:end-4])
                    if elemnumber == nothing || elemnumber == 0
                        println("ERROR: extracting numbering in files resulted in error. Please check!")
                        return
                    end
                    elemstam = i[1:end-j-1]
                    break
                end
            end
            if !isdir(joinpath(pathtofolderall, elemstam))
                mkpath(joinpath(pathtofolderall, elemstam))
            end
            if !occursin(Regex(sliceregex), i)
                mv(joinpath(pathtofolderall, i),
                joinpath(pathtofolderall, elemstam, string(lpad(elemnumber, 7, "0"), ".txt")))
            end
        end
    end
    println("Done")
end

"""
    checkforcompleteness(pathtofolder::String)

Check the sorted folders if there are missing files,
sorttxtfilestofolders needs to be executed before. Not exposed.
"""
function checkforcompleteness(pathtofolder::String)
    println("Checking for completeness and writing start/end idx-file")
    #regex for the 'xx_points_slice.txt'-files
    sliceregex =string(".*", "_points_slice.txt")
    nrsequence = readdir(pathtofolder)
    for x in nrsequence
        if !occursin(Regex(sliceregex), x)
            newpathtofolder = joinpath(pathtofolder,x)
            rm(joinpath(newpathtofolder,"idx.txt"), force=true)
            files = readdir(newpathtofolder)
            nrfiles = length(files)
            startfile = files[1]
            endfile = files[end]
            idxstartfile = parse(Int, startfile[1:end-4])
            idxendfile = parse(Int, endfile[1:end-4])
            nrshould = idxendfile - idxstartfile + 1
            if nrshould != nrfiles
                println("Seq:", x, "  ERROR: File(s) missing")
            else
                file = open(joinpath(newpathtofolder,"idx.txt"), "w")
                write(file, string(lpad(idxstartfile, 7, '0'), "\r\n"))
                write(file, string(lpad(idxendfile, 7, '0'), "\r\n"))
                write(file, "1: index startfile;  2: index endfile\r\n")
                close(file)
            end
        end
    end
    println("Done with checking for completeness and writing start/end idx-file")
end

"""
    sortandcheckforcompleteness(folder::String)

Combine 'sorttxtfilestofolder'- and 'checkforcompleteness'-functions
"""
function sortandcheckforcompleteness(folder::String)
    sorttxtfilestofolders(folder)
    checkforcompleteness(folder)
end

"""
    is_number(a::String)::Bool

Return 'true', if argument is number, otherwise 'false' ('10' vs. '2_')
"""
function is_number(a::String)::Bool
    return tryparse(Int64, a) !== nothing
end

"""
    sequencesin(parent_dir::String)

Create an Int-vector with all folder names (sequences) in the
'parent_dir'-directory
"""
function sequencesin(parent_dir::String)
    stringseqs = filter(isdir, readdir(parent_dir, join=true))
    return stringseqs
end

"""
    fromdirloadpointdata(pathtofile::String, seqnr::Int64)

Read coordinates for excerpt from file and return
"""
function fromdirloadpointdata(pathtofile::String, seqnr::Int64)
    file = open(joinpathpath(pathtofile, string(seqnr, "_points_slice.txt")), "r")
    temp = read(file, String)
    #corners: 1st dim: screen nr (1=left, 2=right)
    #         2nd dim: 1=lefttop, 2=righttop, 3=rightbottom, 4=leftbottom
    #         3rd dim: 1=rowindex, 2=colindex
    corners = zeros(Int, 2, 4, 2)
    #scale: 1st dim: 1=scale-row (cm/pxl), 2=scale-col (cm/pxl)
    scale = zeros(Float16, 2)

    #filling the points
    #left screen
    corners[1, 1, 1]  = parse(Int, temp[1:4])
    corners[1, 1, 2]  = parse(Int, temp[6:9])
    corners[1, 2, 1]  = parse(Int, temp[11:14])
    corners[1, 2, 2]  = parse(Int, temp[16:19])
    corners[1, 3, 1]  = parse(Int, temp[21:24])
    corners[1, 3, 2]  = parse(Int, temp[26:29])
    corners[1, 4, 1]  = parse(Int, temp[31:34])
    corners[1, 4, 2]  = parse(Int, temp[36:39])
    #right screen
    corners[2, 1, 1]  = parse(Int, temp[42:45])
    corners[2, 1, 2]  = parse(Int, temp[47:50])
    corners[2, 2, 1]  = parse(Int, temp[52:55])
    corners[2, 2, 2]  = parse(Int, temp[57:60])
    corners[2, 3, 1]  = parse(Int, temp[62:65])
    corners[2, 3, 2]  = parse(Int, temp[67:70])
    corners[2, 4, 1]  = parse(Int, temp[72:75])
    corners[2, 4, 2]  = parse(Int, temp[77:80])
    #scale
    scale[1] = parse(Float16, temp[83:87])
    scale[2] = parse(Float16, temp[89:93])
    #transitionpoint-col
    tpcol = parse(Int, temp[95:98])

    close(file)
    return corners, scale, tpcol
end

"""
    fromdirloadsupplementary(pathtofile::String)

Reads idxstart and idxend from idx.txt-file created by
'checkforcompleteness'-function
"""
function fromdirloadsupplementary(pathtofile::String)
    file = open(joinpath(pathtofile, "idx.txt"), "r")
    temp = read(file, String)
    idxstart = parse(Int, temp[1:7])
    idxend = parse(Int, temp[10:16])
    close(file)
    return idxstart, idxend
end

"""
    fromIRloadIRpicture(filename::String, skiprows::Int64)::Array

Read single IRBIS-exported ASCII-picture and format the resulting array
"""
function fromIRloadIRpicture(filename::String, skiprows::Int64)::Array
    a=readdlm2(filename, '\t'; decimal=',',
    skipstart=9+skiprows, use_mmap=true, dims=(768-skiprows, 1025))
    a=Array{Float16}(a[:, 1:1024])
    return a
end

"""
    fromIRloadIRpictureIRB2ASCII(filename::String, skiprows::Int64)::Array

Read single IRB2ASCII-exported ASCII-picture and format the resulting array
"""
function fromIRloadIRpictureIRB2ASCII(filename::String, skiprows::Int64)::Array
    a=readdlm2(filename, '\t'; decimal='.',
    skipstart=1+skiprows, use_mmap=true)#, dims=(768-skiprows, 1025))
    a=Array{Float16}(a[:, 1:end-1])
    return a
end

"Read IRBIS-exported sequence of ASCII-pictures with
'fromIRloadIRpicture'-function"
function fromIRloadIRsequence(strseq::String,
    idxstartpic::Int64, idxendpic::Int64, skiprows::Int64)::Array
    println("loading IR sequence ", strseq, " ...")
    idxpics = collect(idxstartpic:idxendpic)
    totalnrpictures = length(idxpics)
    IRdata = zeros(Int16, 768-skiprows, 1024, totalnrpictures)
    Threads.@threads for i in 1:totalnrpictures
        filename = joinpath(strseq, string(lpad(idxpics[i],7,"0"), ".txt"))
        IRdata[:, :, i] = round.(Int, fromIRloadIRpicture(filename, skiprows))
        println(i, "/", totalnrpictures)
    end
    println("Done")
    return IRdata
end

"""
    fromIRloadIRsequenceIRB2ASCII(strseq::String,
    idxstartpic::Int64, idxendpic::Int64, skiprows::Int64)::Array

Read IRB2ASCII-exported sequence of ASCII-pictures with
'fromIRloadIRpictureIRB2ASCII'-function
"""
function fromIRloadIRsequenceIRB2ASCII(strseq::String,
    idxstartpic::Int64, idxendpic::Int64, skiprows::Int64)::Array
    println("loading IR sequence ", strseq, " ...")
    idxpics = collect(idxstartpic:idxendpic)
    totalnrpictures = length(idxpics)
    filename = joinpath(strseq, string(lpad(idxpics[1],7,"0"), ".txt"))
    tmp = fromIRloadIRpictureIRB2ASCII(filename, skiprows)
    IRdata = zeros(Int16, size(tmp,1), size(tmp,2), totalnrpictures)
    IRdata[:,:,1] = round.(Int, tmp.*100)
    tmp = nothing
    Threads.@threads for i in 2:totalnrpictures
        filename = joinpath(strseq, string(lpad(idxpics[i],7,"0"), ".txt"))
        IRdata[:, :, i] = round.(Int, fromIRloadIRpictureIRB2ASCII(filename, skiprows).*100)
        println(i, "/", totalnrpictures)
    end
    println("Done")
    return IRdata
end

#=
"Export the obtained 3-dim array (x,y: picture; z:sequence) to .MAT file"
function exportseqtoMAT(a::Array, stringseq::String, fileprefix::String)
    println("Exporting sequence to .mat file...")
    file = matopen(joinpath("..", string(stringseq, ".mat")), "w")
    write(file, "IRdata", a)
    close(file)
    println("Done")
end
=#

"""
    saveirasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)

Save 3D IRdata-array as NetCDF4
"""
function saveirasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)
    @info("Saving IRdata output to NetCDF4-file")
    ds = NCDataset(target, "c")
    defDim(ds, "row", size(data, 1))
    defDim(ds, "col", size(data, 2))
    defDim(ds, "frame", size(data, 3))
    defVar(ds, name, data, ("row", "col", "frame"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

#=
"Go to directory and create .mat array for a predefined excerpt of the IR-image"
function creatematarray(strseq::String, fileprefix::String)
    println(string("Going to directory ", strseq, " ..."))
    cd(strseq)
    (idxstartpic, idxendpic) = fromdirloadsupplementary(strseq)
    #(starttime, endtime) = importfromfiles.fromdirloadsupptime(pathstring)
    #framespersec =, "data\\"
    #(idxendpic - idxstartpic + 1)/Dates.value(Second(endtime - starttime))
    #importfromfiles.exportsupplementary(evaluationfolder, numseq, idxstartpic,
    #idxendpic, Time(starttime), Time(endtime), framespersec)
    IRdata = fromIRloadIRsequence(strseq, idxstartpic, idxendpic,0)
    exportseqtoMAT(IRdata, strseq, fileprefix)
end
=#

#=
"Go to directory and create .mat array for exported with IRB2ASCII"
function creatematarrayIRB2ASCII(strseq::String, fileprefix::String)
    println(string("Going to directory ", strseq, " ..."))
    cd(strseq)
    (idxstartpic, idxendpic) = fromdirloadsupplementary(strseq)
    #(starttime, endtime) = importfromfiles.fromdirloadsupptime(pathstring)
    #framespersec =, "data\\"
    #(idxendpic - idxstartpic + 1)/Dates.value(Second(endtime - starttime))
    #importfromfiles.exportsupplementary(evaluationfolder, numseq, idxstartpic,
    #idxendpic, Time(starttime), Time(endtime), framespersec)
    IRdata = fromIRloadIRsequenceIRB2ASCII(strseq, idxstartpic, idxendpic,0)
    exportseqtoMAT(IRdata, strseq, fileprefix)
end
=#

"""
    createnetcdfIRB2ASCII(strseq::String, target::String, deflatelvl::Int64=5)


    createnetcdfIRB2ASCII(strseq::String, target::String)

Go to directory and create .nc for exported with IRB2ASCII
"""
function createnetcdfIRB2ASCII(strseq::String, target::String, deflatelvl::Int64=5)
    println(string("Going to directory ", strseq, " ..."))
    cd(strseq)
    (idxstartpic, idxendpic) = fromdirloadsupplementary(strseq)
    #(starttime, endtime) = importfromfiles.fromdirloadsupptime(pathstring)
    #framespersec =, "data\\"
    #(idxendpic - idxstartpic + 1)/Dates.value(Second(endtime - starttime))
    #importfromfiles.exportsupplementary(evaluationfolder, numseq, idxstartpic,
    #idxendpic, Time(starttime), Time(endtime), framespersec)
    IRdata = fromIRloadIRsequenceIRB2ASCII(strseq, idxstartpic, idxendpic,0)
    saveirasnetcdf(IRdata, "irdata", target, deflatelvl)
end

"""
    saveprofiletonetcdf(target::String, data::Array, deflatelvl::Int64=5)

save profiles to NetCDF4
"""
function saveprofiletonetcdf(target::String, data::Array, deflatelvl::Int64=5)
    @show("Saving profile to netCDF")
    ds = NCDataset(target, "c")
    defDim(ds, "pxl", size(data, 1))
    defDim(ds, "frame", size(data, 2))
    defVar(ds, "profile", data, ("pxl", "frame"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

"""
    saveavgprofiletonetcdf(target::String, z::Vector, data_mean::Vector, data_low_quart::Vector, data_up_quart::Vector, deflatelvl::Int64=5)

save time averaged profiles to NetCDF4
"""
function saveavgprofiletonetcdf(target::String, z::Vector, data_mean::Vector, data_low_quart::Vector, data_up_quart::Vector, deflatelvl::Int64=5)
    @show("Saving time-averaged profile to netCDF")
    ds = NCDataset(target, "c")
    defDim(ds, "zdim", size(z, 1))
    defVar(ds, "z", z, ("zdim",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "T_avg", data_mean, ("zdim",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "T_25quart", data_low_quart, ("zdim",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "T_75quart", data_up_quart, ("zdim",); shuffle=true, deflatelevel=deflatelvl)    
    close(ds)
end

"""
    readprofilefromnetcdf(dir::String)

Read profile from netCDF
"""
function readprofilefromnetcdf(dir::String)
    netfile = NCDataset(dir, "r")
    dat = netfile["profile"]
    dattot = dat[:,:]
    close(netfile)
    return dattot
end

"""
    readavgprofilefromnetcdf(dir::String)

Read time averaged profile (average, 25%, and 75% quartile) from netCDF
"""
function readavgprofilefromnetcdf(dir::String)
    netfile = NCDataset(dir, "r")
    z = netfile["z"]
    z = z[:]
    prof_avg = netfile["T_avg"]
    prof_avg = prof_avg[:]
    prof_low_quart = netfile["T_25quart"]
    prof_low_quart = prof_low_quart[:]
    prof_up_quart = netfile["T_75quart"]
    prof_up_quart = prof_up_quart[:]
    close(netfile)
    return z, prof_avg, prof_low_quart, prof_up_quart
end

end #module
