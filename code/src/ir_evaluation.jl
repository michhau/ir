######################################################
###          MODULE FOR IR-DATA-EVALUATION         ###
###            author: Michi Haugeneder            ###
######################################################
module irev

using ReadWriteDlm2, DelimitedFiles, Dates, PyCall,
    ImageFiltering, ProgressMeter, LaTeXStrings, NCDatasets,
    Statistics, ProgressMeter, Printf
#not on HYPERION
if (gethostname() == "slfl29682" || gethostname() == "Michi-T450s" || gethostname() == "x1carbon5")
    using PyPlot
    cm = pyimport("matplotlib.cm")
    cramericm = pyimport("cmcrameri.cm")
    patches = pyimport("matplotlib.patches")
    mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
end
#include("general.jl")
#include("fourier.jl")
#import .gen
#import .ft

export loadfromNetCDF4, loadexcerptfromNetCDF4, loadcolsfromNetCDF4, findfiles, saveirasnetcdf,
    fromIRloadwinddata, fromdirloadsupptime,
    exportsupplementary, subtracttimeaverage, subtractbacklookingtimeaverage, setsnowsurfacetonan,
    spatialgaussfilter, completekagapreprocessing, getdims, showheatmap, showheatmaptrad,
    removesnowsurface, averagewindfield, fetchdependfunc,
    ProfileSpec, findfiles, plot_first_frame, prepare_multiple_profiles, plot_multi_profiles

######################################################
###                LOADING AND I/O                 ###
######################################################
"""
    loadfromNetCDF4(file::String, varname::String)

Load variable from NetCDF4 file
"""
function loadfromNetCDF4(file::String, varname::String)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = Float16.(dat[:,:,:])
    if eltype(dat) == Int16
        dattot ./= 100.0
    end
    header_Dict = Dict{String, Any}()
    #check for header
    try
        sourcefile = netfile["sourcefile"]
        header_Dict["sourcefile"] = sourcefile[:]
        timestamp = netfile["timestamp"]
        header_Dict["timestamp"] = timestamp[:]
        timestamp_ms = netfile["timestamp_ms"]
        header_Dict["timestamp_ms"] = timestamp_ms[:]
        tempshutter = netfile["temp_shutter"]
        header_Dict["temp_shutter"] = tempshutter[:]
        imageindex = netfile["imageindex"]
        header_Dict["imageindex"] = imageindex[:]
    catch e
    end
    close(netfile)
    return dattot, header_Dict
end

"""
    loadexcerptfromNetCDF4(file::String, varname::String,
    rowrange, colrange)

Load excerpt of variable from NetCDF4 file
"""
function loadexcerptfromNetCDF4(file::String, varname::String,
    rowrange, colrange)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = Float16.(dat[rowrange,colrange,:])
    if eltype(dat) == Int16
        dattot ./= 100.0
    end
    close(netfile)
    return dattot
end

"""
    loadexcerptfromNetCDF4(file::String, varname::String,
    rowrange, colrange)

Load excerpt of variable from NetCDF4 file
"""
function loadexcerptfromNetCDF4(file::String, varname::String,
    rowrange, colrange, framerange)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = Float16.(dat[rowrange,colrange,framerange])
    if eltype(dat) == Int16
        dattot ./= 100.0
    end
    close(netfile)
    return dattot
end

"""
    loadcolsfromNetCDF4(file::String, varname::String,
    colrange)

Load excerpt of variable from NetCDF4 file
"""
function loadcolsfromNetCDF4(file::String, varname::String,
    colrange)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = Float16.(dat[:, colrange, :])
    if eltype(dat) == Int16
        dattot ./= 100.0
    end
    close(netfile)
    return dattot
end

"""
    findfiles(folder::String, pat::Regex)::Vector{String}

Return vector containing files in 'folder' that match
Regular Expression 'pat'
"""
function findfiles(folder::String, pat::Regex)::Vector{String}
    files = readdir(folder)
    return files[occursin.(pat, files)]
end

"""
    saveirasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)

Saving 3D IRdata-array as NetCDF4
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

"""
    fromIRloadwinddata(
    filename::String, startdate::DateTime, enddate::DateTime)

Load and format ultrasonic-data from Logger and return DateTime-vector,
data array, and a columnname vector
"""
function fromIRloadwinddata(
    filename::String, startdate::DateTime, enddate::DateTime)
    @info("loading wind data from IR-ultrasonic...")
    df = readdlm(
        string(filename), ',', String, '\n';
        skipstart=4, use_mmap=true)
    #extract first column to be DateTime
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    #convert to Float64
    df = parse.(Float64, df[startdate.<timeofmeasure.<enddate, 5:6])
    #specify a String vector with the column names
    columnnames = ["Winddirection", "Windspeed"]
    return timeofmeasure[startdate.<timeofmeasure.<enddate], df, columnnames
end

"""
    fromdirloadsupptime(pathtofile::String)

Read starttime and endtime from idx.txt file
"""
function fromdirloadsupptime(pathtofile::String)
    file = open(string(pathtofile, "data\\idx.txt"), "r")
    temp = read(file, String)
    dateformat = DateFormat("yyy-mm-dd HH:MM:SS")
    starttime = DateTime.(temp[10:28], dateformat)
    endtime = DateTime.(temp[31:49], dateformat)
    close(file)
    return starttime, endtime
end

"""
    exportsupplementary(
    pathtofolder::String, nrseq::Integer, idxstart::Integer, idxend::Integer, timestart,
    timeend, timestep::AbstractFloat)

Export .dat file with supplementary information
"""
function exportsupplementary(
    pathtofolder::String, nrseq::Integer, idxstart::Integer, idxend::Integer, timestart,
    timeend, timestep::AbstractFloat)
    @info("Exporting supplementary file...")
    sup = ["          " "INDEX" "TIME"
        "START       " string(idxstart) string(timestart)
        "END       " string(idxend) string(timeend)
        "DURATION  " string(idxend - idxstart + 1) string(Second(timeend - timestart))
        "" "" ""
        "FRAMES/S  " string(timestep) ""]
    filename =
        string(string(nrseq), "_supplementary.dat")
    #string(string(nrseq), "_", string(idxstart), "to", string(idxend),"_supplementary.dat")
    open(string(pathtofolder, filename), "w") do io
        writedlm(io, sup)
    end
end

#=
"Returning a vector containing all irdataxx.mat files
belonging to the sequence (stored in the seq folder)."
function filesinseq(folder::String)::Vector{String}
    files = readdir(folder)
    re = r"^irdata.*\.mat"
    return files[occursin.(re, files)]
end
=#
######################################################
###                   PREPROCESSING                ###
######################################################
#=
"""
    subtracttimeaverage(IRdatain::Array, nrelements::Integer)::Array

Perform time-averaging for all points as due to 'Inagaki et al. (2013)'
'nrelements' refers to third dimension of 'IRdata' (time-dimension)
"""
function subtracttimeaverage(IRdatain::Array, nrelements::Integer)::Array
    IRdatareturn = similar(IRdatain)
    @info("Performing time averaging...")
    @sync for i in 1:size(IRdatain, 2)
        Threads.@spawn for j in 1:size(IRdatain, 1)
            IRdatareturn[j, i, :] = IRdatain[j, i, :] - gen.movingaverage(IRdatain[j, i, :], nrelements)
        end
    end
    return IRdatareturn
end

"""
    subtractbacklookingtimeaverage(IRdatain::Array, nrelements::Integer)::Array

Perform back-looking time-averaging for all points as due to 'Inagaki et al. (2013)'
'nrelements' refers to third dimension of 'IRdata' (time-dimension)
"""
function subtractbacklookingtimeaverage(IRdatain::Array, nrelements::Integer)::Array
    IRdatareturn = similar(IRdatain)
    @info("Performing back-looking time averaging...")
    @sync for i in 1:size(IRdatain, 2)
        Threads.@spawn for j in 1:size(IRdatain, 1)
            IRdatareturn[j, i, :] = IRdatain[j, i, :] - gen.backlookingmovavg(IRdatain[j, i, :], nrelements)
        end
    end
    println("Done")
    return IRdatareturn
end

"""
    setsnowsurfacetonan(irdata::Array, thresh::AbstractFloat=1.5; framewise=false)

Set the pixels showing the snow surface (see condition) to NaN
and return row nr. for each column
"""
function setsnowsurfacetonan(irdata::Array, thresh::AbstractFloat=1.5; framewise=false)
    @info("Setting snow-surface to NaN")
    irout = copy(irdata)
    if !framewise
        snowsurf = zeros(Int64, size(irdata, 2))
    else
        snowsurf = zeros(Int64, size(irdata, 2), size(irdata, 3))
    end
    #only do it for the lowest 30% assuming top is row=1
    lim = round(Int, size(irdata, 1) * 0.7)
    if !framewise
        @sync for icol in 1:size(irdata, 2)
            Threads.@spawn for jrow in lim:size(irdata, 1)
                if any(x -> x < thresh, irdata[jrow, icol, :])
                    irout[jrow, icol, :] .= NaN
                end
            end
        end
        for jcol in 1:size(irdata, 2)
            for irow in lim:size(irdata, 1)
                if isnan(irdata[irow, jcol, round(Int, (1 + size(irdata, 3)) / 2)])
                    snowsurf[jcol] = irow
                    break
                end
            end
        end
    else #do it frame-by-frame
        irtmp = irdata[lim:end, :, :]
        irtmp[irtmp.<thresh] .= NaN
        irout[lim:end, :, :] = irtmp
        for fr in 1:size(irout, 3)
            for col in 1:size(irout, 2)
                tmp = maximum(filter(x->!isnothing(x), [0 findfirst(x -> !isnan(x), irout[end:-1:1, col, fr])]))
                snowsurf[col, fr] = size(irout, 1) - tmp + 1
            end
        end
    end
    return irout, snowsurf
end

"""
    spatialgaussfilter(irdata::Array, sigmrow::Integer, sigmcol::Integer)

Apply picture-wise spatial filtering using Gauss-kernel with sigmarow and sigmacol
"""
function spatialgaussfilter(irdata::Array, sigmrow::Integer, sigmcol::Integer)
    :Array
    @info("Picture-wise spatial Gauss-filtering with sigma_row=", sigmrow,
        " and sigma_col=", sigmcol)
    out = similar(irdata)
    for i in 1:size(irdata, 3)
        out[:, :, i] = imfilter(irdata[:, :, i], Kernel.gaussian((sigmrow, sigmcol)), NA())
    end
    return out
end

"""
    completekagapreprocessing(irin::Array, gaussparam::Vector{Integer},
    framespersec, backlookingavgtime, fco)

Apply complete preprocessing to (already cut) 'irin' array.
"""
function completekagapreprocessing(irin::Array, gaussparam::Vector{Integer},
    framespersec, backlookingavgtime, fco)
    @info("Preprocessing IR-data for Kaga-windfield-estimation")
    #setting the snow to NaN and spatial Gauss filtering
    (irtmp, snowsurf) = irev.setsnowsurfacetonan(irin)
    irtmp = irev.spatialgaussfilter(irtmp, gaussparam[1], gaussparam[2])

    #high-pass filtering = subtract backlooking time average
    irtmp2 = irev.subtractbacklookingtimeaverage(irtmp, round(Int, backlookingavgtime * framespersec))

    #low-pass filtering using Fourier transformation
    (lpir,) = ft.fouriercutoffmatrix(irtmp2, round(Int, framespersec * size(irtmp2, 3)), fco)
    return lpir, snowsurf
end
=#
######################################################
###                 PROCESSING IR                  ###
######################################################
"""
    getdims(ncfile::String, varname::String)

get dimensions of variable in NetCDF4 file
"""
function getdims(ncfile::String, varname::String)
    ds = NCDataset(ncfile, "r")
    a = ds[varname]
    return size(a)
end

"""
    showheatmap(data::Array, min, max, title::String="IRdata", cmperpxlu=0.606, cmperpxlw=0.549, xis0pxl=30, yis0pxl=360)

Show a heatmap of a 2D-array
"""
function showheatmap(data::Array, min, max, title::String="IRdata", cmperpxlu=0.606, cmperpxlw=0.549, xis0pxl=30, yis0pxl=360)
    fntsze = 18
    lblsze = fntsze - 5
    patches = pyimport("matplotlib.patches")
    PyPlot.pygui(true)
    fig = PyPlot.figure(figsize=(0.469366 * 16, 5.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title)
    xleft = -(xis0pxl - 1) * cmperpxlu / 100
    xright = size(data, 2) * cmperpxlu / 100 + xleft
    ybottom = (yis0pxl - size(data, 1)) * cmperpxlw / 100
    ytop = size(data, 1) * cmperpxlw / 100 + ybottom
    #rect = patches.Rectangle((xleft+302*cmperpxlu/100, (size(data,1)-289)*cmperpxlw/100), 0.103, 0.0933, edgecolor="black", facecolor="None")
    #rect.set_lw(2)
    #ax.add_patch(rect)
    hm = ax.imshow(data, extent=[xleft, xright, ybottom, ytop], cmap="turbo", vmin=min, vmax=max)
    ax.set_xlabel(L"x~\mathrm{[m]}", fontsize=fntsze)
    ax.set_ylabel(L"h~\mathrm{[m]}", fontsize=fntsze)
    ax.tick_params(labelsize=lblsze)
    cbar = fig.colorbar(hm, ax=ax)
    cbar.set_label(L"T~\mathrm{[^\circ C]}", fontsize=fntsze)
    cbar.ax.tick_params(labelsize=lblsze)
    PyPlot.show()
end

"""
    showheatmaptrad(data::Array, min, max, title::String="IRdata")

Show a heatmap of a 2D-array traditional without meter stuff
"""
function showheatmaptrad(data::Array, min, max, title::String="IRdata")
    fntsze = 22
    lblsze = fntsze - 5
    patches = pyimport("matplotlib.patches")
    PyPlot.pygui(true)
    fig = PyPlot.figure(figsize=(16, 5.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title)
    cmap = PyPlot.get_cmap("tab10")
    #rect1 = patches.Rectangle((227+90,98+425), 17, 17, edgecolor=cmap(0), facecolor="None")
    #rect2 = patches.Rectangle((260,190), 10, 10, edgecolor=cmap(1), facecolor="None")
    #rect3 = patches.Rectangle((12,325), 10, 10, edgecolor=cmap(2), facecolor="None")
    #rect4 = patches.Rectangle((800,190), 10, 10, edgecolor=cmap(3), facecolor="None")
    #rect5 = patches.Rectangle((800,300), 10, 10, edgecolor=cmap(4), facecolor="None")
    #rect6 = patches.Rectangle((900,335), 10, 10, edgecolor=cmap(5), facecolor="None")
    #rect1.set_lw(1.5)
    #rect2.set_lw(1.5)
    #rect3.set_lw(1.5)
    #rect4.set_lw(1.5)
    #rect5.set_lw(1.5)
    #rect6.set_lw(1.5)
    #ax.add_patch(rect1)
    #ax.add_patch(rect2)
    #ax.add_patch(rect3)
    #ax.add_patch(rect4)
    #ax.add_patch(rect5)
    #ax.add_patch(rect6)
    #xleft = -(xis0pxl-1)*cmperpxl/100
    #xright = size(data, 2)*cmperpxl/100+xleft
    hm = ax.imshow(data, cmap="turbo", vmin=min, vmax=max)
    PyPlot.axis("off")
    #ax.tick_params(labelsize=lblsze)
    divider = mpl_axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = fig.colorbar(hm, cax=cax)
    cbar.set_label(L"T~\mathrm{[^\circ C]}", fontsize=fntsze)
    cbar.ax.tick_params(labelsize=lblsze)
    PyPlot.show()
end

#=
"""
    removesnowsurface(irin::Array)::Array

Remove snow surface (set to NaN before)
"""
function removesnowsurface(irin::Array)::Array
    #only do it for the lowest 30% assuming top is row=1
    lim = round(Int, size(irin, 1) * 0.7)
    idxsnowsurf = zeros(Int64, size(irin, 2))
    picpos = round(Int, size(irin, 3) / 2)
    for icol in 1:size(irin, 2)
        idxsnowsurf[icol] = findfirst(x -> isnan(x), irin[lim:end, icol, picpos]) + lim - 1
        if isnothing(idxsnowsurf[icol])
            idxsnowsurf[icol] = size(irin, 1)
        end
    end
    highestsnowsurf = minimum(idxsnowsurf)
    irout = zeros(Float64, highestsnowsurf - 1, size(irin, 2), size(irin, 3))
    for icol in 1:size(irin, 2)
        for jpic in 1:size(irin, 3)
            irout[:, icol, jpic] = irin[(-highestsnowsurf+1:-1).+idxsnowsurf[icol], icol, jpic]
        end
    end
    return irout
end
=#
######################################################
###                 POST-PROCESSING                ###
######################################################
#=
"Gliding average for the windfield, 'nrelements' refers to 3rd dim of IRdata"
function averagewindfield(windfield::Array, nrelements::Integer)::Array
    @info("Gliding average for windfield")
    avgwind = similar(windfield)
    for icol in 1:size(windfield, 2)
        for jcol in 1:size(windfield, 1)
            avgwind[jcol, icol, :, 1] = gen.movingaverage(windfield[jcol, icol, :, 1], nrelements)
            avgwind[jcol, icol, :, 2] = gen.movingaverage(windfield[jcol, icol, :, 2], nrelements)
        end
    end
    return avgwind
end
=#
######################################################
###                   EVALUATION                   ###
######################################################
"""
    fetchdependfunc(f::Function, irsize, snowsurf::Vector, fetchis0::Integer, cmperpxl::AbstractFloat)::Array

Create a value table for a fetch-dependend function 'f'. 'fetchis0' is the 0 point;
fetch increases to the right. Return array with 1st col x, 2nd col y.
"""
function fetchdependfunc(f::Function, irsize, snowsurf::Vector, fetchis0::Integer, cmperpxl::AbstractFloat)::Array
    xleftmeter = (1 - fetchis0) * cmperpxl * 0.01
    xrightmeter = (irsize[2] - fetchis0) * cmperpxl * 0.01
    xinmeter = collect(xleftmeter:cmperpxl*0.01:xrightmeter)
    for i in 1:length(xinmeter)
        if xinmeter[i] < 0.0
            xinmeter[i] = NaN
        end
    end
    res = f.(xinmeter)
    for i in 1:length(res)
        if snowsurf[i] > 270
            res[i] = (irsize[1] - snowsurf[i] - 1) * cmperpxl / 100 + res[i]
        else
            res[i] = NaN
        end
    end
    return hcat(xinmeter, res)
end

######################################################
###                    CONTRASTS                   ###
######################################################

"""
    struct ProfileSpec
        row_range :: UnitRange{Int}   # rows that form the vertical axis
        col_range :: UnitRange{Int}   # columns that will be collapsed (median)
        label     :: String           # legend entry (e.g. “Site A”)
    end
"""
struct ProfileSpec
    row_range :: UnitRange{Int}
    col_range :: UnitRange{Int}
    label     :: String
end

"""Convert a folder name like `250412_083015` → DateTime (UTC)."""
function foldername_to_dt(name::String)::DateTime
    parsed_date = DateTime(1900,01,01,00,00,00)
    try
        parsed_date = DateTime(name, dateformat"yymmdd_HHMMSS") + Year(2000)
    catch e
         @debug "Folder $(name) does not follow the expected naming scheme."
    end
    return parsed_date
end

"""
    file_index(fname::String) → Int

`fname` is something like `irdata_0012.nc`.  Returns the integer index
(`12`).  Throws an error if the pattern does not match.
"""
function fileindex(fname::String)
    m = match(r"irdata_(\d{2,})\.nc$", fname)
    m === nothing && error("Unexpected NetCDF filename: $fname")
    return parse(Int, m.captures[1])
end

"""
    find_files(root_dir::String,
               t_start::DateTime,
               t_end::DateTime) → Vector{String}

Return a **sorted** vector of absolute paths to NetCDF files that contain at
least one frame whose timestamp lies inside `[t_start, t_end]`.

Assumptions that the algorithm exploits
----------------------------------------
* Each folder is named `yyyymmdd_HHMMSS` – the start time of the *first*
  NetCDF file (`irdata_0000.nc`) inside that folder.
* Inside a folder the files are named `irdata_XXXX.nc` where `XXXX` is a
  zero‑padded counter that increases monotonically with time (≈ 16–17 s
  between successive files).
* The NetCDF variable `timestamp` stores “days since 1900‑01‑01”.
"""
function find_files(root_dir::String, t_start::DateTime, t_end::DateTime)::Vector{String}
    @info "Scanning $(root_dir) …"

    # ------------------------------------------------------------------
    # 1️⃣  Gather candidate folders (those that start early enough)
    # ------------------------------------------------------------------
    candidate_folders = String[]
    for (folder_root, dirs, _) in walkdir(root_dir; follow_symlinks=false)
        for d in dirs
            dt = foldername_to_dt(d)
            if t_start - Day(2) ≤ dt ≤ t_end
                push!(candidate_folders, joinpath(folder_root, d))
            end
        end
    end

    if isempty(candidate_folders)
        @warn "No candidate folders found under $(root_dir)."
        return String[]
    end

    # Sort folders chronologically (helps deterministic processing)
    sort!(candidate_folders)

    # ------------------------------------------------------------------
    # 2️⃣  Walk through the candidate folders in chronological order
    # ------------------------------------------------------------------
    matched_paths = String[]

    for folder in candidate_folders
        folder_start = foldername_to_dt(basename(folder))
        @debug "Processing folder $(folder) (starts $(folder_start))"

        # ----------------------------------------------------------------
        # 2a️⃣  List the NetCDF files inside the folder (already sorted)
        # ----------------------------------------------------------------
        nc_files = filter(f -> endswith(f, ".nc"),
                          sort(readdir(folder; join=true)))

        # Early‑out flag: once we have seen a file whose *first* timestamp >
        # t_end we can break out of the inner loop (no later file can be useful).
        stop_folder = false

        for nc_path in nc_files
            # ----------------------------------------------------------------
            # 2b️⃣  Quick numeric check using the file name index
            # ----------------------------------------------------------------
            file_index = fileindex(nc_path)

            # Approximate timestamp of the *first* frame in this file:
            #   folder_start + file_index * Δt_per_file
            # Δt_per_file ≈ 500 frames / 30 Hz ≈ 16.666… s
            Δt_per_file = Second(round(Int, 500 / 30))   # 16 s (rounded)
            approx_file_start = folder_start + file_index * Δt_per_file

            if approx_file_start > t_end
                @debug "Skipping remaining files in $(folder) – start $(approx_file_start) > t_end"
                stop_folder = true
                break
            end

            # ----------------------------------------------------------------
            # 2c️⃣  Read ONLY the first timestamp from the NetCDF file
            # ----------------------------------------------------------------
            ds = Dataset(nc_path, "r")
            try
                ts = ds["timestamp"][:][1] 
                if t_start ≤ ts ≤ t_end
                    push!(matched_paths, nc_path)
                elseif ts > t_end
                    # This file already starts after the interval → stop scanning
                    stop_folder = true
                    break
                end
                # If ts_dt < t_start we simply continue – a later file may match.
            finally
                close(ds)
            end
        end   # ← file loop

        stop_folder && break   # no later folder can contain relevant data
    end   # ← folder loop

    sort!(matched_paths)
    @info "Found $(length(matched_paths)) matching NetCDF files."
    return matched_paths
end

# ----------------------------------------------------------------------
# Load a single frame (as Float64) and apply the conversion factor
# ----------------------------------------------------------------------
"""
    read_frame(path::String, frame_idx::Int) -> Matrix{Float64}

Loads `irdata[:,:,frame_idx]` from `path` and divides by 100 to obtain °C.
Only the requested slice is read, keeping memory usage minimal.
"""
function read_frame(path::String, frame_idx::Int)::Matrix{Float64}
    ds = Dataset(path, "r")
    try
        # NetCDF stores data in C order (row‑major); NCDatasets returns an
        # Array with the same ordering, so slicing works directly.
        raw = ds["irdata"][:, :, frame_idx]   # 2‑D slice
        return raw ./ 100.0                    # convert to °C
    finally
        close(ds)
    end
end

# ----------------------------------------------------------------------
# Plot the first matching frame for interactive selection
# ----------------------------------------------------------------------
function plot_first_frame(filelist::Vector{String}; tmin=-99.9, tmax=99.9)
    @info "Loading first frame of the earliest file …"
    first_file = filelist[1]
    frame = read_frame(first_file, 1)   # frame index 1 (NetCDF is 1‑based)
    (fig, ax) = PyPlot.subplots()
    if tmin != -99.9 || tmax != 99.9
        im = ax.imshow(frame, origin="upper", cmap=cramericm.batlow, vmin=tmin, vmax=tmax)#; cmap=cmap, origin="upper")
    else
        im = ax.imshow(frame, origin="upper", cmap=cramericm.batlow)#; cmap=cmap, origin="upper")
    end
    ax.set_title("First frame (file: $(basename(first_file)))")
    fig.colorbar(im, ax=ax, label="Temperature (°C)")
    return fig, ax, frame
end

"""
    count_matching_frames_opt(filelist::Vector{String},
                             t_start::DateTime,
                             t_end::DateTime) → Int

Counts how many frames lie inside `[t_start, t_end]` **without opening
most files**.  It assumes that every file except the *last file in its
folder* (and the *last file of the whole selection*) contains exactly
`FRAMES_PER_FILE` frames.  Those two exceptional files are opened to read
their actual `timestamp` length.

The function returns the total number of frames that intersect the interval.
"""
function count_matching_frames(filelist::Vector{String},
                                  t_start::DateTime,
                                  t_end::DateTime,
                                  secs_per_file::Real,
                                  frames_per_file::Int)::Int
    # ------------------------------------------------------------------
    # 1️⃣  Group files by folder (they are already sorted alphabetically)
    # ------------------------------------------------------------------
    sort!(filelist)                                 # ensure deterministic order
    folder_groups = Dict{String,Vector{String}}()

    for p in filelist
        folder = dirname(p)                         # parent directory
        push!(get!(folder_groups, folder, String[]), p)
    end

    total_frames = 0

    # ------------------------------------------------------------------
    # 2️⃣  Process each folder
    # ------------------------------------------------------------------
    for (folder, files) in sort(collect(folder_groups); by=first)
        folder_dt = foldername_to_dt(folder)
        folder_dt === nothing && continue           # skip malformed folders

        # Files inside a folder are already sorted because we sorted `filelist`
        # The *last* file in this folder may be truncated.
        last_file_in_folder = files[end]

        for (i, fpath) in enumerate(files)
            idx = fileindex(basename(fpath))

            # --------------------------------------------------------------
            # Compute the *theoretical* start‑time of this file
            # --------------------------------------------------------------
            file_start = folder_dt + Millisecond(round(Int,
                                 idx * secs_per_file * 1000))

            # By default we assume a full‑size file
            n_frames_this_file = frames_per_file

            # --------------------------------------------------------------
            # Detect the two “exceptional” files that need real inspection
            # --------------------------------------------------------------
            is_last_in_folder = (fpath == last_file_in_folder)
            is_last_overall   = (fpath == filelist[end])

            if (is_last_in_folder) .&& .!(is_last_overall)
                # Open the file and read the actual length of the timestamp vector
                ds = Dataset(fpath, "r")
                try
                    n_frames_this_file = length(ds["timestamp"][:])
                finally
                    close(ds)
                end
            end
            if is_last_overall
                ds = Dataset(fpath, "r")
                try
                    timestamps = ds["timestamp"][:]

                    # Identify indices within the interval
                    idx_in = findall(t -> t ≤ t_end, timestamps)
                    n_frames_this_file = length(idx_in)
                finally
                    close(ds)
                end
            end

            frames_here = max(0, n_frames_this_file)
            total_frames += frames_here
        end
    end

    @info "Counted $total_frames frames that intersect the requested interval."
    return total_frames
end

"""
    read_and_collapse!(
        profile_data::Matrix{Float64}
        filelist::Vector{String},
        t_start::DateTime,
        t_end::DateTime,
        row_range::UnitRange{Int},
        col_range::UnitRange{Int},
    ) → Nothing

Fills `profile_data` column‑wise.  Each column corresponds to one frame.
If `length(col_range) > 1` the function first computes the **row‑wise median**
across the selected columns; otherwise the single column is copied unchanged.
Only the rows in `row_range` are ever read from disk.
"""
function read_and_collapse!(profile_data::AbstractArray{Float64},
                           filelist::Vector{String},
                           t_start::DateTime,
                           t_end::DateTime,
                           row_range::UnitRange{Int},
                           col_range::UnitRange{Int})::Nothing
    col_len = length(col_range)
    cur_col = 1                               # which column of profile_data we are filling

    total_frames = size(profile_data, 2)       # number of columns we allocated
    prog = Progress(total_frames; desc="Reading frames", dt=0.5)

    for path in filelist
        ds = Dataset(path, "r")
        try
            ts_vec = ds["timestamp"][:]

            # Find the indices of frames that belong to the interval
            idxs = findall(t -> t_start .≤ t .≤ t_end, ts_vec)

            # If there are no matching frames in this file, skip it
            isempty(idxs) && continue

            # ----------------------------------------------------------------
            # Load **only** the rows/columns we need, for *all* matching frames
            # ----------------------------------------------------------------
            # NetCDF slicing: var[:, :, idx]  →  (row, col, frame)
            # We request a view of shape (nRows, nCols, nFrames)
            # NCDatasets allows us to pass a tuple of ranges.
            # NOTE: Julia uses 1‑based indexing, same as NetCDF.
            # ----------------------------------------------------------------
            ir = ds["irdata"]
            sub = Float32.(ir[row_range, col_range, idxs])   # 3‑D Array{Float32,3} (or Float64)

            # Convert to °C once (divide by 100)
            sub ./= 100.0

            n_frames = size(sub, 3)

            if col_len == 1
                # No column collapsing – just copy the single column per frame
                # `sub` has dimensions (nRows, 1, nFrames)
                for f = 1:n_frames
                    profile_data[:, cur_col] = view(sub, :, 1, f)
                    cur_col += 1
                    next!(prog)
                end
            else
                # Collapse columns by **row‑wise median** (fast, vectorised)
                # `mapslices(median, sub; dims=2)` returns an (nRows, 1, nFrames) array
                collapsed = mapslices(median, sub; dims=2)   # median across cols
                for f = 1:n_frames
                    profile_data[:, cur_col] = view(collapsed, :, 1, f)
                    cur_col += 1
                    next!(prog)
                end
            end
        finally
            close(ds)
        end
    end

    @assert cur_col - 1 == size(profile_data, 2) "Mismatch in counted frames vs. filled columns."
    @info "Finished loading and collapsing $((cur_col-1)) frames into profile_data."
    return nothing
end

"""
    compute_profile_statistics(profile_data::Matrix{Float64}) →
        NamedTuple{(:median, :q25, :q75, :min, :max), Tuple{Vector{Float64},...}}

`profile_data` is (n_rows × n_total_frames).  Each row contains all the
collapsed temperature values for a given height.  The function returns the
per‑row median, the 25‑th and 75‑th quantiles (IQR), and the min / max.
"""
function compute_profile_statistics(profile_data::AbstractArray{Float64})
    # Quantiles are cheap when applied row‑wise via `mapslices`
    median_vec    = mapslices(x -> mean(x), profile_data; dims=2)[:]
    qlow_vec    = mapslices(x -> quantile(x, 0.25), profile_data; dims=2)[:]
    qhigh_vec    = mapslices(x -> quantile(x, 0.75), profile_data; dims=2)[:]
    min_vec    = mapslices(minimum, profile_data; dims=2)[:]
    max_vec    = mapslices(maximum, profile_data; dims=2)[:]

    return (median = median_vec,
            qlow   = qlow_vec,
            qhigh  = qhigh_vec,
            min    = min_vec,
            max    = max_vec)
end

"""
    prepare_multiple_profiles(root_dir::String,
                              t_start::DateTime,
                              t_end::DateTime,
                              specs::Vector{ProfileSpec}) → Vector{NamedTuple}

For each `ProfileSpec` in `specs` the function:
   1. finds the matching NetCDF files,
   2. counts how many frames belong to the interval,
   3. allocates a reduced matrix,
   4. fills it with the collapsed data,
   5. computes median/IQR/min/max per row.

The returned vector has the same order as `specs`; each element is a
`NamedTuple` with fields `median, q25, q75, min, max, row_range, label`.
"""
function prepare_multiple_profiles(files::Vector{String},
                                   t_start::DateTime,
                                   t_end::DateTime,
                                   specs::Vector{ProfileSpec},
                                   frames_per_file::Int,
                                   sample_rate::Float64)

    secs_per_frame = 1.0 / sample_rate
    secs_per_file   = frames_per_file * secs_per_frame   # ≈ 16.666 s

    # Count total frames (shared by all profiles)
    n_frames = count_matching_frames(files, t_start, t_end, secs_per_file, frames_per_file)

    # Allocate a *single* buffer that will be reused for every profile.
    n_rows_max = maximum(length.(getfield.(specs, :row_range)))   # biggest row span
    profile_buffer = Matrix{Float64}(undef, n_rows_max, n_frames)

    results = NamedTuple[]   # will hold one NamedTuple per spec

    @info "Loading $(size(specs,1)) profiles"

    for spec in specs
        # -----------------------------------------------------------------
        # Resize the buffer to the exact number of rows needed for this spec
        # (the buffer is larger than necessary for some specs – that’s fine).
        # -----------------------------------------------------------------
        rows = spec.row_range
        n_rows = length(rows)
        view_buf = @view profile_buffer[1:n_rows, :]   # cheap view, no copy

        # -----------------------------------------------------------------
        # Fill the view with the collapsed data for this spec
        # -----------------------------------------------------------------
        read_and_collapse!(view_buf, files, t_start, t_end,
                           spec.row_range, spec.col_range)

        # -----------------------------------------------------------------
        # Compute statistics on the *filled* view
        # -----------------------------------------------------------------
        stats = compute_profile_statistics(view_buf)

        # -----------------------------------------------------------------
        # Pack everything we need for plotting
        # -----------------------------------------------------------------
        push!(results, (median = stats.median,
                        qlow    = stats.qlow,
                        qhigh    = stats.qhigh,
                        min    = stats.min,
                        max    = stats.max,
                        row_range = rows,
                        label  = spec.label))
    end

    return results
end

function plot_multi_profiles(first_frame::Matrix{Float64},
                             results::Vector{NamedTuple},
                             specs::Vector{ProfileSpec},
                             t_start::DateTime,
                             t_end::DateTime,
                             custom_title::String = "";
                             tmin::Real = -2.0,
                             tmax::Real = 3.0)
    # --------------------------------------------------------------
    # Create the figure with two sub‑plots (profile | reference frame)
    # --------------------------------------------------------------
    fig, (ax_prof, ax_ref) = PyPlot.subplots(1, 2; figsize=(13, 6))

    # --------------------------------------------------------------
    # Colour handling – use Matplotlib's default cycle unless you supply yours
    # --------------------------------------------------------------
    colors = PyPlot.get_cmap("tab10")

    # --------------------------------------------------------------
    # 1️⃣  Plot each profile (left panel)
    # --------------------------------------------------------------
    for (i, res) in enumerate(results)
        col = colors(i-1)

        # y‑axis = pixel rows (the same for every profile – they may differ,
        # but Matplotlib will align them automatically because we plot against the
        # *same* row numbers stored in `res.row_range`.)
        y = collect(res.row_range)

        # Median line
        ax_prof.plot(res.median, y; color=col, lw=2, label="$(res.label) – mean")

        # IQR shaded region
        ax_prof.fill_betweenx(y, res.qlow, res.qhigh; color=col, alpha=0.25, label="$(res.label) – IQR")

        # Min / Max dashed lines
        #ax_prof.plot(res.min, y; color=col, ls="--", lw=1, label="$(res.label) – min/max")
        #ax_prof.plot(res.max, y; color=col, ls="--", lw=1)   # same label not needed
    end

    ax_prof.invert_yaxis()
    ax_prof.set_xlabel("Temperature (°C)")
    ax_prof.set_ylabel("Pixel row (height)")
    ax_prof.set_title("Vertical temperature profiles")
    ax_prof.legend(fontsize=9, framealpha=0.8)
    ax_prof.grid(true)

    # --------------------------------------------------------------
    # 2️⃣  Plot the reference frame (right panel) and draw coloured rectangles
    # --------------------------------------------------------------
    im = ax_ref.imshow(first_frame; cmap=cramericm.batlow, origin="upper", vmin=tmin, vmax=tmax)
    ax_ref.set_title("Sample frame with profile locations")

    # One rectangle per profile – colour matches the curve
    for (i, spec) in enumerate(specs)
        col = colors(i-1)

        # Rectangle geometry (Matplotlib expects lower‑left corner)
        rect = patches.Rectangle((first(spec.col_range) - 0.5, first(spec.row_range) - 0.5),
                length(spec.col_range), length(spec.row_range), linewidth = 2,
                edgecolor = col, facecolor = "none", zorder = 10, alpha = 0.4)
        ax_ref.add_patch(rect)
    end

    # Colour‑bar for the thermal image
    fig.colorbar(im, ax=ax_ref, label="Temperature (°C)")

    # --------------------------------------------------------------
    # 3️⃣  Overall title (interval + optional custom part)
    # --------------------------------------------------------------
    interval_str = @sprintf("%s – %s",
                            Dates.format(t_start, "yyyy‑mm‑dd HH:MM:SS"),
                            Dates.format(t_end,   "yyyy‑mm‑dd HH:MM:SS"))
    suptitle = "IR Temperature Profiles  [$interval_str]"
    if !isempty(custom_title)
        suptitle *= "  –  $custom_title"
    end
    fig.suptitle(suptitle, fontsize=12)

    PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig
end

end #module