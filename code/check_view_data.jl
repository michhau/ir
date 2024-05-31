######################################################
###           CHECK & VIEW IR-SCREEN DATA          ###
###            author: Michi Haugeneder            ###
######################################################
#=
Checking and viewing of the already imported data (.nc-files)
=#
using PyCall, Statistics, LaTeXStrings, Dates #DelimitedFiles, Dates, Random, FFTW, Dates
#StatsBase,ImageFiltering,
import PyPlot
GridSpec = pyimport("matplotlib.gridspec")
mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")

importdir = joinpath(@__DIR__)
pathtofile = "/home/haugened/Documents/data/ir/"

include(joinpath(importdir, "src", "ir_evaluation.jl"))
#include(joinpath(importdir, "src", "general.jl"))
import .irev
#import .gen
PyPlot.pygui(true)

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################
#completing the filepath
pathtofile = joinpath(pathtofile, "210531_143457/")
#which sequence should be loaded (string)
ncfile = "irdata_09_to_11.nc"
#for plotting
cmperpxl = 0.6
cmperpxlu = 0.606
cmperpxlw = 0.549
framespersec = 29.9
######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###                   IMPORT STUFF                 ###
######################################################
IRdata = irev.loadfromNetCDF4(joinpath(pathtofile, ncfile), "irdata")

irev.showheatmaptrad(IRdata[:, :, 3250], 4, 8, "")

#=
#median for all pixels in a period
period = 3241:3270

irmedian = fill(NaN, size(IRdata, 1), size(IRdata, 2))
irvar = fill(NaN, size(IRdata, 1), size(IRdata, 2))
for i in 1:size(IRdata, 2)
    for j in 1:size(IRdata, 1)
        irmedian[j, i] = median(IRdata[j, i, period])
        irvar[j, i] = std(IRdata[j, i, period])
    end
end

irev.showheatmap(irmedian, 4, 8, "")
=#
######################################################
###                    PLOT DATA                   ###
######################################################
#animation of IRdata
starttime = Time(14,40,57)

"Update-function for the animation"
function animupdate(frame::Int64)
    fig.suptitle(string("Duerrboden 31.05.2021 ", starttime + Millisecond(round(Int, 1000*frame/framespersec))," LT"))
    hm.set_data(IRdata[:, :, frame])
    return hm
end

xis0pxl = 73
fig = PyPlot.figure(figsize=(16, 9))
ax = fig.subplots(1, 1)
fig.suptitle(string("Duerrboden 31.05.2021 ", starttime," LT"))
xleft = -(xis0pxl - 1) * cmperpxlu / 100
xright = size(IRdata, 2) * cmperpxlu / 100 + xleft
ybottom = 0
ytop = size(IRdata, 1) * cmperpxlw / 100
hm = ax.imshow(IRdata[:, :, 1], extent=[xleft, xright, ybottom, ytop], cmap="turbo", vmin=7.5, vmax=11)
ax.set_xlabel("fetch distance over snow x [m]")
ax.set_ylabel("height [m]")
hm_cb = fig.colorbar(hm, ax=ax)
hm_cb.set_label(L"T~\mathrm{[^\circ C]}")
ani = animation.FuncAnimation(fig, animupdate, frames=collect(1:size(IRdata, 3)), interval=5, repeat_delay=500)
#uncomment following line to save animation
#ani.save(string(pathtofile, "video_", ncfile, ".mp4"), fps=framespersec)

######################################################
#interactive plotting

#function to update heatmap color minimum
function submit_vmin(value)
    hm.set_clim(vmin=value)
    PyPlot.draw()
end

#function to update heatmap color maximum
function submit_vmax(value)
    hm.set_clim(vmax=value)
    PyPlot.draw()
end

#updatefunction for the sliders
function update(val)
    global framei
    framei = round(Int64, slider_time.val)
    hm.set_data(IRdata[:, :, framei])
    fig.suptitle(string("210215_ref - frame ", framei))
    fig.canvas.draw_idle()
end

function onclick(event)
    global framei
    if event.key == "left" && framei != 1
        framei = framei - 1
    elseif event.key == "right" && framei != size(IRdata, 3)
        framei = framei + 1
    end
    hm.set_data(IRdata[:, :, framei])
    slider_time.set_val(framei)
    fig.canvas.draw_idle()
end

framei = 1
colormin = -0.5
colormax = 0.5
fig = PyPlot.figure(figsize=(16, 9))
gs = GridSpec.GridSpec(1, 2, width_ratios=[3, 1], left=0.05)
ax1 = fig.add_subplot(get(gs, 0))
fig.suptitle("210215_ref - frame 1")
ax1.set_xlabel("spatial axis - column [pxl]")
ax1.set_ylabel("spatial axis - row [pxl]")
#ax1.yaxis.set_ticklabels([])
hm = ax1.imshow(IRdata[:, :, 1], cmap="turbo", vmin=colormin, vmax=colormax)
#ax2 = fig.add_subplot(get(gs,1))
axslider_time = PyPlot.axes([0.67, 0.85, 0.25, 0.01])
slider_time = mpwidgets.Slider(axslider_time, "frame [pics]", 1, size(IRdata, 3), valinit=1, valstep=1, valfmt="%1.0f")

axbox_vmin = PyPlot.axes([0.67, 0.78, 0.03, 0.02])
axbox_vmax = PyPlot.axes([0.75, 0.78, 0.03, 0.02])
text_box_vmin = mpwidgets.TextBox(axbox_vmin, "vmin")
text_box_vmax = mpwidgets.TextBox(axbox_vmax, "vmax")
text_box_vmin.on_submit(submit_vmin)
text_box_vmax.on_submit(submit_vmax)
text_box_vmin.set_val(colormin)
text_box_vmax.set_val(colormax)
cb = fig.colorbar(hm, ax=ax1, orientation="horizontal", fraction=0.03)
cb.set_label("temperature [C]")

cid = fig.canvas.mpl_connect("key_press_event", onclick)

slider_time.on_changed(update)

PyPlot.show()
