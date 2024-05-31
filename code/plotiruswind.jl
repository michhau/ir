######################################################
###          PLOT 2D-US WIND (IR-SCREEN)           ###
###            author: Michi Haugeneder            ###
######################################################
#=
Plot 2D ultrasonic wind-measurement for given source file and time interval
For use in terminal (argument parsing)
arguments:
 1. sourcefile
 2. starttime (yyyymmddhhmmss)
 3. endtime (yyymmddhhmmss)
=#
using PyCall, ArgParse, Dates
#StatsBase,ImageFiltering,
import PyPlot
importdir = @__DIR__
pathtofile = "/home/haugened/Documents/data/2d_us_ir/"
include(joinpath(importdir, "src", "ir_evaluation.jl"))
import .irev

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################

######################################################

######################################################
###                ARGUMENT PARSING                ###
######################################################
function parse_terminal()
    s = ArgParseSettings()
    @add_arg_table s begin
        "sourcefile"
        help = "path to sourcefile containing measured data"
        arg_type = String
        required = true
        "starttime"
        help = "starttime (yyyymmddhhmmss)"
        required = true
        "endtime"
        help = "endtime (yyyymmddhhmmss)"
        required = true
    end

    return parse_args(s)
end

######################################################

######################################################
###                       MAIN                     ###
######################################################
parsed_args = parse_terminal()

#converting the parsed arguments to desired filetypes
source = parsed_args["sourcefile"]
dateformat = DateFormat("yyyymmddHHMMSS")
starttime = DateTime(parsed_args["starttime"], dateformat)
endtime = DateTime(parsed_args["endtime"], dateformat)

println()
println("-----------S-T-A-R-T-------------")

#load data
(timeaxis, winddata,) = irev.fromIRloadwinddata(source, starttime, endtime)


######################################################
###                      PLOT                      ###
######################################################
colorax1 = "#1f77b4"
colorax2 = "#ff7f0e"
PyPlot.pygui(true)
fig = PyPlot.figure(figsize=(12, 7))
ax1 = fig.add_subplot(111)
ax1.set_title("2D-ultrasonic wind measurement @ IR screen")
ax1.set_xlabel("time")
ax1.set_ylabel("speed [m/s]", color=colorax1)
ax1.tick_params(axis="y", labelcolor=colorax1)
ax1.grid()
ax1.plot(timeaxis, winddata[:, 2], color=colorax1)
ax2 = ax1.twinx()
ax2.set_ylabel("direction [deg]", color=colorax2)
ax2.tick_params(axis="y", labelcolor=colorax2)
ax2.plot(timeaxis, winddata[:, 1], color=colorax2, alpha=0.7)
fig.tight_layout()
PyPlot.show()
######################################################

######################################################


######################################################

println("------------D-O-N-E---------------")
