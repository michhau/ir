#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conversion from IRB files to NetCDF files.

@author: Marcus Huntemann
"""

IMFIELD_SIZE=2*1024*768 #bytes*imagex*imagey

from os import listdir, path
import pandas as pd
import numpy as np
from glob import glob
import re
from netCDF4 import Dataset
from datetime import datetime


def read_cam_data(fn,read_ascii=False):
    """
    reads in a single "irb" file and stores the content in a struct
    
    Parameters
    ----------
    fn : irb filename 
    read_ascii : bool, optional
        if True it will read the image parameter from corresponding csv files 
        assming that all frames has been exported by IRBIS 3.1

    Returns
    -------
    outdict : Dictionary containing all parameters in the irb file

    """
    with open(fn,"rb") as o:
        data=o.read()
    
    img_idxs=[(m.start(0)-IMFIELD_SIZE, m.start(0)) for m in re.finditer(b'\\[GPS\\]', data)]

    lats=np.array([float(i)/60 for i in re.findall(b"Latitude=([-0-9.]{,20})",data)])
    lons=np.array([float(i)/60 for i in re.findall(b"Longitude=([-0-9.]{,20})",data)])
    alti=np.array([float(i) for i in re.findall(b"Altitude=([-0-9.]{,20})",data)])
    speed=np.array([float(i)*0.514444 for i in re.findall(b"Speed_Knot=([-0-9.]{,20})",data)])
    course=np.array([float(i) for i in re.findall(b"Course=([-0-9.]{,20})",data)])
    hdop=np.array([float(i) for i in re.findall(b"HDOP=([-0-9.]{,20})",data)])

    #utctimes=[i.decode() for i in re.findall(b"UTCTime=([0-9.]{,20})",data)]
    utcdates=[i.decode() for i in re.findall(b"UTCDate=([0-9.]{,20})",data)]
    utctimes=[]
    for i in re.findall(b"UTCTime=([0-9.]{,20})",data):
        if i[-6:]==b"60.000":
            utctimes.append((i[:-6]+b"59.999").decode())
        else:
            utctimes.append(i.decode())
            
    imgs=[]
    dates=[]
    if read_ascii:
        afns=sorted(glob(fn[:-4]+"*.csv"))

    for i, idx in enumerate(img_idxs):
        if read_ascii:
            imgs.append(pd.read_csv(afns[i],skiprows=23,delimiter=";",header=None).values[:,:-1])
        else:
            imgs.append(np.frombuffer(data[idx[0]:idx[1]],dtype=np.uint16).reshape(768,1024))
            #imgs.append((np.random.rand(480,640)*10+np.arange(640).T/640/10+270.0)) #dummy data for testing
        
        utctimes[i]="0"*(6-len(utctimes[i].split(".")[0]))+utctimes[i]
        cdate=datetime.strptime(utcdates[i]+':'+utctimes[i],"%d%m%y:%H%M%S.%f") #seperator eingefügt, da vorher komische Zeitstempel
        dates.append(cdate)
        
#     import IPython
#     IPython.embed()
    
    startdate=str(datetime.fromordinal(dates[0].toordinal()))

    filename=path.dirname(fn)+"/IRdata_PS149-CONTRASTS_"+startdate[0:10]+".nc"

    imgs=np.array(imgs)#.swapaxes(0,2)
    times=[dates[i]-datetime.fromordinal(dates[0].toordinal()) for i in range(len(img_idxs))]
    times=[t.seconds+t.days*60*60*24 for t in times]
    outdict={}
    outdict["ncfilename"]=filename #path.dirname(fn)+"_data.nc"
    outdict["infilename"]=fn
    outdict["startdate"]=startdate
    outdict["TBs"]=imgs
    outdict["lats"]=lats
    outdict["lons"]=lons
    outdict["altitude"]=alti
    outdict["speed"]=speed
    outdict["course"]=course
    outdict["hdop"]=hdop
    outdict["times"]=times
    return outdict

def prepare_netcdf(od):
    """
    prepares a skeleton of a netcdf4 file to store the irb files in
    
    Parameters
    ----------
    od : Dictory as returned by the `read_cam_data` routine

    Returns
    -------
    nothing

    """
    with Dataset(od["ncfilename"], "w", format="NETCDF4") as rootgrp:
        rootgrp.title="UB IR Camera data"
        rootgrp.institution="University of Bremen"
        rootgrp.product_name="IR Camera brightness temperatures dump"
        rootgrp.PI_name="Gunnar Spreen"
        rootgrp.conventions = "CF-1.7"
        rootgrp.createDimension('X', 640)
        rootgrp.createDimension('Y', 480)
        rootgrp.createDimension('time',None)#subject to change
        #time
        time=rootgrp.createVariable("time","f4",("time",),zlib=True)
        time.units = "s"
        time.stardard_name = "time"
        time.long_name="seconds since "+od["startdate"]
        
        #latitude
        latitude=rootgrp.createVariable("latitude","f4",("time",),zlib=True)
        latitude.units = "degree_north"
        latitude.stardard_name = "latitude"
        #longitude
        longitude=rootgrp.createVariable("longitude","f4",("time",),zlib=True)
        longitude.units = "degree_east"
        longitude.stardard_name = "longitude"
        #altitude
        altitude=rootgrp.createVariable("altitude","f4",("time",),zlib=True)
        altitude.units = "m"
        altitude.standard_name = "altitude"
        altitude.long_name = "height above the geoid from GPS"
        #speed
        speed=rootgrp.createVariable("speed","f4",("time",),zlib=True)
        speed.units = "m/s"
        speed.standard_name = "platform_speed"
        speed.long_name = "platform speed as reported by GPS"
        
        #course
        course=rootgrp.createVariable("course","f4",("time",),zlib=True)
        course.units = "degree"
        course.standard_name = "platform_azimuth_angle"
        #hdop
        hdop=rootgrp.createVariable("hdop","f4",("time",),zlib=True)
        hdop.units = "1"
        hdop.standard_name = "hdop"
        hdop.long_name="horizontal dilution of precision"
        #TBs
        TBs=rootgrp.createVariable("temps","f4",("time","Y","X"),zlib=True,least_significant_digit=2)
        TBs.standard_name = "brightness_temerature"
        TBs.units = "K"

def write_into_NC(od,offset=0):
    """
    writes the data into a prepared netCDF4 file (created by `prepare_netcdf`)

    Parameters
    ----------
    od : Dictionary as prepared by `read_cam_data`
    offset : int, optional
        starting frame in the dataset

    Returns
    -------
    nothing

    """
    with Dataset(od["ncfilename"], "a") as D:
        D["hdop"][offset:]= od["hdop"][:]
        D["course"][offset:] = od["course"][:]
        D["speed"][offset:] = od["speed"][:]
        D["altitude"][offset:] = od["altitude"][:]
        D["longitude"][offset:] = od["lons"][:]
        D["latitude"][offset:] = od["lats"][:]
        D["time"][offset:]=od["times"][:]
        D["temps"][offset:,:,:] = od["TBs"][:,:,:]

def irb_fold_to_NC(foldername,read_ascii=True,outfile=False):
    """
    converts all data of an IR camera recorded folder into a NetCDF file
    read ascii can be set in order to use the IR camera data from exported
    ascii files instead of the irb files this is required for the compressed
    irb files
    """
    fns=sorted(glob(path.join(foldername,"*.irb")))
    od=read_cam_data(fns[0],read_ascii=read_ascii)
    if outfile:
        od["ncfilename"]=outfile
    prepare_netcdf(od)
    N_FRAMES_PER_FILE=len(od["times"])
    for i,fn in enumerate(fns):
        od=read_cam_data(fn,read_ascii=read_ascii)
        if outfile:
            od["ncfilename"]=outfile
        write_into_NC(od,i*N_FRAMES_PER_FILE)
        
'''        
if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser(description=
    "stores data from IRB files in to a NetCDF file")

    p.add_argument("irbfolder", type=str,help="location of the irb files to store in one NC file, should stem from one time series or flight")
    p.add_argument("--csv",help="uses the IR camera data from exported CSV files (stored in irbfolder) instead of from the IRB files",action='store_true')
    p.add_argument("--outfile", type=str, help="name of the NetCDF outfile, defualt: 'irbfolder'_data.nc at the location of 'irbfolder'") 

    pargs=p.parse_args()
    if not path.exists(pargs.irbfolder):
        raise FileNotFoundError("\ndirectory "+pargs.irbfolder+" not found")
    print (pargs)
    if pargs.outfile:
        irb_fold_to_NC(pargs.irbfolder,read_ascii=pargs.csv,outfile=pargs.outfile)
    else:
        irb_fold_to_NC(pargs.irbfolder,read_ascii=pargs.csv)
'''

if __name__ == '__main__':
    irb_fold_to_NC('/home/michi/Documents/slf/CONTRASTS25/data/ir/tmp_3d/250826_111801/')