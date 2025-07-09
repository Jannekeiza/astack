#!/usr/bin/env python
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------------#
# - Import the useful packages ------------------------------------------------------#

import obspy
import numpy as np
import matplotlib.pyplot as plt
import h5py
import obspyh5
import os
import datetime
import sys

from obspy import UTCDateTime

from obspy import Stream, Trace, read, read_inventory
from obspy.signal.filter import bandpass

from obspy.core import read
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy.signal.trigger import z_detect
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import carl_sta_trig
from obspy.signal.trigger import delayed_sta_lta
from obspy.signal.trigger import pk_baer
from obspy.signal.trigger import ar_pick

from obspy.signal.trigger import coincidence_trigger
from obspy.signal.cross_correlation import xcorr_pick_correction

#from mpl_toolkits.basemap import Basemap

#------------------------------------------------------------------------------------#
# - Set values ----------------------------------------------------------------------#

#Defined from input values
fmin = float(sys.argv[1])
fmax = float(sys.argv[2])
testdir = sys.argv[3]

overwrite = False
if len(sys.argv) > 4:
    if sys.argv[4] == 'True':
        overwrite = True
    elif sys.argv[4] == 'False':
        overwrite = False
    else:
        print("Invalid argument for overwrite. Use True or False.")
        sys.exit(1)

st=Stream()

sampling_rate='H'
channel='Z'

sec=60
lsec=30
usec=30

snr_threshold = 2.5

sample_rate = 0.05 #In seconds

phase_type = "P"

st=Stream()

savedir='/projects/prjs1435/test_waveforms/Figures/P_arrival_plots/'
ev_writedir='/projects/prjs1435/test_waveforms/Astack/'+testdir+'/Input_data/'
maindir='/projects/prjs1435/test_waveforms/seismograms_'+sampling_rate+'H_resp'

eventdir = maindir+'/2020/01/'
# Get the list of events from the directory
events = [d for d in os.listdir(eventdir) if os.path.isdir(os.path.join(eventdir, d))]
events=sorted(events)
print(events)

errorfile=open(os.path.join(ev_writedir, 'error_log.txt'),'w')

#------------------------------------------------------------------------------------#
# - Subroutines ---------------------------------------------------------------------#

def calculate_snr(tr,taupy_time):
    noise_window = [taupy_time - 30, taupy_time - 5]
    signal_window = [taupy_time - 2, taupy_time + 10]

    noise_data = tr.slice(starttime=noise_window[0], endtime=noise_window[1]).data
    signal_data = tr.slice(starttime=signal_window[0], endtime=signal_window[1]).data

    # Compute RMS (Root Mean Square)
    rms_noise = np.sqrt(np.mean(noise_data**2))
    rms_signal = np.sqrt(np.mean(signal_data**2))

    # Calculate SNR
    snr = rms_signal / rms_noise
    #print(f"SNR (RMS): {snr:.2f}")
    
    tr.stats.snr=f"{snr:.2f}"

    return snr

def write_event_file(event, station_count, evlon, evlat, evdep,evortime,ds, phase_type, ev_writedir):
    event_time = datetime.datetime.fromisoformat(evortime[:-1])

    os.makedirs(ev_writedir, exist_ok=True)
    filename = os.path.join(ev_writedir, f"{str(event)}_{str(fmin)}-{str(fmax)}Hz.aq") #
    with open(filename, 'w') as f:
        f.write(f"{station_count}\n") #Number of stations (initially written as 0, then corrected)
        f.write(f"{evlat:.4f} {evlon:.4f} {evdep:.2f}\n") #Event location
        f.write(f"{event_time.year} {event_time.month} {event_time.day}\n") #Event time
        f.write(f"{event_time.hour} {event_time.minute} {event_time.second} 0\n") #Event time + trace start time (I've set to 0... doesn't seem to matter?)
        f.write(f"{ds} {phase_type}\n") #Sample rate and Phase being stacked

    filename_snr = os.path.join(ev_writedir, f"SNR_{str(event)}_{str(fmin)}-{str(fmax)}Hz.txt") #
    with open(filename_snr, 'w') as f:
        for trace in st:
            f.write(f"{trace.stats.station} SNR: {trace.stats.snr}\n")
        
def write_trace_data(st,ev_writedir,event):
    for trace in st:
        station=trace.stats.station
        filename = os.path.join(ev_writedir, f"{str(event)}_{str(fmin)}-{str(fmax)}Hz.aq")
        with open(filename, 'a') as f:
            f.write(f"1 {len(trace.data)} {0} {station}\n") #Trace exist/ have data = 1, number of samples, realtive shift (stacked agasint ak135, not another station ergo = 0. Can change if wanetd though to stack differently), station name
            for sample_value in trace.data:
                f.write(f"{sample_value:.5f} ") #write the trace
            f.write("\n")
    print(filename)

#------------------------------------------------------------------------------------#
# - MAIN ----------------------------------------------------------------------------#

for event in events:
    # check if path exists
    aqpath = os.path.join(ev_writedir, f"{str(event)}_{str(fmin)}-{str(fmax)}Hz.aq")
    if os.path.exists(aqpath):
        print(f"Event {event} already processed")
        if overwrite == True:
            print(f"Overwriting event {event}")
        else:
            print(f"Skipping event {event}")
            continue
    else:
        print(f"Processing event {event}")

    st=Stream()
    for station in ['NE301','NE302','NE303','NE304','NE305','NE306','NE307','NE308','NE310','NE311','NE312','NE317','NE318']:
        file=maindir+'/2020/01/'+event+'/NR/'+station+'/'+station+'.'+event+'.hdf5'
        if os.path.exists(file):
            print(file)
        else:
            print('File not found')
            continue
        
        with h5py.File(file, "r") as f:
         itemname = "/Waveforms"
         if itemname not in f:
          print(f"{file} doesnt have an item /Waveforms")
          errorfile.writelines(f"{file} doesnt have an item /Waveforms \n")
           
         else:

          item=f["/Waveforms"]

          for item2 in item:
            #print(item2)
            #print(item[item2].keys())
            for key in item[item2].keys():
              #print(key)
              dataset=item[item2][key]
              if dataset.attrs['channel'] == 'D'+channel:
                waveform_data = dataset[:]
                start_time = UTCDateTime(dataset.attrs['starttime'])
                trace = Trace(data=waveform_data)
                trace.stats.network = dataset.attrs['network']
                trace.stats.station = dataset.attrs['station']
                trace.stats.location = dataset.attrs['location']
                trace.stats.channel = dataset.attrs['channel']
                trace.stats.starttime = start_time
                trace.stats.distance = f.attrs['distance']
                trace.stats.latitude = f.attrs['station latitude']
                trace.stats.longitude = f.attrs['station longitude']
                trace.stats.event_latitude = f.attrs['event latitude']
                trace.stats.event_longitude = f.attrs['event longitude']
                trace.stats.event_depth=f.attrs['event depth']
                trace.stats.ev_ortime=f.attrs['event origin time']

                trace.stats.starttime = start_time
                if sampling_rate == 'L':
                    trace.stats.sampling_rate = 1
                elif sampling_rate == 'H':
                    trace.stats.sampling_rate = 100

                starttime=start_time+(20)*60-lsec
                endtime=start_time+(20)*60+usec
            
                trace.trim(starttime=starttime, endtime=endtime)
                trace.detrend('linear')
                trace.taper(max_percentage=0.05, type='cosine')

                trace.filter("bandpass", freqmin=fmin, freqmax=fmax, zerophase=True)
                sample_seconds = 1 / sample_rate
                trace.resample(sample_seconds)
                trace.data=trace.data*(10.**9)
                taupy_time = starttime + lsec

                snr=calculate_snr(trace,taupy_time)

                if snr > snr_threshold:
                    print(f"* SNR = {snr:.2f}, SNR passed threshold - saved")
                    st += trace
                else:
                    print(f"* SNR = {snr:.2f}, SNR didn't pass threshold - skipped")  
    f.close()

    #------------------------------------------------------------------------------------#
    # - Write traces to file ------------------------------------------------------------#

    station_count=len(st)
    evlon=trace.stats.event_longitude
    evlat=trace.stats.event_latitude
    evdep=trace.stats.event_depth
    evortime=trace.stats.ev_ortime
    df = trace.stats.sampling_rate
    ds=1/df

    print(st)
    
    if len(st) > 2:
        ev_writedir = "/projects/prjs1435/test_waveforms/Astack/"+testdir+"/Input_data"
        print(event, "passes SNR for enough stations, nr stations = ",len(st))
        write_event_file(event, station_count, evlon, evlat, evdep,evortime,sample_rate, phase_type, ev_writedir)
        write_trace_data(st,ev_writedir,event)

    elif len(st) > 0 and len(st) < 2:
        print(event, "doesn't pass SNR for enough stations, nr stations = ",len(st))
        ev_writedir = "/projects/prjs1435/test_waveforms/Astack/"+testdir+"/Input_data/Unused_data"

        write_event_file(event, station_count, evlon, evlat, evdep,evortime,sample_rate, phase_type, ev_writedir)
        write_trace_data(st,ev_writedir,event)
    
    else:
        print(event, "doesn't pass SNR for any station")
    
    # flush print statements
    sys.stdout.flush()
    
    

