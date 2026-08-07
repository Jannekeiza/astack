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
sample_rate = sys.argv[3]
base_dir = sys.argv[4]
waveforms_root = sys.argv[5]

overwrite = False
if len(sys.argv) > 6:
    if sys.argv[6] == 'True':
        overwrite = True
    elif sys.argv[6] == 'False':
        overwrite = False
    else:
        print("Invalid argument for overwrite. Use True or False.")
        sys.exit(1)

# Which source trees to search, e.g. "FDSN,Groningen" to skip Dictum.
# Defaults to all three. Names are matched case-insensitively.
included_sources = ['FDSN', 'Dictum', 'Groningen']
if len(sys.argv) > 7 and sys.argv[7].strip():
    included_sources = [name.strip() for name in sys.argv[7].split(',')]

# The three source trees to search for each event. Each has its own network
# code and station list (same event IDs are shared across all three trees,
# so matching traces found in more than one source get merged into a single
# combined stream/output file per event).
all_sources = [
    {
        'name': 'FDSN',
        'dir': os.path.join(waveforms_root, 'FDSN', 'seismograms'),
        'network': 'NL',
        'stations': ['ARCN','DBN','G81B','G82B','G83B','G84B','HGN','HRKB','MAME','NE05','NE427','OPLO','TERZ','VKB','WTSB'],
    },
    {
        'name': 'Dictum',
        'dir': os.path.join(waveforms_root, 'Dictum', 'seismograms_20Hz'),
        'network': 'NR',
        'stations': ['NE400','NE401','NE403','NE405','NE406','NE407','NE408','NE409','NE410','NE411','NE412','NE413','NE414','NE416','NE418','NE420','NE421','NE422','NE424','NE425','NE426','NE427'],
    },
    {
        'name': 'Groningen',
        'dir': os.path.join(waveforms_root, 'Groningen', 'seismograms_USGS_HH'),
        'network': 'NR',
        'stations': ['NE301','NE302','NE303','NE304','NE305','NE306','NE307','NE308','NE309','NE310','NE311','NE312','NE317','NE318'],
    },
]

valid_names = {s['name'].lower() for s in all_sources}
unknown = [name for name in included_sources if name.lower() not in valid_names]
if unknown:
    print(f"Unknown source(s) {unknown}. Valid options: {[s['name'] for s in all_sources]}")
    sys.exit(1)

included_lower = {name.lower() for name in included_sources}
sources = [s for s in all_sources if s['name'].lower() in included_lower]
if not sources:
    print("No sources selected - nothing to do.")
    sys.exit(1)

print(f"Using sources: {[s['name'] for s in sources]}")

st=Stream()

if sample_rate == '20':
    sample_type = 'H'
elif sample_rate == '1':
    sample_type = 'L'

channel='Z'

sec=60
lsec=30
usec=30

snr_threshold = 2.5

phase_type = "P"

st=Stream()

ev_writedir=base_dir+'/Input_data/'

errorfile=open(os.path.join(ev_writedir, 'error_log.txt'),'w')

#------------------------------------------------------------------------------------#
# - Subroutines ---------------------------------------------------------------------#

def calculate_snr(tr,taupy_time):
    noise_window = [taupy_time - 30, taupy_time - 5]
    signal_window = [taupy_time - 2, taupy_time + 7]

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
    
    if not os.path.exists(ev_writedir):
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

for year in ['2022','2023','2024','2025']:
    for mon in ['01','02','03','04','05','06','07','08','09','10','11','12']:
        # Gather the union of events found across all three source trees for this year/month
        events = set()
        for source in sources:
            eventdir = os.path.join(source['dir'], year, mon)
            if os.path.exists(eventdir):
                events.update(d for d in os.listdir(eventdir) if os.path.isdir(os.path.join(eventdir, d)))
        events = sorted(events)

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
            for source in sources:
                network = source['network']
                evdir = os.path.join(source['dir'], year, mon, event)

                for station in source['stations']:
                    stdir = os.path.join(evdir, network, station)
                    file = os.path.join(stdir, f"{station}.{event}.hdf5")
                    if os.path.exists(file):
                        print(f"Processing file: {file}")
                    else:
                        print('File not found')
                        #remove the directory if it's empty
                        if os.path.exists(os.path.join(evdir, network)):
                            if os.path.exists(os.path.join(evdir, network, station)):
                                if len(os.listdir(stdir)) == 0:
                                    os.rmdir(stdir)
                                    print(f"Removed empty station directory: {stdir}")
                            else:
                                print(f"Directory {stdir} does not exist")
                                if len(os.listdir(os.path.join(evdir, network))) == 0:
                                    os.rmdir(os.path.join(evdir, network))
                                    os.rmdir(evdir)
                                    print(f"Removed empty event directory: {evdir}")
                        else:
                            print(f"Directory {evdir} does not exist")

                        continue

                    try:
                        with h5py.File(file, "r") as f:
                            itemname = "/Waveforms"
                            if itemname not in f:
                                print(f"{file} doesnt have an item /Waveforms")
                                errorfile.writelines(f"{file} doesnt have an item /Waveforms \n")
                        
                            else:

                                item=f["/Waveforms"]

                                for item2 in item:
                                  print(item2)
                                  for item3 in item[item2]:
                                       name3=item3
                                       print(name3)
                                    #if name3 ==  'waveforms':
                                       #for item3 in item[item2]:
                                       print(item3)
                                       print(item[item2][item3].keys())
                                       #print(item[item2].keys())
                                       for key in item[item2][item3].keys():
                                        print(key)
                                        dataset=item[item2][item3][key]
                                        #print(dataset)
                                        print(dataset.attrs['channel'])
                                        if 'channel' in dataset.attrs and dataset.attrs['channel'] == 'D' + channel:
                                            print("Processing HHZ channel for station ", dataset.attrs['station'])
                                            
                                            #if dataset.attrs['channel'] == 'D'+channel:
                                            waveform_data = dataset[:]
                                            #print('waveform data=' + str(waveform_data))
                                            start_time = UTCDateTime(dataset.attrs['starttime'])
                                            #end_time = UTCDateTime(dataset.attrs['endtime'])
                                            trace = Trace(data=waveform_data)
                                            print('Trace: ', trace)

                                            trace.stats.network = dataset.attrs['network']
                                            trace.stats.station = dataset.attrs['station']
                                            trace.stats.location = dataset.attrs['location']
                                            trace.stats.channel = dataset.attrs['channel']
                                            trace.stats.distance = f.attrs['distance']
                                            trace.stats.latitude = f.attrs['station latitude']
                                            trace.stats.longitude = f.attrs['station longitude']
                                            trace.stats.event_latitude = f.attrs['event latitude']
                                            trace.stats.event_longitude = f.attrs['event longitude']
                                            trace.stats.event_depth=f.attrs['event depth']
                                            trace.stats.ev_ortime=f.attrs['event origin time']
                                            #trace.stats.sampling_rate = f.attrs['sampling_rate']

                                            trace.stats.starttime = start_time
                                            #trace.stats.endtime = end_time
                                            print(start_time)

                                            num_samples = dataset.shape[0]  # Access the first (and only) dimension
                                            print(f"Number of samples: {num_samples}")
                                            sampling_rate=(num_samples-1)/(2*60*60)

                                            trace.stats.sampling_rate = sampling_rate

                                            #print(sampling_rate)

                                            #print(trace)
                                            if trace.stats.distance > 11000:
                                                print(f"Event is too far from station {trace.stats.station}, distance = {trace.stats.distance} km")
                                                errorfile.writelines(f"Event {event} is too far from station {trace.stats.station}, distance = {trace.stats.distance} km \n")
                                                continue

                                            starttime=start_time+(20)*60-lsec
                                            endtime=start_time+(20)*60+usec

                                            trace.trim(starttime=starttime, endtime=endtime)
                                            #trace.detrend('linear')
                                            #trace.taper(max_percentage=0.05, type='cosine')

                                            #print(starttime, endtime)
                                            #print(trace)
                                            print(f"Applying bandpass filter: {fmin} - {fmax} Hz")
                                            trace.filter("bandpass", freqmin=fmin, freqmax=fmax, zerophase=True)
                                            #sample_seconds = 1 / sample_rate
                                            # resample to 20 Hz if not already at 20 Hz (to match ak135)
                                            #if sampling_rate < 0.05:
                                            print(f"Resampling to {sample_rate} Hz")
                                            trace.resample(int(sample_rate))
                                            print("trace.stats.sampling_rate after resampling: ", trace.stats.sampling_rate)
                                            num_samples = dataset.shape[0]  # Access the first (and only) dimension
                                            print(f"Nr of samples after resampling: {num_samples}")

                                            trace.data=trace.data*(10.**9)
                                            taupy_time = starttime + lsec

                                            snr=calculate_snr(trace,taupy_time)
                                            
                                            print('trace processed: ', trace)

                                            if snr > snr_threshold:
                                                print(f"* SNR = {snr:.2f}, SNR passed threshold - saved")
                                                st += trace
                                            else:
                                                print(f"* SNR = {snr:.2f}, SNR didn't pass threshold - skipped")

                                            # check if distance of event is within 90 degree radius from station
                                  else:
                                    #print(item[item2].keys())
                                    for key in item[item2].keys():
                                        print(key)
                                        dataset=item[item2][key]
                                        #print(dataset)
                                        #print(dataset.attrs.keys())
                                        if 'channel' in dataset.attrs and dataset.attrs['channel'] == 'D' + channel:
                                            #if dataset.attrs['channel'] == 'D'+channel:
                                            waveform_data = dataset[:]
                                            #print('waveform data=' + str(waveform_data))
                                            start_time = UTCDateTime(dataset.attrs['starttime'])
                                            #end_time = UTCDateTime(dataset.attrs['endtime'])
                                            trace = Trace(data=waveform_data)
                                            print('trace=',trace)

                                            trace.stats.network = dataset.attrs['network']
                                            trace.stats.station = dataset.attrs['station']
                                            trace.stats.location = dataset.attrs['location']
                                            trace.stats.channel = dataset.attrs['channel']
                                            trace.stats.distance = f.attrs['distance']
                                            trace.stats.latitude = f.attrs['station latitude']
                                            trace.stats.longitude = f.attrs['station longitude']
                                            trace.stats.event_latitude = f.attrs['event latitude']
                                            trace.stats.event_longitude = f.attrs['event longitude']
                                            trace.stats.event_depth=f.attrs['event depth']
                                            trace.stats.ev_ortime=f.attrs['event origin time']
                                            #trace.stats.sampling_rate = f.attrs['sampling_rate']

                                            trace.stats.starttime = start_time
                                            #trace.stats.endtime = end_time
                                            print(start_time)

                                            num_samples = dataset.shape[0]  # Access the first (and only) dimension
                                            print(f"Number of samples: {num_samples}")
                                            #sampling_rate=(num_samples-1)/(2*60*60)

                                            #trace.stats.sampling_rate = sampling_rate

                                            #print(sampling_rate)

                                            #print(trace)
                                            if trace.stats.distance > 11000:
                                                print(f"Event is too far from station {trace.stats.station}, distance = {trace.stats.distance} km")
                                                errorfile.writelines(f"Event {event} is too far from station {trace.stats.station}, distance = {trace.stats.distance} km \n")
                                                continue

                                            starttime=start_time+(20)*60-lsec
                                            endtime=start_time+(20)*60+usec

                                            trace.trim(starttime=starttime, endtime=endtime)
                                            #trace.detrend('linear')
                                            #trace.taper(max_percentage=0.05, type='cosine')

                                            #print(starttime, endtime)
                                            #print(trace)
                                            print(f"Applying bandpass filter: {fmin} - {fmax} Hz")
                                            trace.filter("bandpass", freqmin=fmin, freqmax=fmax, zerophase=True)
                                            #sample_seconds = 1 / sample_rate
                                            # resample to 20 Hz if not already at 20 Hz (to match ak135)
                                            #if sampling_rate < 0.05:
                                            print(f"Resampling to {sample_rate} Hz")
                                            trace.resample(int(sample_rate))
                                            print("trace.stats.sampling_rate after resampling: ", trace.stats.sampling_rate)
                                            num_samples = dataset.shape[0]  # Access the first (and only) dimension
                                            print(f"Nr of samples after resampling: {num_samples}")

                                            trace.data=trace.data*(10.**9)
                                            taupy_time = starttime + lsec

                                            snr=calculate_snr(trace,taupy_time)
                                            
                                            print('trace processed: ', trace)

                                            if snr > snr_threshold:
                                                print(f"* SNR = {snr:.2f}, SNR passed threshold - saved")
                                                st += trace
                                            else:
                                                print(f"* SNR = {snr:.2f}, SNR didn't pass threshold - skipped")

                                f.close()
                    except:
                        print(f"Error processing file: {file}")
            #------------------------------------------------------------------------------------#
            # - Write traces to file ------------------------------------------------------------#

            station_count=len(st)
            if station_count == 0:
                print(event, "doesn't pass SNR for any station")
                continue

            evlon=trace.stats.event_longitude
            evlat=trace.stats.event_latitude
            evdep=trace.stats.event_depth
            evortime=trace.stats.ev_ortime
            df = trace.stats.sampling_rate
            ds=1/df

            print('Stream of traces saved: ', st)

            if station_count > 1:
                ev_writedir = base_dir+"/Input_data"
                print(event, "passes SNR for enough stations, nr stations = ",len(st))
                write_event_file(event, station_count, evlon, evlat, evdep,evortime,ds, phase_type, ev_writedir)
                write_trace_data(st,ev_writedir,event)

            elif station_count > 0 and station_count < 2:
                print(event, "doesn't pass SNR for enough stations, nr stations = ",station_count)
                ev_writedir = base_dir+"/Input_data/Unused_data"

                write_event_file(event, station_count, evlon, evlat, evdep,evortime,ds, phase_type, ev_writedir)
                write_trace_data(st,ev_writedir,event)

                
                # flush print statements
                sys.stdout.flush()
                    

                

