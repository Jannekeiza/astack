"""""
mseed_to_tcas.py
Reads in an events.csv file that contains the time, location  and id of an event.
For each station in station_list.csv (contains station name, lat, long, elev) and calculates
the predicted arrival of a desired teleseismic phase and searches an archive to see if 
mseed exists for that event and station. Filters the st and writes to an event .aq file 
needed for tcas. The predicted arrival of the phase is centred on the lower and upper crop 
time... e.g. if lower_crop = 20 and upper_crop=40, then the predicted arrival would be 
at 20 seconds of the stream. 

The run_all.py script will run this parallelised for speed (incorporated in workflow). 
mseed_to_tcas.py is more of an indvual/ check for the  mseed to tcas step. 

-T.O'H 12/24.

"""""

import pandas as pd
import datetime
import os
import glob
import numpy as np
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy import read, read_inventory, UTCDateTime

def main():
    phase_type = "P"
    events_df = read_csv('events.csv')
    stations_df = read_csv('station_list.csv')
    mseed_path = '/raid5/Iceland/data/ARCHIVE/'
    inv = None #read_inventory("/raid4/tpo21/scripts/combined_inventory.xml")
    event_save_path = "./"
    
    sample_rate = 0.020 #In seconds
    
    # List of frequency bands for filtering
    bandpass_freq_min = 0.5
    bandpass_freq_max = 2.0
    
    lower_crop_time = 20
    upper_crop_time = 40
    
    for event_row in events_df.itertuples(index=False):
        process_event(event_row, stations_df, phase_type, sample_rate, mseed_path, lower_crop_time, upper_crop_time, bandpass_freq_min, bandpass_freq_max, event_save_path, inv)

def read_csv(filename):
    return pd.read_csv(filename, header=None)

def write_event_file(event_number, station_count, event_coords, event_time, sample_interval, phase_type, event_save_path):
    # Create the directory if it doesn't exist
    os.makedirs(event_save_path, exist_ok=True)
    filename = os.path.join(event_save_path, f"event{event_number}.aq")
    with open(filename, 'w') as f:
        f.write(f"{station_count}\n") #Number of stations (initially written as 0, then corrected)
        f.write(f"{event_coords[0]:.4f} {event_coords[1]:.4f} {event_coords[2]:.2f}\n") #Event location
        f.write(f"{event_time.year} {event_time.month} {event_time.day}\n") #Event time
        f.write(f"{event_time.hour} {event_time.minute} {event_time.second} 0\n") #Event time + trace start time (I've set to 0... doesn't seem to matter?)
        f.write(f"{sample_interval} {phase_type}\n") #Sample rate and Phase being stacked

#Updates the .aq file station count (first line) once all known stations have been searched       
def update_station_count(event_number, station_count, event_save_path):
    filename = os.path.join(event_save_path, f"event{event_number}.aq")
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines[0] = f"{station_count}\n"
    with open(filename, 'w') as f:
        f.writelines(lines)        

def read_mseed(base_directory,phase_type, station, arrival_time, event_number, sample_rate, upper_crop_time, lower_crop_time, bandpass_freq_min, bandpass_freq_max, event_save_path, inv):
    try:
        if not isinstance(arrival_time, str):
            arrival_time = str(arrival_time)

        #Predicted phase arrival time
        arrival_time = datetime.datetime.strptime(arrival_time, "%Y-%m-%d %H:%M:%S.%f")

        #Longer window for detrending etc.
        readstarttime = UTCDateTime(arrival_time) - 3600
        readendtime = UTCDateTime(arrival_time) + 3600
        #Wanted stack time. THe stack will be centred on what ever the time is below
        starttime = UTCDateTime(arrival_time) - lower_crop_time 
        endtime = UTCDateTime(arrival_time) + upper_crop_time

        print(f'Trim-start= {starttime}. Trim-end = {endtime}')

        ##########
        #This is for Highland data archive .... change if stored in a different format
        year = arrival_time.year
        julian_day = arrival_time.timetuple().tm_yday

        directory_path = os.path.join(base_directory, str(year), f"{julian_day:03d}")
        if not os.path.exists(directory_path):
            print(f"Directory {directory_path} does not exist.")
            return None, []
        
        #Only reading Z... add other comps if needed/ wanted for different phase stack
        pattern_z = os.path.join(directory_path, f"{year}{julian_day:03d}_*_{station}_Z2.m") 
        file_list_z = glob.glob(pattern_z)
        ##########
        
        sample_seconds = 1 / sample_rate #From seconds to Hz

        #If mseed not found for station print
        if not file_list_z:
            print(f"File {pattern_z} does not exist.")
            return None, []

        # Read and preprocess the MiniSEED file
        st = read(file_list_z[0], starttime=readstarttime, endtime=readendtime)  # Read the full MiniSEED file. hour either side of event.
        # st.remove_response(inventory=inv)
        st.detrend("linear")  # Apply detrending
        st.taper(max_percentage=0.05, type="cosine")  # Apply tapering
        st.filter("bandpass", freqmin=bandpass_freq_min, freqmax=bandpass_freq_max, zerophase=True)  # Apply bandpass filter
        st.resample(sample_seconds)  # Resample the entire stream

        # Slice the stream for different time windows
        st = st.slice(starttime=starttime, endtime=endtime)  # Main window

        print(f"{phase_type} phase completed.")
        print(st)
        
        #Write the st as text to the .aq file
        for trace in st:
            event_file_path = os.path.join(event_save_path, f"event{event_number}.aq") 
            with open(event_file_path, 'a') as f:
                f.write(f"1 {len(trace.data)} {0} {station}\n") #Trace exist/ have data = 1, number of samples, realtive shift (stacked agasint ak135, not another station ergo = 0. Can change if wanetd though to stack differently), station name
                for sample_value in trace.data:
                    f.write(f"{sample_value} ") #write the trace
                f.write("\n")

        return st, []
    except Exception as e:
        print(f"Failed to read MiniSEED file for station {station} and event {event_number}: {e}")
        return None, []

def process_event(event_row, stations_df, phase_type, sample_rate, mseed_path, lower_crop_time, upper_crop_time, bandpass_freq_min, bandpass_freq_max, event_save_path, inv):
    #From events.csv ... assuming an online download, e.g., from IRIS to csv. Remove csv header berfore running
    event_coords = (event_row[1], event_row[2], event_row[3])
    event_number = event_row[0]
    event_time = datetime.datetime.fromisoformat(event_row[4][:-1])  # Remove 'Z' at the end from the time in the event.csv file

    stations_with_data = []
    
    ## Write_event_file header
    write_event_file(event_number, 0, event_coords, event_time, sample_rate, phase_type, event_save_path)

    print("\n" + "#" * 100)
    print(f"Event: {event_time}")
    print("#" * 100 + "\n")

    # For inital alignment... can use any reference TauPy model
    m = TauPyModel(model="ak135")
    
    for station_row in stations_df.itertuples(index=False):
        station_coords = (station_row[1], station_row[2])
        distance = locations2degrees(event_coords[0], event_coords[1], station_coords[0], station_coords[1])

        ## Make sure not splitting phase type e.g., read full PKiKP not P, K .....
        if type(phase_type) != list:
            phase_type = [phase_type]

        ## Calc times for desired phase
        arrivals = m.get_travel_times(source_depth_in_km=event_coords[2], distance_in_degree=distance, phase_list=phase_type)
        distance_azi, azimuth, back_azimuth = gps2dist_azimuth(event_coords[0], event_coords[1], station_coords[0], station_coords[1])
        print(f"For {station_row[0]} and {event_number}, backazimuth = {back_azimuth}")

        ## Loop through for each arrival phase
        for arrival in arrivals:
            ## arrival = event time + expected in seconds.
            arrival_time = event_time + datetime.timedelta(seconds=arrival.time)
            print(f"Station: {station_row[0]}, Phase: {arrival.name}, Predicted Arrival Time: {arrival_time}")
            
            st, _ = read_mseed(mseed_path, phase_type, station_row[0], arrival_time, event_number, sample_rate, upper_crop_time, lower_crop_time, bandpass_freq_min, bandpass_freq_max, event_save_path, inv)
            # Only add stations that have valid waveform data
            if st and len(st) > 0:
                stations_with_data.append(station_row[0])
            
    station_count = len(stations_with_data)    
    update_station_count(event_number, station_count, event_save_path)

if __name__ == "__main__":
    main()