#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import glob
import sys

# Define the parameters
npoints = 100  # Number of points in the trace
sampr = 0.05   # Sampling rate
tshft = 0   # Time shift
swpol = 1      # Polarity switch
maxd = 1.0e-5  # Maximum amplitude

# get list of files in the directory
# Define the directory containing the files

fmin = sys.argv[1]
fmax = sys.argv[2]
plots= sys.argv[3] #'arrivals' or 'both' or 'waveforms'
ptype='asf'

if plots != 'arrivals' and ptype == 'initial':
    directory = "/projects/prjs1435/test_waveforms/Astack/SNR_test/Input_data"
else:
    directory = "/projects/prjs1435/test_waveforms/Astack/SNR_test/Output_data"

savedir='/projects/prjs1435/test_waveforms/Astack/SNR_test/Figures/'

print("Directory:", directory, "savedir:", savedir, "plots:", plots)
# check if the director exists and otherwise create it
if not os.path.exists(savedir):
    os.makedirs(savedir)

if plots == 'arrivals' or plots == 'both':
    print("Plotting arrivals")    
    # Define the pattern to match files
    ttr_pattern = os.path.join(directory, "ts??????_??????_"+fmin+"-"+fmax+"Hz.ttr")

    # List all files matching the pattern
    ttr_files = glob.glob(ttr_pattern)
    # Print the list of files
    for file in ttr_files:
        print(file)
        filename = os.path.basename(file)
        event = filename[2:15]
        print(event)
        # Initialize lists to store data
        names = []
        values = []
        errors = []

        station_file = '/projects/prjs1435/test_waveforms/files/Input_files/deepNL_station_locations.txt'
        
        with open(file, "r") as file:
            nr_stat=file.readline()
            ev_lld=file.readline()
            evla=ev_lld.split()[0]
            evlo=ev_lld.split()[1]
            evd=ev_lld.split()[2]

            lines = file.readlines()[6:]  # Start from line 9 (index 8 in Python)
            for line in lines:
                parts = line.split()
                if len(parts) >= 4:  # Ensure the line has enough columns
                    names.append(parts[1])  # Second column: station names
                    values.append(float(parts[2]))  # Third column: values
                    errors.append(float(parts[3]))  # Fourth column: error bars

        # Convert values and errors to NumPy arrays for plotting
        values = np.array(values)
        errors = np.array(errors)
        
        # Plot the data
        plt.figure(figsize=(10, 6))
        plt.errorbar(names, values, yerr=errors, fmt='o', capsize=5, color='blue')

        # Customize the plot
        plt.xlabel("Station Names")
        plt.ylabel("Values")
        plt.title("P-wave arrival times "+event)
        plt.xticks(rotation=45, ha='right')  # Rotate X-axis labels for better readability
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.savefig(savedir+'Astack_'+event+'_arrivals.png')
        print(savedir+'Astack_'+event+'_arrivals.png')
        station_latitudes = []
        station_longitudes = []
        distev = []

        for name in names:
            # get station names, latitude and longitude from station file
            with open(station_file, "r") as file:
                lines = file.readlines()
                for line in lines:
                    parts = line.split()
                    if len(parts) >= 3:  # Ensure the line has enough columns
                        if name == parts[0]:
                            lat = float(parts[1])
                            lon = float(parts[2])
                            station_latitudes.append(lat)
                            station_longitudes.append(lon)
                            #calculate distance between event and station
                            dist = np.sqrt((float(evla)-lat)**2 + (float(evlo)-lon)**2)
                            distev.append(dist)
                            break
                else:
                    print(f"Station {name} not found in the station file")
                    continue
        
        fig = plt.figure(figsize=(7, 7))
        m = Basemap(projection='merc', llcrnrlat=52.4, urcrnrlat=54, llcrnrlon=5.2, urcrnrlon=8, resolution='i')
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        m.drawmapboundary(fill_color='lightblue')
        m.fillcontinents(color='lightgreen', lake_color='aqua')
        fig.show()
        norm = plt.Normalize(np.min(values), np.max(values))
        cmap = plt.cm.coolwarm_r

        # Convert latitudes and longitudes to map coordinates
        x, y = m(station_longitudes, station_latitudes)
        sc = m.scatter(x, y, c=values, cmap=cmap, norm=norm, s=100, edgecolor='k', zorder=5)

        cb = m.colorbar(sc, location='right', pad='5%')
        cb.set_label('dt values')

        # Add station names as labels
        for i, name in enumerate(names):
                x_name, y_name = x[i]-5000, y[i]+1000
                plt.text(x_name, y_name, name, fontsize=10, ha='right', color='blue')

        event_x, event_y = m(float(evlo), float(evla))

        # Find the closest station
        closest_station_index = np.argmin(distev)

        if x[closest_station_index] > float(evlo):
            value=-1
            txt_length = +50000
        else:
            value=1
            txt_length = -5000

        closest_station_index = np.argmin(distev)
        print(closest_station_index)
        closest_station_x = x[closest_station_index]+50000*value
        closest_station_y = y[closest_station_index]-50000

        x_ref, y_ref = m(x[closest_station_index]+0.1, y[closest_station_index]-0.1)

        dx = event_x - closest_station_x
        dy = event_y - closest_station_y
        # Normalize direction vector
        length = np.sqrt(dx**2 + dy**2)
        dx /= length
        dy /= length
        arrow_length = 50000 
        # Plot arrow in the direction of the event
        plt.annotate('', xy=(closest_station_x + dx * arrow_length, closest_station_y + dy * arrow_length), xytext=(closest_station_x, closest_station_y),
                        arrowprops=dict(arrowstyle="<-", color='darkgreen',lw=3))

        evlor=round(float(evlo),2)
        evlar=round(float(evla),2)
        evtext = f"{evlor}, {evlar}"
        print(evtext)
        #txt_length = +60000 #-10000 #
        plt.annotate(evtext, (closest_station_x + dx * txt_length, closest_station_y + dy * txt_length), fontsize=10)


        plt.title("P-wave relative arrival times "+event)
        plt.legend()
        plt.savefig(savedir+'Astack_'+event+'_map.png')
        print(savedir+'Astack_'+event+'_map.png')

# plot seismograms
if plots == 'waveforms' or plots == 'both':
    print("Plotting waveforms")
    if ptype == 'initial':
        aq_pattern = os.path.join(directory, "??????_??????_"+fmin+"-"+fmax+"Hz.aq")
    elif ptype == 'asi':
        aq_pattern = os.path.join(directory, "asi??????_??????_"+fmin+"-"+fmax+"Hz.aq")
    else:
        aq_pattern = os.path.join(directory, "asf??????_??????_"+fmin+"-"+fmax+"Hz.aq")
    print(aq_pattern)
    # List all files matching the pattern
    aq_files = glob.glob(aq_pattern)
    print("aq_files:", aq_files)
    # Print the list of files
    for aqfile in aq_files:
        print(aqfile)
        filename = os.path.basename(aqfile)
        if ptype == 'initial':
            event = filename[0:13]
        else:
            event = filename[3:16]
        print(event)

        with open(aqfile, "r") as file:
            # Read the number of stations
            ns = int(file.readline().strip())
            
            # Read the event latitude, longitude, and depth
            elat, elon, depth = map(float, file.readline().strip().split())
            
            # Read the event date (year, month, day)
            year, month, day = map(int, file.readline().strip().split())
            
            # Read the event time (hour, minute, second)
            if ptype == 'initial':
                hour, minute, second, zero = map(float, file.readline().strip().split())
            else:
                hour, minute, second = map(float, file.readline().strip().split())

            # Read the sampling rate
            if ptype == 'initial':
                sampr, phase = file.readline().strip().split()
            else:
                sampr = file.readline().strip()

            sampr = float(sampr)
            
            # Format the label (similar to the Fortran `write` statement)
            label = f"{aqfile} {elat:7.2f} {elon:7.2f} {depth:6.1f}km  {day:02}/{month:02}/{year:04} {int(hour):02}:{int(minute):02}:{second:05.2f}"
            #print("Label:", label)
            
            # loop over the stations and get the header and the data for each
            plt.figure(figsize=(10, 6))
            for i in range(ns):
                # Read the header of the first station
                line = file.readline().strip().split()
                swpol = float(line[0])  # Polarity switch
                npoints = int(line[1])  # Number of points
                tshft = float(line[2])  # Time shift
                sname = line[3]  # Station name (string)
                
                # Print the station header information
                print(f"Station Header: swpol={swpol}, npoints={npoints}, tshft={tshft}, sname={sname}")

                #data=file.readline().strip()
                #print(data)
                #get values from one line of numbers in an array from the file
                data = np.array([float(x) for x in file.readline().strip().split()])
                if sname == "zssl":
                    data=data*20000
                elif sname == "zscp":
                    data=data*10000
                else:
                    data=data*10
                
                if maxd > 1.0e-6:
                    amp = swpol * data / (1.333 * maxd) + np.arange(1, npoints + 1)
                else:
                    amp = swpol * data + np.arange(1, npoints + 1)

                # Apply time shift
                ptim = np.arange(npoints) * sampr + tshft

                # Set color for each station
                #color_map = {
                #    'NE301': 'red',     'NE302': 'orange',  'NE303': 'pink',
                #   'NE304': 'brown',   'NE305': 'purple',  'NE306': 'magenta',
                #    'NE307': 'olive',  'NE308': 'lime',    'NE309': 'green',
                #    'NE310': 'teal',    'NE311': 'cyan',    'NE312': 'indigo',
                #    'NE317': 'blue',    'NE318': 'navy'
                #    }
                color_map = {
                    'NE301': 'blue',        'NE302': 'deepskyblue', 'NE303': 'lightskyblue',
                    'NE304': 'pink',        'NE305': 'hotpink',     'NE306': 'deeppink',
                    'NE307': 'crimson',     'NE308': 'red',         'NE309': 'orangered',
                    'NE310': 'darkorange',  'NE311': 'orange',      'NE312': 'olive',
                    'NE317': 'green',       'NE318': 'lime'
                }
                # Offset each waveform by its station number (i+1) for separation                
                if sname == "zssl":
                    plt.plot(ptim, amp - (i+1)*4*1e8, label=sname, color='grey')
                elif sname == "zscp":
                    plt.plot(ptim, amp - (i+10)*5*1e8, label=sname, color='black')
                else:
                    plt.plot(ptim, amp + 10e9 - i*1e9, label=sname, color=color_map[sname])

        # Customize the plot
        plt.xlabel("Time (s)")
        #plt.ylabel("Station Number")
        plt.title("P-wave arrivals for event "+event+", frequency band: "+fmin+'-'+fmax+'Hz')
        plt.grid(True)

        # Show the plot
        plt.legend(loc='upper right', framealpha=1)
        plt.savefig(savedir+'Astack_'+event+'_seismograms_'+ptype+'.png')
        print(savedir+'Astack_'+event+'_seismograms_'+ptype+'.png')