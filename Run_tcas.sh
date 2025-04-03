#!/bin/bash

# Directory containing the .aq event files
EVENT_DIR="/sns/seismoml/test_waveforms/Astack/Input_data"
TCAS_CMD_FILE="tcas.cmd"         # Path to the tcas.cmd file
TSAC_EXECUTABLE="tcas"           # Path to the tsac executable

# Check if the event directory exists
if [ ! -d "$EVENT_DIR" ]; then
    echo "Error: Directory $EVENT_DIR does not exist."
    exit 1
fi

# Check if the tcas.cmd file exists
if [ ! -f "$TCAS_CMD_FILE" ]; then
    echo "Error: File $TCAS_CMD_FILE does not exist."
    exit 1
fi

fmins=(0.02 0.05 0.10 0.20) #, 0.1, 0.2, 0.5]
fmaxs=(1.0 2.0) #, 0.2, 0.5, 1.0, 2.0]
plottype="both" #or both or map

for fmin in ${fmins[@]}; do
    echo $fmin
    for fmax in ${fmaxs[@]}; do
        echo "Processing frequency band: $fmin - $fmax Hz"
        if [ $fmax > $fmin ] ; then

            echo "Frequency band: $fmin - $fmax Hz"

            # if aq file doesn exist then
            echo "Converting HDF5 files to txt files..."
            ./HDF5_to_txt.py $fmin $fmax
            echo "Done"

            # Iterate over all .aq files in the directory
            for EVENT_FILE in "$EVENT_DIR"/?????????????_"$fmin"-"$fmax"Hz.aq; do
                # Extract the event name from the file name
                EVENT_NAME=$(basename "$EVENT_FILE")

                # Update the tcas.cmd file with the new event name
                oldname=`grep aq $TCAS_CMD_FILE | awk '{print $1}'`
                sed -i "s/$oldname/$EVENT_NAME/g" "$TCAS_CMD_FILE"
                #olddir=`grep astack $TCAS_CMD_FILE`
                #sed -i "s/$olddir/$EVENT_DIR/g" "$TCAS_CMD_FILE"

                # Execute the tsac command
                echo "Processing event: $EVENT_NAME"
                ./$TSAC_EXECUTABLE

                # Check if tsac executed successfully
                if [ $? -ne 0 ]; then
                    echo "Error: tsac failed for event $EVENT_NAME"
                    exit 1
                fi
            done

            echo "Plotting the results..."
            ./Plot_astack.py $fmin $fmax $plottype
            echo "Done"

            echo "All events processed successfully."
        fi
    done
done
