#!/bin/bash

# ------------------------------------------------------------------------------- #
# Input values ------------------------------------------------------------------ #

fmin=0.02
fmax=1.0
plottype="both" #or both or map
testdir="Test_012020"

# ------------------------------------------------------------------------------- #
# main -------------------------------------------------------------------------- #

# Directory containing the .aq event files
EVENT_DIR="/projects/prjs1435/test_waveforms/Astack/$testdir/Input_data"
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

echo "Processing frequency band: $fmin - $fmax Hz"
if [ $fmax > $fmin ] ; then

    echo "Frequency band: $fmin - $fmax Hz"

    # **** HDF5 to txt conversion ****
    echo "Converting HDF5 files to txt files..."
    ./HDF5_to_txt.py $fmin $fmax $testdir
    echo "Done"

    # Iterate over all .aq files in the directory
    for EVENT_FILE in "$EVENT_DIR"/?????????????_"$fmin"-"$fmax"Hz.aq; do
        # Extract the event name from the file name
        EVENT_NAME=$(basename "$EVENT_FILE")

        # Update the tcas.cmd file with the new event name
        oldname=`grep aq $TCAS_CMD_FILE | awk '{print $1}'`
        sed -i "s/$oldname/$EVENT_NAME/g" "$TCAS_CMD_FILE"

        # **** Execute the tsac command ****
        echo "Processing event: $EVENT_NAME"
        ./$TSAC_EXECUTABLE

        # Check if tsac executed successfully
        if [ $? -ne 0 ]; then
            echo "Error: tsac failed for event $EVENT_NAME"
            exit 1
        fi
    done

    # **** Plotting the results ****
    echo "Plotting the results..."
    ./Plot_astack.py $fmin $fmax $plottype $testdir
    echo "Done"

    echo "All events processed successfully."
fi