#!/bin/bash

# ------------------------------------------------------------------------------- #
# Input values ------------------------------------------------------------------ #

fmin=0.02
fmax=1.0
plottype="both" #or both or map
sample_rate=20
overwrite=False
sources="FDSN,Dictum"  # comma-separated subset of FDSN,Dictum,Groningen to search

# ------------------------------------------------------------------------------- #
# main -------------------------------------------------------------------------- #

# Directory of this script, so it can be run from anywhere
#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="/home/laatji/bin/astack"

# Directory containing the .aq event files
BASE_DIR="/datasets/itc/gaia/deepnl/Waveforms/Astack_combined/Astack_FDSN_Dictum"  # Base directory for input and output data
WAVEFORMS_ROOT="/datasets/itc/gaia/deepnl/Waveforms"  # Root directory containing the FDSN/Dictum/Groningen seismogram trees
EVENT_DIR="$BASE_DIR/Input_data"  # Directory containing the .aq event files

TCAS_CMD_FILE="tcas.cmd"         # Path to the tcas.cmd file
TSAC_EXECUTABLE="$SCRIPT_DIR/tcas"           # Path to the tsac executable

# Check if the event directory exists

for dir in "Input_data" "Output_data" "Figures" ; do
    if [ ! -d "$BASE_DIR/$dir" ]; then
        echo "Directory $BASE_DIR/$dir does not exist."
        if mkdir -p "$BASE_DIR/$dir"; then
            echo "Directory $BASE_DIR/$dir created successfully."
        else
            echo "Error: Failed to create directory $BASE_DIR/$dir."
            exit 1
        fi
    fi
done

# Check if the tcas.cmd file exists
if [ ! -f "$TCAS_CMD_FILE" ]; then
    echo "Error: File $TCAS_CMD_FILE does not exist."
    exit 1
fi

echo "Processing frequency band: $fmin - $fmax Hz"
if (( $(echo "$fmax > $fmin" | bc -l) )); then

    echo "Frequency band: $fmin - $fmax Hz"

    # **** HDF5 to txt conversion ****
    echo "Converting HDF5 files to txt files..."
    #"$SCRIPT_DIR/HDF5_to_txt_combined.py" $fmin $fmax $sample_rate $BASE_DIR $WAVEFORMS_ROOT $overwrite $sources
    echo "Done"

    # Iterate over all .aq files in the directory
    echo "event dir = $EVENT_DIR"
    echo "base dir = $BASE_DIR"

    '''
    for EVENT_FILE in "$EVENT_DIR"/?????????????_"$fmin"-"$fmax"Hz.aq; do
        echo "Processing event file: $EVENT_FILE"
        # Extract the event name from the file name
        EVENT_NAME=$(basename "$EVENT_FILE")
        echo "Event name: $EVENT_NAME"

        # Update the tcas.cmd file with the new event name
        oldname=`grep aq $TCAS_CMD_FILE | awk '{print $1}'`
        sed -i "s/$oldname/$EVENT_NAME/g" "$TCAS_CMD_FILE"
        echo $TCAS_CMD_FILE $TSAC_EXECUTABLE
        echo "Updated tcas.cmd from $oldname to $EVENT_NAME"

        # **** Execute the tsac command ****
        echo "Processing event: $EVENT_NAME"
        "$TSAC_EXECUTABLE"

        # Check if tsac executed successfully
        if [ $? -ne 0 ]; then
            echo "Error: tsac failed for event $EVENT_NAME"
            exit 1
        fi
    done
    '''
    
    # **** Plotting the results ****
    echo "Plotting the results..."
    "$SCRIPT_DIR/Plot_astack.py" $fmin $fmax $plottype $sample_rate $BASE_DIR
    echo "Done"

    echo "All events processed successfully."
fi