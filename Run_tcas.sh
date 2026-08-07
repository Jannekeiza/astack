#!/bin/bash

# ------------------------------------------------------------------------------- #
# Input values ------------------------------------------------------------------ #

fmin=0.02
fmax=0.2
plottype="both" #or both or map
sample_rate=20
dataset="Dictum"  # or FDSN or Groningen
run_scripts="123"  # 1 = HDF5 to txt conversion, 2 = tsac execution, 3 = plotting, 12 = HDF5 to txt conversion and tsac execution, 13 = HDF5 to txt conversion and plotting, 23 = tsac execution and plotting, 123 = all three steps

# ------------------------------------------------------------------------------- #
# main -------------------------------------------------------------------------- #
SCRIPT_DIR="/home/laatji/bin/astack"

# Directory containing the .aq event files
ASTACK_DIR="/datasets/itc/gaia/deepnl/Waveforms/$dataset/Astack/HH_data_0.02-0.2Hz"  # Base directory for input and output data
SEIS_DIR="/datasets/itc/gaia/deepnl/Waveforms/$dataset/seismograms_20Hz"  # Directory containing the HDF5 files
EVENT_DIR="$ASTACK_DIR/Input_data"  # Directory containing the .aq event files

TCAS_CMD_FILE="$SCRIPT_DIR/tcas.cmd"         # Path to the tcas.cmd file
TSAC_EXECUTABLE="$SCRIPT_DIR/tcas"           # Path to the tsac executable

# Check if the event directory exists

for dir in "Input_data" "Output_data" "Figures" ; do
    if [ ! -d "$ASTACK_DIR/$dir" ]; then
        echo "Directory $ASTACK_DIR/$dir does not exist."
        if mkdir -p "$ASTACK_DIR/$dir"; then
            echo "Directory $ASTACK_DIR/$dir created successfully."
        else
            echo "Error: Failed to create directory $ASTACK_DIR/$dir."
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
if [ $fmax > $fmin ] ; then

    echo "Frequency band: $fmin - $fmax Hz"

    if [[ $run_scripts == *"1"* ]]; then
        # **** HDF5 to txt conversion ****
        echo "Converting HDF5 files to txt files..."
        $SCRIPT_DIR/HDF5_to_txt.py $fmin $fmax $sample_rate $ASTACK_DIR $SEIS_DIR $dataset
        echo "Done"
        success_1=True
    fi

    # Iterate over all .aq files in the directory
    #echo "event dir = $EVENT_DIR"
    #echo "astack dir = $ASTACK_DIR"

    if [[ $run_scripts == *"2"* ]]; then
        # replace line 9 of tcas.cmd with $ASTACK_DIR
        sed -i "9s|.*|$ASTACK_DIR|" "$TCAS_CMD_FILE"

        for EVENT_FILE in "$EVENT_DIR"/?????????????_"$fmin"-"$fmax"Hz.aq; do
            # Extract the event name from the file name
            EVENT_NAME=$(basename "$EVENT_FILE")

            # Update the tcas.cmd file with the new event name
            oldname=`grep aq $TCAS_CMD_FILE | awk '{print $1}'`
            sed -i "s/$oldname/$EVENT_NAME/g" "$TCAS_CMD_FILE"

            # **** Execute the tsac command ****
            echo "Processing event: $EVENT_NAME"
            $SCRIPT_DIR/$TSAC_EXECUTABLE

            # Check if tsac executed successfully
            if [ $? -ne 0 ]; then
                echo "Error: tsac failed for event $EVENT_NAME"
                exit 1
            fi
        done
        success_2=True
    fi

    if [[ $run_scripts == *"3"* ]]; then
        # **** Plotting the results ****
        echo "Plotting the results..."
        $SCRIPT_DIR/Plot_astack.py $fmin $fmax $plottype $sample_rate $ASTACK_DIR $dataset
        echo "Done"
        success_3=True
    fi

    if [[ $success_1 == True && $success_2 == True && $success_3 == True ]]; then
        echo "All steps completed successfully."
    elif [[ $success_1 == True && $success_2 == True ]]; then
        echo "HDF5 to txt conversion and tsac execution completed successfully."
    elif [[ $success_1 == True && $success_3 == True ]]; then
        echo "HDF5 to txt conversion and plotting completed successfully."
    elif [[ $success_2 == True && $success_3 == True ]]; then
        echo "tsac execution and plotting completed successfully."
    elif [[ $success_1 == True ]]; then
        echo "HDF5 to txt conversion completed successfully."
    elif [[ $success_2 == True ]]; then
        echo "tsac execution completed successfully."
    elif [[ $success_3 == True ]]; then
        echo "Plotting completed successfully."
    else
        echo "No steps were executed."
    fi

fi