#!/bin/bash

# Initialize variables
args=""
feature_file=""
vid_file=""

# Parse other arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        *) args="$args $1" ;;  # Collect other arguments
    esac
    shift
done

# Define paths inside the container
container_input_path="/wkdir/input/data.rds"
container_output_dir="/wkdir/output"
container_marker_path="/wkdir/input/markers.txt"
container_feature_path="/wkdir/input/features.txt"
container_vid_path="/wkdir/input/vid.pkl"

# Check if the feature file is mounted (i.e., -v is used to mount a file)
if [[ -f "$container_feature_path" ]]; then
    feature_file="--feature_dir $container_feature_path"
fi


# Check if the feature file is mounted (i.e., -v is used to mount a file)
if [[ -f "$container_vid_path" ]]; then
    vid_file="--vidmodel_dir $container_vid_path"
fi


# Run the main script with the provided arguments
bash /opt/VID/run_vid.sh "$container_input_path" --output_dir "$container_output_dir" --marker_dir "$container_marker_path" $feature_file $args $vid_file 

echo "VID process complete! Outputs written to $output_dir"

