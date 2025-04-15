#!/bin/bash

# Initialize variables
args=""
feature_file=""
marker_file=""
label_file=""

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
container_label_path="/wkdir/input/labels.txt"

# Check if the feature file is mounted (i.e., -v is used to mount a file)
if [[ -f "$container_feature_path" ]]; then
    feature_file="--feature_dir $container_feature_path"
fi

# Check if the feature file is mounted (i.e., -v is used to mount a file)
if [[ -f "$container_marker_path" ]]; then
    marker_file="--marker_dir $container_marker_path"
fi

# Check if the feature file is mounted (i.e., -v is used to mount a file)
if [[ -f "$container_label_path" ]]; then
    label_file="--label_dir $container_label_path"
fi

# remove matplotlib cache
rm -rf ~/.cache/matplotlib

# Run the main script with the provided arguments
bash /opt/VID/bin/run_vid "$container_input_path" --output_dir "$container_output_dir"  $feature_file $args $marker_file $label_file

echo "VID process complete!"

