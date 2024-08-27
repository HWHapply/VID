#!/bin/bash

# Define the list of optional arguments
optional_args_list=(
    "h5ad_dir" 
    "data_dir" 
    "meta_dir"
    "output_dir"
    "marker_dir"
    "feature_dir" 
    "clinical_column"
    "batch_column" 
    "sample_column"
    "test_ratio"
    "num_split"
    "metamodel"
    "threshold"
    "average"
    "random_state"
    "n_jobs"
    "verbose"
)

# Initialize an array to hold the optional arguments for Python
optional_args=()

# Function to display usage
usage() {
    echo "Usage: $0 <seuratobj_dir> [--h5ad_dir H5AD_DIR] [--data_dir DATA_DIR] [--meta_dir META_DIR] [--output_dir OUTPUT_DIR] [--marker_dir MARKER_DIR]
                  [--feature_dir FEATURE_DIR] [--clinical_column CLINICAL_COLUMN] [--batch_column BATCH_COLUMN] [--sample_column SAMPLE_COLUMN]
                  [--test_ratio TEST_RATIO] [--num_split NUM_SPLIT] [--metamodel METAMODEL] [--threshold THRESHOLD] [--average AVERAGE]
                  [--random_state RANDOM_STATE] [--n_jobs N_JOBS] [--verbose VERBOSE] [--help]"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --help)
            usage
            ;;
        -*)
            # Argument is in the --name format
            arg_name="${1#--}"  # Remove the leading '--'
            if [[ " ${optional_args_list[@]} " =~ " ${arg_name} " ]]; then
                # Argument is in the list of known optional arguments
                if [[ -n "$2" && ! "$2" =~ ^- ]]; then
                    optional_args+=("$1 $2")
                    shift
                else
                    echo "Error: Value for $1 is missing"
                    usage
                fi
            else
                echo "Error: Unknown argument $1"
                usage
            fi
            ;;
        *)
            # Positional argument
            if [[ -z "$seuratobj_dir" ]]; then
                seuratobj_dir="$1"
            else
                echo "Error: Unexpected positional argument $1"
                usage
            fi
            ;;
    esac
    shift
done

# Check if the positional argument is provided
if [[ -z "$seuratobj_dir" ]]; then
    echo "Error: Positional argument is required"
    usage
fi

# Create a folder named as the current time
timestamp=$(date +"%Y%m%d_%H%M%S")
rawdata_dir="./${timestamp}/data"
output_dir="./${timestamp}/output"
mkdir -p "$rawdata_dir"
mkdir -p "$output_dir"

# Run the R script with the positional argument and the output directory
rscript Input_Preparation_R.R "$seuratobj_dir" "$rawdata_dir"

# Collect files in the output directory
h5ad_file=$(find "$rawdata_dir" -maxdepth 1 -type f -name "data.h5ad")
dmatrix_file=$(find "$rawdata_dir" -maxdepth 1 -type f -name "dmatrix.csv")
metadata_file=$(find "$rawdata_dir" -maxdepth 1 -type f -name "metadata.csv")

# Pass the appropriate arguments to the Python script based on the output
if [[ -n "$h5ad_file" && -z "$dmatrix_file" && -z "$metadata_file" ]]; then
    optional_args+=("--h5ad_dir $h5ad_file")
elif [[ -n "$dmatrix_file" && -n "$metadata_file" ]]; then
    optional_args+=("--data_dir $dmatrix_file" "--meta_dir $metadata_file")
else
    echo "Error: Unexpected files in output directory"
    exit 1
fi

# Add --output_dir to optional_args if it's not already present
if [[ ! " ${optional_args[@]} " =~ " --output_dir " ]]; then
    optional_args+=("--output_dir $output_dir")
fi

# Prepare arguments for the Python script
python_args=$(printf "%s " "${optional_args[@]}")

# Run the Python script with the optional arguments
python run_vid.py $python_args

# Save the predicted infection status to seurat object
Rscript Save_Result.R "$seuratobj_dir" "$rawdata_dir"