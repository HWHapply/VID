# Step 1: Create a new environment with the provided .yml file
echo "Creating new environment from vid_env.yml..."
conda env create -f vid_env.yml -n vid_env

# Step 2: Activate the environment
echo "Activating the environment..."
conda activate vid_env

echo "Environment setup complete!"

# Step 3: Create directory and download demo data
echo "Creating demo directory and downloading data..."
mkdir -p ./demo/data
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/bdkv2napos1md1uca2wg8/demo.rds?rlkey=bhe5deyz2o6kenj2s2fypxkzv&st=8armlfka&dl=1' -O ./demo/data/demo.rds

echo "Demo dataset installed!
