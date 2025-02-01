#!/bin/bash

# Exit immediately if a command fails
set -e

echo "Starting VID setup..."

# Detect the operating system
OS_TYPE=$(uname)
case "$OS_TYPE" in
    "Linux")
        ENV_FILE="./env/vid_env_Linux.yml"
        SHELL_CONFIG="$HOME/.bashrc"
        ;;
    "Darwin")
        ENV_FILE="./env/vid_env_macOS.yml"
        SHELL_CONFIG="$HOME/.zshrc"
        ;;
    "CYGWIN"* | "MINGW"* | "MSYS"*)
        ENV_FILE="./env/vid_env_Windows.yml"
        SHELL_CONFIG="$HOME/.bashrc" # Adjust if using a different shell in Windows
        ;;
    *)
        echo "Unsupported OS: $OS_TYPE"
        exit 1
        ;;
esac

echo "Detected OS: $OS_TYPE"
echo "Using environment file: $ENV_FILE"

# Create Conda environment without activating it
if [ -f "$ENV_FILE" ]; then
    conda env create -f "$ENV_FILE" -n vid_env
else
    echo "Error: $ENV_FILE not found!"
    exit 1
fi

# Add VID/bin to PATH if not already added
VID_BIN="$(pwd)/bin"
if ! grep -q "$VID_BIN" "$SHELL_CONFIG"; then
    echo "Updating $SHELL_CONFIG to add VID/bin to PATH..."
    echo "export PATH=\"$VID_BIN:\$PATH\"" >> "$SHELL_CONFIG"
fi

# Reload shell configuration
source "$SHELL_CONFIG"

echo "Setup complete!"
echo "To use the environment, run: conda activate vid_env"

