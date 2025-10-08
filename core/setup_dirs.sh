
#!/bin/bash
# GRF Directory Setup Script
# This script ensures all necessary build directories exist
set -e
echo "GRF Build Directory Setup"
echo "========================="
# Check current directory
if [ ! -f "Makefile" ] || [ ! -d "src" ] || [ ! -d "test" ]; then
    echo "ERROR: This script must be run from the core directory containing Makefile, src/, and test/"
    exit 1
fi
echo "Creating build directory structure..."
# Create main build directories
mkdir -p build/obj
echo "Scanning source tree for subdirectories..."
# Find all subdirectories in src and create corresponding obj directories
find src -type d | while read dir; do
    objdir="build/obj/$dir"
    mkdir -p "$objdir"
    echo "Created: $objdir"
done
echo "Scanning test tree for subdirectories..."
# Find all subdirectories in test and create corresponding obj directories
find test -type d | while read dir; do
    objdir="build/obj/$dir"
    mkdir -p "$objdir"
    echo "Created: $objdir"
done
echo ""
echo "Directory structure created successfully!"
echo ""
echo "Directory tree:"
find build -type d | sort | sed 's/^/  /'
echo ""
echo "Now you can run 'make test' or 'make library' safely."
