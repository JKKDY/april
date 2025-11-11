#!/bin/bash
#
# Installs all necessary build, test, and coverage tools
# for the April project on an Ubuntu-based WSL.
#
# This script will:
# 1. Install base utilities (cmake, gcovr, cpio, dos2unix, wget).
# 2. Install GCC 14 + G++ 14 from the PPA.
# 3. Install Clang 21 + Clang++ 21 from apt.llvm.org.
# 4. Configure clang-21 as the primary, default C/C++ compiler.
#

# Stop the script if any command fails
set -e

echo "--- Starting April dev environment setup ---"

# 1. Update package list
sudo apt-get update

# 2. Install Base Utilities
echo "--- Installing base utilities (cmake, gcovr, cpio, dos2unix, wget)... ---"
sudo apt-get install -y \
    cmake \
    gcovr \
    cpio \
    dos2unix \
    wget \
    lsb-release \
    build-essential

# 3. Install GCC 14 (for secondary testing and gcov-14)
echo "--- Installing GCC 14 ---"
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y gcc-14 g++-14

# 4. Install Clang 21 (primary compiler)
echo "--- Installing Clang 21 ---"
# Download the official LLVM installer script
wget https://apt.llvm.org/llvm.sh
chmod +x llvm.sh
# Run the script to install version 21
sudo ./llvm.sh 21
# Clean up the script
rm llvm.sh

# 5. Configure Default Compilers (update-alternatives)
echo "--- Configuring clang-21 as primary default compiler ---"

# --- Install all compilers as "alternatives" ---
# C Compilers
sudo update-alternatives --install /usr/bin/cc cc /usr/bin/clang-21 100
sudo update-alternatives --install /usr/bin/cc cc /usr/bin/gcc-14 50

# C++ Compilers
sudo update-alternatives --install /usr/bin/cxx cxx /usr/bin/clang++-21 100
sudo update-alternatives --install /usr/bin/cxx cxx /usr/bin/g++-14 50

# --- Set the highest priority (clang-21) as the default ---
sudo update-alternatives --set cc /usr/bin/clang-21
sudo update-alternatives --set cxx /usr/bin/clang++-21

# 6. Verification
echo ""
echo "--- Setup Complete! ---"
echo "Verifying default versions (should be clang-21):"
echo -n "Default C Compiler (cc): "
cc --version | head -n 1
echo -n "Default C++ Compiler (cxx): "
cxx --version | head -n 1
echo ""
echo "Verifying all installed tools:"
gcc-14 --version | head -n 1
gcov-14 --version | head -n 1
gcovr --version
cpio --version | head -n 1
cmake --version | head -n 1
dos2unix --version