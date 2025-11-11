#!/bin/bash
#
# Builds and runs the benchmark in Release mode.
#
# This script is location-aware. If run from a Windows
# mount (e.g., /mnt/d/) inside WSL, it will automatically
# copy the required files to a native Linux directory
# (e.g., ~/april_release_bench) for a faster, more stable build.
#
# 100% vibe coded
#


# --- Configuration ---
set -e

# Directory on the native WSL file system
DEST_DIR="~/april_release_bench"

# Compilers
CC_TOOL="clang"
CXX_TOOL="clang++"

# --- Initialization ---
# Expand the tilde (~) in the path
DEST_DIR=$(eval echo $DEST_DIR)
CURRENT_DIR=$(pwd)
RUN_DIR=$CURRENT_DIR
NEEDS_COPY=false

# --- Dependency Check ---
echo "🔎 Checking dependencies..."
if ! command -v $CC_TOOL &> /dev/null; then
    echo "❌ Error: $CC_TOOL is not found." >&2; exit 1;
fi
if ! command -v $CXX_TOOL &> /dev/null; then
    echo "❌ Error: $CXX_TOOL is not found." >&2; exit 1;
fi

# --- Location Check ---
# Check if we're in WSL and on a /mnt/ drive
if grep -qi "microsoft" /proc/version && [[ "$CURRENT_DIR" == /mnt/* ]]; then
    echo "⚠️  WSL + /mnt/ drive detected. Staging files to native filesystem for performance."
    NEEDS_COPY=true
    RUN_DIR=$DEST_DIR
fi
echo "---------------------------------"

# --- 1. Copy Files (if needed) ---
if $NEEDS_COPY; then
    echo "🧹 Cleaning and creating destination directory: $RUN_DIR"
    rm -rf $RUN_DIR
    mkdir -p $RUN_DIR

    echo "CI: copying required files (.h, .hpp, .cpp, CMakeLists.txt)..."
    
    # 1. Copy top-level CMakeLists.txt
    cp "$CURRENT_DIR/CMakeLists.txt" "$RUN_DIR/"
    
    # 2. Use find + cpio to copy filtered files from subdirs
    echo "CI: Copying source files from include/ and benchmark/..."
    
    # This runs 'find' from the source directory to get relative paths,
    # then pipes them to 'cpio' to recreate the structure at the destination.
    (cd "$CURRENT_DIR" && \
     find include benchmark \
     \( -name "*.h" -o -name "*.hpp" -o -name "*.cpp" -o -iname "cmakelists.txt" \) \
     -print | cpio -pdm "$RUN_DIR")
fi

# --- 2. Enter Build Directory ---
# All subsequent commands happen in the RUN_DIR
cd $RUN_DIR
echo "⚙️  Operating in: $(pwd)"

# --- 3. Modify CMakeLists ---
# We must do this every time to ensure we ONLY build the benchmark.
# This makes the script safe to run in a directory without test/ or examples/.
echo "🔧 Modifying CMakeLists.txt to build only benchmarks..."
sed -i '/add_subdirectory(examples)/s/^/#/' CMakeLists.txt
sed -i '/add_subdirectory(test)/s/^/#/' CMakeLists.txt

# --- 4. Configure ---
echo "⚙️  Configuring CMake (Build Type: Release)..."
export CC=$CC_TOOL
export CXX=$CXX_TOOL

rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# --- 5. Build ---
echo "🔨 Building project (using all cores)..."
cmake --build build -j $(nproc)

# --- 6. Run Benchmark ---
echo ""
echo "🚀 --- Running Benchmark (cube_benchmarkLL) ---"
./build/benchmark/cube_benchmarkLL
echo "--- Benchmark complete ---"