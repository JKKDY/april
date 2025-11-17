#!/bin/bash
#
# Generates a code coverage report for the 'April' library.
# This script must be run from the project's root directory.
#
# Dependencies:
#   - cmake
#   - gcov-14 (or your matching gcov tool)
#   - gcovr (Run: sudo apt install gcovr)
#
# completely vibe coded
#

# --- Configuration ---
# Stop the script if any command fails
set -e

# The directory to build in
BUILD_DIR="build-coverage"

# The specific gcov executable to match your compiler
# (g++-14 uses gcov-14)

C_COMPILER="gcc-14"    
CXX_COMPILER="g++-14"  
GCOV_TOOL="gcov-14"    

# --- Dependency Check ---
echo "🔎 Checking dependencies..."
if ! command -v cmake &> /dev/null; then
    echo "❌ Error: cmake is not installed." >&2
    exit 1
fi
if ! command -v $GCOV_TOOL &> /dev/null; then
    echo "❌ Error: $GCOV_TOOL is not found." >&2
    echo "   (Did you install gcc-14?)" >&2
    exit 1
fi
if ! command -v gcovr &> /dev/null; then
    echo "❌ Error: gcovr is not installed." >&2
    echo "   (Run: sudo apt install gcovr)" >&2
    exit 1
fi
echo "✅ Dependencies found."


# --- 1. Configure ---
echo "⚙️  Configuring CMake (Build Type: Debug, Coverage: ON)..."
rm -rf $BUILD_DIR
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=$C_COMPILER \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
    -DAPRIL_ENABLE_COVERAGE=ON

# --- 2. Build ---
echo "🔨 Building project..."
# Use $(nproc) to build using all available CPU cores
cmake --build $BUILD_DIR -j $(nproc)

# --- 3. Run Tests ---
echo "🧪 Running tests (to generate coverage data)..."
# We run ctest from inside the build directory
(cd $BUILD_DIR && ctest --output-on-failure)

# --- 4. Generate Report ---
echo "📊 Generating coverage report..."
gcovr --gcov-executable $GCOV_TOOL \
      --root . \
      --filter 'include/april/.*' \
      --html-details "$BUILD_DIR/coverage_report.html" \
      --print-summary \
      "$BUILD_DIR"

# --- 5. View Report ---
REPORT_PATH="file://$(pwd)/$BUILD_DIR/coverage_report.html"
echo "✅ Report generation complete."
echo "   View the report: $REPORT_PATH"

# Optional: automatically open the report
if command -v xdg-open &> /dev/null; then
    xdg-open $REPORT_PATH
fi