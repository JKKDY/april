#!/bin/bash
#
# Generates an LLVM/Clang-based coverage report for the 'April' library.
# Run this from the project root.
#
# Dependencies:
#   - clang / clang++
#   - llvm-cov
#   - llvm-profdata
#   - cmake
#

set -e  # Exit on error

BUILD_DIR="build-coverage"

C_COMPILER="clang-21"
CXX_COMPILER="clang++-21"
PROFDATA_TOOL="llvm-profdata-21"
COV_TOOL="llvm-cov-21"

# --- Dependency Check ---
echo "🔎 Checking dependencies..."
for cmd in cmake $C_COMPILER $CXX_COMPILER llvm-cov llvm-profdata; do
    if ! command -v $cmd &> /dev/null; then
        echo "❌ Error: $cmd not found." >&2
        exit 1
    fi
done
echo "✅ Dependencies found."

# --- 1. Configure ---
echo "⚙️ Configuring CMake (Debug + Coverage)..."
rm -rf "$BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=$C_COMPILER \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
    -DAPRIL_ENABLE_COVERAGE=ON

# # --- 2. Build ---
echo "🔨 Building..."
cmake --build "$BUILD_DIR" -j "$(nproc)"

# --- 3. Run Tests ---
echo "🧪 Running tests..."

# GET ABSOLUTE PATH
# We use $(pwd) to ensure the path is absolute, regardless of where CTest runs from.
PROJECT_ROOT=$(pwd)
export LLVM_PROFILE_FILE="$PROJECT_ROOT/$BUILD_DIR/test_april-%p.profraw"

# Run CTest
ctest --test-dir "$BUILD_DIR" --output-on-failure

# --- 4. Merge Profile Data ---
echo "🔄 Merging profile data..."

# Check if files exist before merging to avoid confusing errors
count=$(ls "$BUILD_DIR"/*.profraw 2>/dev/null | wc -l)
if [ "$count" -eq 0 ]; then
    echo "❌ Error: No .profraw files generated. Did the tests actually run?"
    exit 1
fi

$PROFDATA_TOOL merge -sparse "$BUILD_DIR"/*.profraw -o "$BUILD_DIR/coverage.profdata"

# --- 5. Generate Reports ---
TEST_EXECUTABLE="$BUILD_DIR/test/test_april"

if [ ! -f "$TEST_EXECUTABLE" ]; then
    echo "❌ Error: Test executable not found at $TEST_EXECUTABLE"
    exit 1
fi

echo "📊 Generating coverage report..."

IGNORE_REGEX="test/|external/|build-coverage/"

# 5a. Terminal Summary
echo "--- Line Coverage Summary ---"
$COV_TOOL report \
    "$TEST_EXECUTABLE" \
    -instr-profile="$BUILD_DIR/coverage.profdata" \
    -ignore-filename-regex="$IGNORE_REGEX" \
    -use-color

# 5b. HTML Report
REPORT_DIR="$BUILD_DIR/coverage-report"
$COV_TOOL show \
    "$TEST_EXECUTABLE" \
    -instr-profile="$BUILD_DIR/coverage.profdata" \
    -ignore-filename-regex="$IGNORE_REGEX" \
    -format=html \
    -output-dir="$REPORT_DIR"

echo "✅ HTML Report generated at: $REPORT_DIR/index.html"