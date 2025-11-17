#!/bin/bash
#
# Copies the essential source files (CMake, include, test, scripts, benchmark, examples)
# from the WSL /mnt/ drive to a native Linux directory (~/april)
# to ensure stable, high-performance builds and coverage runs.
#

# --- Configuration ---
set -e

# Source path on your Windows drive
SOURCE_REPO_PATH="/mnt/d/dev/april"

# Destination path on your native Linux filesystem
DEST_DIR="~/april"

# --- Setup ---
# Expand the tilde (~) in the destination path
DEST_DIR=$(eval echo "$DEST_DIR")

echo "Staging project for native build..."
echo "  Source: $SOURCE_REPO_PATH"
echo "  Target: $DEST_DIR"

# 1. Check if source exists
if [ ! -d "$SOURCE_REPO_PATH" ]; then
    echo "❌ Error: Source directory $SOURCE_REPO_PATH not found." >&2
    exit 1
fi

# 2. Clean and create the destination directory
echo "🧹 Cleaning destination..."
rm -rf "$DEST_DIR"
mkdir -p "$DEST_DIR"

# 3. Copy top-level CMakeLists.txt
echo "CI: copying top-level CMakeLists.txt..."
cp "$SOURCE_REPO_PATH/CMakeLists.txt" "$DEST_DIR/"

# 4. Copy filtered contents of subdirectories
# We 'cd' into the source repo to run 'find', which gives us
# the relative paths that 'cpio' needs to recreate the
# directory structure at the destination.
echo "CI: copying sources from 'include/', 'test/', 'scripts/', 'benchmark/', and 'examples/'..."
(cd "$SOURCE_REPO_PATH" && \
 { \
   find scripts -name "*.sh" ; \
   find include test examples \
     \( -name "*.h" -o -name "*.hpp" -o -name "*.cpp" -o -iname "cmakelists.txt" \) \
     -print ; \
   find benchmark \
     -type d -name env -prune -o \
     \( -name "*.h" -o -name "*.hpp" -o -name "*.cpp" -o -iname "cmakelists.txt" \) \
     -print ; \
 } | cpio -pdm "$DEST_DIR")

echo ""
echo "✅ Project files successfully staged."
echo "You can now run your scripts from the native directory:"
echo "cd $DEST_DIR"
echo "./scripts/run_coverage.sh"