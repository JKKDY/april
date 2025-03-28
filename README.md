

# Build Instructions



### Step 1: Install vcpkg
Make sure you have vcpkg on your system. If not, you can install it with:

##### On Windows (PowerShell or Command Prompt):
```sh
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```

##### On Linux/macOS:
```bash
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
```

Optional (recommended for Visual Studio users):
```sh
vcpkg integrate install
```

### Step 2: Clone APRIL

Download repository with 
```bash
git clone https://github.com/JKKDY/april.git
cd april
```

### Step 3: Build using vcpkg and CMake
##### On Windows:
```powershell
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```

##### On Linux/macOS:
```bash
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=~/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build
```
