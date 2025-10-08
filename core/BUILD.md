
# GRF C++ Library Build System
This directory contains the C++ core implementation of the Generalized Random Forest (GRF) library with a comprehensive build system.
## Quick Start
### Linux/macOS
```bash
# Make build script executable
chmod +x build.sh
# Build and test
./build.sh check
# Or use make directly
make test
```
### Windows
```cmd
# Build and test
build.bat check
# Or use make directly (requires MSYS2/MinGW)
make test
```
## Build Targets
### Primary Targets
- `make all` - Build library and tests
- `make library` - Build static library only
- `make test` - Build and run tests
- `make clean` - Clean build artifacts
### Development Targets
- `make debug` - Build debug version
- `make release` - Build optimized release
- `make coverage` - Build with coverage and run tests
- `make benchmark` - Run performance benchmarks
### Testing Targets
- `make run-tests` - Build and run all tests
- `make run-tests-verbose` - Run tests with verbose output
- `make run-tests-specific TEST=pattern` - Run specific tests
- `make list-tests` - List available tests
### Code Quality Targets
- `make format` - Format code (requires clang-format)
- `make lint` - Lint code (requires cpplint)
- `make analyze` - Static analysis (requires cppcheck)
### Installation Targets
- `make install [PREFIX=/path]` - Install library and headers
- `make uninstall [PREFIX=/path]` - Uninstall
## Prerequisites
### Required
- **C++ Compiler**: g++, clang++, or MSVC with C++11 support
- **Build Tools**: make, ar (archiver)
- **Threading**: pthread (Linux/macOS) or Windows threading
### Optional (for additional features)
- **clang-format**: Code formatting
- **cppcheck**: Static analysis
- **cpplint**: Code linting
- **gcov/lcov**: Code coverage
- **doxygen**: Documentation generation
- **valgrind**: Memory debugging (Linux)
## Platform-Specific Setup
### Ubuntu/Debian
```bash
sudo apt-get install build-essential g++ make
# Optional tools
sudo apt-get install clang-format cppcheck valgrind lcov doxygen
```
### CentOS/RHEL/Fedora
```bash
sudo yum install gcc-c++ make
# Or for newer versions:
sudo dnf install gcc-c++ make
```
### macOS
```bash
# Install Xcode Command Line Tools
xcode-select --install
# Or use Homebrew
brew install gcc make
```
### Windows
#### Option 1: MSYS2 (Recommended)
1. Install [MSYS2](https://www.msys2.org/)
2. Open MSYS2 terminal and run:
```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-make
# Add C:\msys64\mingw64\bin to your PATH
```
#### Option 2: Visual Studio
1. Install Visual Studio 2017 or later with C++ workload
2. Use Developer Command Prompt
3. May require CMake instead of make
## Directory Structure
```
core/
├── Makefile              # Main build file
├── config.mk             # Build configuration
├── build.sh              # Build script (Linux/macOS)
├── build.bat             # Build script (Windows)
├── src/                  # Source code
│   ├── analysis/         # Analysis components
│   ├── commons/          # Common utilities
│   ├── forest/           # Forest implementations
│   ├── prediction/       # Prediction strategies
│   ├── relabeling/       # Relabeling strategies
│   ├── sampling/         # Sampling methods
│   ├── splitting/        # Splitting rules
│   └── tree/             # Tree components
├── test/                 # Test files
│   ├── setup.cpp         # Test main (Catch framework)
│   └── */                # Test categories
├── third_party/          # Third-party dependencies
│   ├── Eigen/            # Linear algebra library
│   ├── optional/         # Optional header
│   └── random/           # Random number generation
└── build/                # Build output (created during build)
    ├── obj/              # Object files
    ├── libgrf.a          # Static library
    └── grf_test          # Test executable
```
## Configuration
### Environment Variables
- `CXX`: C++ compiler (default: g++)
- `CXXFLAGS`: Additional compiler flags
- `PREFIX`: Installation prefix (default: /usr/local)
- `JOBS`: Number of parallel build jobs
### Custom Configuration
Edit `config.mk` to customize build settings:
- Compiler paths
- Platform-specific flags
- Additional include/library paths
- Feature toggles
## Testing
The library uses the [Catch2](https://github.com/catchorg/Catch2) testing framework.
### Running Tests
```bash
# All tests
make run-tests
# Verbose output
make run-tests-verbose
# Specific tests
make run-tests-specific TEST="*Regression*"
# With coverage
make coverage
# List available tests
make list-tests
```
### Test Categories
- **Analysis**: Split frequency computation
- **Commons**: Utility classes and data structures
- **Forest**: Forest training and prediction
- **Prediction**: Prediction strategy implementations
- **Relabeling**: Sample relabeling strategies
- **Sampling**: Random sampling methods
- **Splitting**: Tree splitting rules
- **Tree**: Tree construction and operations
## Debugging
### Debug Build
```bash
make debug
gdb ./build/grf_test
```
### Memory Checking (Linux)
```bash
make valgrind
```
### Sanitizers
```bash
# Address sanitizer
CXXFLAGS="-fsanitize=address" make test
# Thread sanitizer
CXXFLAGS="-fsanitize=thread" make test
```
## Performance
### Optimized Build
```bash
make release
```
### Benchmarking
```bash
make benchmark
```
### Profiling
```bash
CXXFLAGS="-pg" make test
./build/grf_test
gprof ./build/grf_test gmon.out > profile.txt
```
## Integration
### Using the Library
After building, you can link against the static library:
```cpp
#include <forest/RegressionForest.h>
// ... your code
```
```bash
g++ -std=c++11 your_code.cpp -L./build -lgrf -I./src -I./third_party
```
### Installation
```bash
# System-wide installation
sudo make install
# Custom location
make install PREFIX=/opt/grf
# In your projects
g++ your_code.cpp -lgrf -I/opt/grf/include
```
## Troubleshooting
### Common Issues
1. **Compiler not found**
   - Install g++ or clang++
   - Check PATH environment variable
2. **Make not found (Windows)**
   - Install MSYS2 or use build.bat
   - Ensure make is in PATH
3. **Directory creation errors**
   - Run the directory setup script first:
     - Linux/macOS: `chmod +x setup_dirs.sh && ./setup_dirs.sh`
     - Windows: `setup_dirs.bat`
   - Or manually create: `mkdir -p build/obj/src/prediction/collector`
4. **Threading errors**
   - Ensure pthread is available
   - On Windows, use MinGW-w64
4. **Missing headers**
   - Headers are included in third_party/
   - No external dependencies required
5. **Test failures**
   - Run with verbose output: `make run-tests-verbose`
   - Check compiler version (need C++11)
   - Try debug build: `make debug`
### Getting Help
- Check build output for specific error messages
- Ensure all prerequisites are installed
- Try a clean rebuild: `make clean && make test`
- Use the issue tracker for bugs
## License
This build system is part of the GRF project and follows the same license terms.
