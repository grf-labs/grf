
#!/bin/bash
# GRF Build Script
# Convenient wrapper around the Makefile for common tasks
set -e
# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color
# Print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}
print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}
print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}
print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}
# Check if we're in the right directory
if [ ! -f "Makefile" ] || [ ! -d "src" ] || [ ! -d "test" ]; then
    print_error "This script must be run from the core directory containing Makefile, src/, and test/"
    exit 1
fi
# Function to show help
show_help() {
    echo "GRF Build Script"
    echo "================"
    echo ""
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  build          Build the library and tests"
    echo "  test           Build and run all tests"
    echo "  test-fast      Build and run tests with minimal output"
    echo "  test-verbose   Build and run tests with verbose output"
    echo "  test-specific  Run specific tests (requires TEST environment variable)"
    echo "  library        Build only the static library"
    echo "  debug          Build debug version"
    echo "  release        Build optimized release version"
    echo "  clean          Clean all build artifacts"
    echo "  coverage       Build with coverage and generate report"
    echo "  benchmark      Run performance benchmarks"
    echo "  format         Format all source code"
    echo "  lint           Run code linting"
    echo "  analyze        Run static analysis"
    echo "  install        Install library and headers"
    echo "  uninstall      Uninstall library and headers"
    echo "  check          Run comprehensive checks (build, test, lint)"
    echo "  help           Show this help message"
    echo ""
    echo "Environment variables:"
    echo "  CXX            C++ compiler to use (default: g++)"
    echo "  PREFIX         Installation prefix (default: /usr/local)"
    echo "  TEST           Test pattern for test-specific command"
    echo "  JOBS           Number of parallel jobs (default: auto-detect)"
    echo ""
    echo "Examples:"
    echo "  $0 build                    # Build library and tests"
    echo "  $0 test                     # Run all tests"
    echo "  TEST='*Regression*' $0 test-specific  # Run regression tests only"
    echo "  PREFIX=/opt/grf $0 install # Install to /opt/grf"
    echo "  JOBS=4 $0 build           # Build with 4 parallel jobs"
}
# Detect number of CPU cores for parallel builds
detect_jobs() {
    if [ -n "$JOBS" ]; then
        echo "$JOBS"
    elif command -v nproc >/dev/null 2>&1; then
        nproc
    elif [ -f /proc/cpuinfo ]; then
        grep -c ^processor /proc/cpuinfo
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu 2>/dev/null || echo 1
    else
        echo 1
    fi
}
# Check dependencies
check_dependencies() {
    local missing=0
    
    if ! command -v ${CXX:-g++} >/dev/null 2>&1; then
        print_error "C++ compiler '${CXX:-g++}' not found"
        missing=1
    fi
    
    if ! command -v make >/dev/null 2>&1; then
        print_error "make not found"
        missing=1
    fi
    
    if ! command -v ar >/dev/null 2>&1; then
        print_error "ar (archiver) not found"
        missing=1
    fi
    
    if [ $missing -eq 1 ]; then
        print_error "Missing required dependencies. Please install them first."
        exit 1
    fi
}
# Run make with appropriate number of jobs
run_make() {
    local jobs=$(detect_jobs)
    local target="$1"
    shift
    
    print_info "Building with $jobs parallel jobs..."
    make -j"$jobs" "$target" "$@"
}
# Main command processing
case "${1:-help}" in
    build)
        print_info "Building GRF library and tests..."
        check_dependencies
        run_make all
        print_success "Build completed successfully!"
        ;;
        
    test)
        print_info "Building and running tests..."
        check_dependencies
        run_make run-tests
        print_success "All tests passed!"
        ;;
        
    test-fast)
        print_info "Running tests with minimal output..."
        check_dependencies
        print_info "Building tests..."
        run_make test
        if [ $? -ne 0 ]; then
            print_error "Build failed"
            exit 1
        fi
        print_info "Running test executable..."
        ./build/grf_test --reporter compact
        if [ $? -ne 0 ]; then
            print_error "Tests failed"
            exit 1
        fi
        print_success "Tests completed!"
        ;;
        
    test-verbose)
        print_info "Running tests with verbose output..."
        check_dependencies
        run_make run-tests-verbose
        print_success "Verbose tests completed!"
        ;;
        
    test-specific)
        if [ -z "$TEST" ]; then
            print_error "TEST environment variable must be set for specific tests"
            echo "Example: TEST='*Regression*' $0 test-specific"
            exit 1
        fi
        print_info "Running specific tests matching: $TEST"
        check_dependencies
        run_make run-tests-specific TEST="$TEST"
        print_success "Specific tests completed!"
        ;;
        
    library)
        print_info "Building static library only..."
        check_dependencies
        run_make library
        print_success "Library built successfully!"
        ;;
        
    debug)
        print_info "Building debug version..."
        check_dependencies
        run_make debug
        print_success "Debug build completed!"
        ;;
        
    release)
        print_info "Building optimized release version..."
        check_dependencies
        run_make release
        print_success "Release build completed!"
        ;;
        
    clean)
        print_info "Cleaning build artifacts..."
        make clean
        print_success "Clean completed!"
        ;;
        
    coverage)
        print_info "Building with coverage and generating report..."
        check_dependencies
        if ! command -v gcov >/dev/null 2>&1; then
            print_warning "gcov not found, coverage may not work properly"
        fi
        run_make coverage
        print_success "Coverage report generated!"
        ;;
        
    benchmark)
        print_info "Running performance benchmarks..."
        check_dependencies
        run_make benchmark
        print_success "Benchmarks completed!"
        ;;
        
    format)
        print_info "Formatting source code..."
        if ! command -v clang-format >/dev/null 2>&1; then
            print_error "clang-format not found. Please install clang-format first."
            exit 1
        fi
        make format
        print_success "Code formatting completed!"
        ;;
        
    lint)
        print_info "Running code linting..."
        if command -v cpplint >/dev/null 2>&1; then
            make lint
        elif command -v cppcheck >/dev/null 2>&1; then
            make analyze
        else
            print_warning "No linting tools found (cpplint or cppcheck)"
            print_info "Installing cppcheck is recommended for static analysis"
        fi
        print_success "Linting completed!"
        ;;
        
    analyze)
        print_info "Running static analysis..."
        if ! command -v cppcheck >/dev/null 2>&1; then
            print_error "cppcheck not found. Please install cppcheck first."
            exit 1
        fi
        make analyze
        print_success "Static analysis completed!"
        ;;
        
    install)
        print_info "Installing library and headers..."
        if [ $(id -u) -ne 0 ] && [ -z "$PREFIX" ]; then
            print_warning "Installing to system directories may require sudo"
            print_info "Use PREFIX=/path/to/install to install to a custom location"
        fi
        check_dependencies
        run_make library
        make install PREFIX="${PREFIX:-/usr/local}"
        print_success "Installation completed!"
        ;;
        
    uninstall)
        print_info "Uninstalling library and headers..."
        make uninstall PREFIX="${PREFIX:-/usr/local}"
        print_success "Uninstallation completed!"
        ;;
        
    check)
        print_info "Running comprehensive checks..."
        check_dependencies
        
        print_info "Step 1/4: Building..."
        run_make all
        
        print_info "Step 2/4: Running tests..."
        run_make run-tests
        
        print_info "Step 3/4: Checking code format..."
        if command -v clang-format >/dev/null 2>&1; then
            make format
        else
            print_warning "clang-format not found, skipping format check"
        fi
        
        print_info "Step 4/4: Running static analysis..."
        if command -v cppcheck >/dev/null 2>&1; then
            make analyze
        else
            print_warning "cppcheck not found, skipping static analysis"
        fi
        
        print_success "All checks completed successfully!"
        ;;
        
    help|--help|-h)
        show_help
        ;;
        
    *)
        print_error "Unknown command: $1"
        echo ""
        show_help
        exit 1
        ;;
esac
