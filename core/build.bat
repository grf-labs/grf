
@echo off
REM GRF Build Script for Windows
REM Convenient wrapper around the Makefile for common tasks
setlocal enabledelayedexpansion
REM Check if we're in the right directory
if not exist "Makefile" (
    echo [ERROR] Makefile not found. This script must be run from the core directory.
    exit /b 1
)
if not exist "src" (
    echo [ERROR] src directory not found. This script must be run from the core directory.
    exit /b 1
)
if not exist "test" (
    echo [ERROR] test directory not found. This script must be run from the core directory.
    exit /b 1
)
REM Function to show help
if "%1"=="" goto :show_help
if "%1"=="help" goto :show_help
if "%1"=="--help" goto :show_help
if "%1"=="-h" goto :show_help
REM Check dependencies
where g++ >nul 2>&1
if errorlevel 1 (
    echo [ERROR] g++ compiler not found. Please install MinGW-w64 or Visual Studio Build Tools.
    exit /b 1
)
where make >nul 2>&1
if errorlevel 1 (
    echo [ERROR] make not found. Please install MSYS2 or MinGW-w64.
    exit /b 1
)
REM Set default number of jobs
if "%JOBS%"=="" (
    set JOBS=%NUMBER_OF_PROCESSORS%
)
REM Process commands
if "%1"=="build" goto :build
if "%1"=="test" goto :test
if "%1"=="test-fast" goto :test_fast
if "%1"=="test-verbose" goto :test_verbose
if "%1"=="test-specific" goto :test_specific
if "%1"=="library" goto :library
if "%1"=="debug" goto :debug
if "%1"=="release" goto :release
if "%1"=="clean" goto :clean
if "%1"=="coverage" goto :coverage
if "%1"=="benchmark" goto :benchmark
if "%1"=="format" goto :format
if "%1"=="lint" goto :lint
if "%1"=="analyze" goto :analyze
if "%1"=="install" goto :install
if "%1"=="uninstall" goto :uninstall
if "%1"=="check" goto :check
echo [ERROR] Unknown command: %1
echo.
goto :show_help
:build
echo [INFO] Building GRF library and tests...
make -j%JOBS% all
if errorlevel 1 (
    echo [ERROR] Build failed
    exit /b 1
)
echo [SUCCESS] Build completed successfully!
goto :end
:test
echo [INFO] Building and running tests...
make -j%JOBS% run-tests
if errorlevel 1 (
    echo [ERROR] Tests failed
    exit /b 1
)
echo [SUCCESS] All tests passed!
goto :end
:test_fast
echo [INFO] Running tests with minimal output...
make -j%JOBS% test
if errorlevel 1 (
    echo [ERROR] Build failed
    exit /b 1
)
build\grf_test.exe --reporter compact
if errorlevel 1 (
    echo [ERROR] Tests failed
    exit /b 1
)
echo [SUCCESS] Tests completed!
goto :end
:test_verbose
echo [INFO] Running tests with verbose output...
make -j%JOBS% run-tests-verbose
if errorlevel 1 (
    echo [ERROR] Tests failed
    exit /b 1
)
echo [SUCCESS] Verbose tests completed!
goto :end
:test_specific
if "%TEST%"=="" (
    echo [ERROR] TEST environment variable must be set for specific tests
    echo Example: set TEST=*Regression* ^& %0 test-specific
    exit /b 1
)
echo [INFO] Running specific tests matching: %TEST%
make -j%JOBS% run-tests-specific TEST="%TEST%"
if errorlevel 1 (
    echo [ERROR] Tests failed
    exit /b 1
)
echo [SUCCESS] Specific tests completed!
goto :end
:library
echo [INFO] Building static library only...
make -j%JOBS% library
if errorlevel 1 (
    echo [ERROR] Library build failed
    exit /b 1
)
echo [SUCCESS] Library built successfully!
goto :end
:debug
echo [INFO] Building debug version...
make -j%JOBS% debug
if errorlevel 1 (
    echo [ERROR] Debug build failed
    exit /b 1
)
echo [SUCCESS] Debug build completed!
goto :end
:release
echo [INFO] Building optimized release version...
make -j%JOBS% release
if errorlevel 1 (
    echo [ERROR] Release build failed
    exit /b 1
)
echo [SUCCESS] Release build completed!
goto :end
:clean
echo [INFO] Cleaning build artifacts...
make clean
echo [SUCCESS] Clean completed!
goto :end
:coverage
echo [INFO] Building with coverage and generating report...
where gcov >nul 2>&1
if errorlevel 1 (
    echo [WARNING] gcov not found, coverage may not work properly
)
make -j%JOBS% coverage
if errorlevel 1 (
    echo [ERROR] Coverage build failed
    exit /b 1
)
echo [SUCCESS] Coverage report generated!
goto :end
:benchmark
echo [INFO] Running performance benchmarks...
make -j%JOBS% benchmark
if errorlevel 1 (
    echo [ERROR] Benchmark failed
    exit /b 1
)
echo [SUCCESS] Benchmarks completed!
goto :end
:format
echo [INFO] Formatting source code...
where clang-format >nul 2>&1
if errorlevel 1 (
    echo [ERROR] clang-format not found. Please install clang-format first.
    exit /b 1
)
make format
echo [SUCCESS] Code formatting completed!
goto :end
:lint
echo [INFO] Running code linting...
where cpplint >nul 2>&1
if not errorlevel 1 (
    make lint
) else (
    where cppcheck >nul 2>&1
    if not errorlevel 1 (
        make analyze
    ) else (
        echo [WARNING] No linting tools found (cpplint or cppcheck)
        echo [INFO] Installing cppcheck is recommended for static analysis
    )
)
echo [SUCCESS] Linting completed!
goto :end
:analyze
echo [INFO] Running static analysis...
where cppcheck >nul 2>&1
if errorlevel 1 (
    echo [ERROR] cppcheck not found. Please install cppcheck first.
    exit /b 1
)
make analyze
echo [SUCCESS] Static analysis completed!
goto :end
:install
echo [INFO] Installing library and headers...
if "%PREFIX%"=="" set PREFIX=C:\Program Files\grf
make -j%JOBS% library
if errorlevel 1 (
    echo [ERROR] Library build failed
    exit /b 1
)
make install PREFIX="%PREFIX%"
if errorlevel 1 (
    echo [ERROR] Installation failed
    exit /b 1
)
echo [SUCCESS] Installation completed to %PREFIX%!
goto :end
:uninstall
echo [INFO] Uninstalling library and headers...
if "%PREFIX%"=="" set PREFIX=C:\Program Files\grf
make uninstall PREFIX="%PREFIX%"
echo [SUCCESS] Uninstallation completed!
goto :end
:check
echo [INFO] Running comprehensive checks...
echo [INFO] Step 1/4: Building...
make -j%JOBS% all
if errorlevel 1 (
    echo [ERROR] Build failed
    exit /b 1
)
echo [INFO] Step 2/4: Running tests...
make -j%JOBS% run-tests
if errorlevel 1 (
    echo [ERROR] Tests failed
    exit /b 1
)
echo [INFO] Step 3/4: Checking code format...
where clang-format >nul 2>&1
if not errorlevel 1 (
    make format
) else (
    echo [WARNING] clang-format not found, skipping format check
)
echo [INFO] Step 4/4: Running static analysis...
where cppcheck >nul 2>&1
if not errorlevel 1 (
    make analyze
) else (
    echo [WARNING] cppcheck not found, skipping static analysis
)
echo [SUCCESS] All checks completed successfully!
goto :end
:show_help
echo GRF Build Script for Windows
echo ============================
echo.
echo Usage: %0 [command] [options]
echo.
echo Commands:
echo   build          Build the library and tests
echo   test           Build and run all tests
echo   test-fast      Build and run tests with minimal output
echo   test-verbose   Build and run tests with verbose output
echo   test-specific  Run specific tests (requires TEST environment variable)
echo   library        Build only the static library
echo   debug          Build debug version
echo   release        Build optimized release version
echo   clean          Clean all build artifacts
echo   coverage       Build with coverage and generate report
echo   benchmark      Run performance benchmarks
echo   format         Format all source code
echo   lint           Run code linting
echo   analyze        Run static analysis
echo   install        Install library and headers
echo   uninstall      Uninstall library and headers
echo   check          Run comprehensive checks (build, test, lint)
echo   help           Show this help message
echo.
echo Environment variables:
echo   CXX            C++ compiler to use (default: g++)
echo   PREFIX         Installation prefix (default: C:\Program Files\grf)
echo   TEST           Test pattern for test-specific command
echo   JOBS           Number of parallel jobs (default: %%NUMBER_OF_PROCESSORS%%)
echo.
echo Examples:
echo   %0 build                           # Build library and tests
echo   %0 test                            # Run all tests
echo   set TEST=*Regression* ^& %0 test-specific  # Run regression tests only
echo   set PREFIX=C:\grf ^& %0 install    # Install to C:\grf
echo   set JOBS=4 ^& %0 build             # Build with 4 parallel jobs
echo.
echo Prerequisites:
echo   - MinGW-w64 or Visual Studio Build Tools
echo   - MSYS2 (recommended) or standalone make
echo   - Optional: cppcheck, clang-format for additional features
:end
endlocal
