
@echo off
REM GRF Directory Setup Script for Windows
REM This script ensures all necessary build directories exist
setlocal enabledelayedexpansion
echo GRF Build Directory Setup
echo =========================
REM Check current directory
if not exist "Makefile" (
    echo ERROR: Makefile not found. This script must be run from the core directory.
    exit /b 1
)
if not exist "src" (
    echo ERROR: src directory not found. This script must be run from the core directory.
    exit /b 1
)
if not exist "test" (
    echo ERROR: test directory not found. This script must be run from the core directory.
    exit /b 1
)
echo Creating build directory structure...
REM Create main build directories
if not exist "build" mkdir build
if not exist "build\obj" mkdir build\obj
echo Scanning source tree for subdirectories...
REM Create obj directories for src subdirectories
for /d /r "src" %%i in (*) do (
    set "srcpath=%%i"
    set "relpath=!srcpath:%CD%\src\=!"
    set "objpath=build\obj\src\!relpath!"
    if not exist "!objpath!" (
        mkdir "!objpath!" 2>nul
        echo Created: !objpath!
    )
)
REM Also create for the root src directory
if not exist "build\obj\src" mkdir "build\obj\src"
echo Created: build\obj\src
echo Scanning test tree for subdirectories...
REM Create obj directories for test subdirectories
for /d /r "test" %%i in (*) do (
    set "testpath=%%i"
    set "relpath=!testpath:%CD%\test\=!"
    set "objpath=build\obj\test\!relpath!"
    if not exist "!objpath!" (
        mkdir "!objpath!" 2>nul
        echo Created: !objpath!
    )
)
REM Also create for the root test directory
if not exist "build\obj\test" mkdir "build\obj\test"
echo Created: build\obj\test
echo.
echo Directory structure created successfully!
echo.
echo Directory tree:
dir /s /b /a:d build | sort
echo.
echo Now you can run 'make test' or 'make library' safely.
endlocal
