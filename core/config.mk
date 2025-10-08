
# GRF Build Configuration
# This file contains build configuration settings for the GRF C++ library
# Compiler settings
CC=gcc
CXX=g++
# Optional: Custom compiler paths (uncomment to override)
# CXX=/usr/bin/g++-9
# CC=/usr/bin/gcc-9
# Build flags
DEBUG_FLAGS=-g -DDEBUG -fsanitize=address
RELEASE_FLAGS=-O3 -DNDEBUG -march=native
PROFILE_FLAGS=-pg -O2
# Include paths (additional to standard ones)
# EXTRA_INCLUDES=-I/usr/local/include
# Library paths (additional to standard ones)
# EXTRA_LIBS=-L/usr/local/lib
# Platform-specific settings
ifeq ($(OS),Windows_NT)
    # Windows-specific flags
    PLATFORM_FLAGS=-D_WIN32_WINNT=0x0601 -DWIN32_LEAN_AND_MEAN
    EXTRA_LIBS += -lws2_32
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        # Linux-specific flags
        PLATFORM_FLAGS=-D_GNU_SOURCE
        EXTRA_LIBS += -lpthread
    endif
    ifeq ($(UNAME_S),Darwin)
        # macOS-specific flags
        PLATFORM_FLAGS=-D_DARWIN_C_SOURCE
        EXTRA_LIBS += -framework CoreFoundation
    endif
endif
# Testing framework settings
CATCH_VERSION=2.13.7
TEST_FLAGS=-DCATCH_CONFIG_FAST_COMPILE
# Coverage settings
COVERAGE_FLAGS=--coverage -fprofile-arcs -ftest-coverage
COVERAGE_LIBS=-lgcov
# Sanitizer options
ASAN_FLAGS=-fsanitize=address -fno-omit-frame-pointer
TSAN_FLAGS=-fsanitize=thread
UBSAN_FLAGS=-fsanitize=undefined
# Documentation settings
DOXYGEN_CONFIG=Doxyfile
DOC_OUTPUT_DIR=docs
# Installation settings
DEFAULT_PREFIX=/usr/local
HEADER_INSTALL_DIR=$(PREFIX)/include/grf
LIB_INSTALL_DIR=$(PREFIX)/lib
