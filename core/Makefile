
# GRF C++ Library Makefile
# ======================================================================================
# Build configuration
# ======================================================================================
# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -pthread
CXXFLAGS_DEBUG = -std=c++11 -Wall -g -pthread -DDEBUG
CXXFLAGS_FAST = -std=c++11 -Wall -O3 -pthread -DNDEBUG -march=native -flto
INCLUDES = -Isrc -Itest -Ithird_party
# Build optimization flags
FAST_BUILD_FLAGS = -pipe -ffast-math
CCACHE_FLAGS = 
# Check if ccache is available for faster compilation
ifneq ($(shell which ccache 2>/dev/null),)
    CXX := ccache $(CXX)
    CCACHE_FLAGS = -DUSING_CCACHE
endif
# Add build optimization flags
CXXFLAGS += $(FAST_BUILD_FLAGS) $(CCACHE_FLAGS)
CXXFLAGS_DEBUG += $(FAST_BUILD_FLAGS) $(CCACHE_FLAGS)
CXXFLAGS_FAST += $(FAST_BUILD_FLAGS) $(CCACHE_FLAGS)
# Directories
SRCDIR = src
TESTDIR = test
BUILDDIR = build
OBJDIR = $(BUILDDIR)/obj
TESTBUILDDIR = $(BUILDDIR)/test
THIRD_PARTY_DIR = third_party
# Find source files
SRCFILES = $(shell find $(SRCDIR) -name "*.cpp")
TESTFILES = $(shell find $(TESTDIR) -name "*.cpp")
TESTHEADERFILES = $(shell find $(TESTDIR) -name "*.h")
SRCHEADERFILES = $(shell find $(SRCDIR) -name "*.h")
THIRD_PARTY_HEADERS = $(shell find $(THIRD_PARTY_DIR) -name "*.h*")
# Object files
SRCOBJS = $(SRCFILES:$(SRCDIR)/%.cpp=$(OBJDIR)/src/%.o)
TESTOBJS = $(TESTFILES:$(TESTDIR)/%.cpp=$(OBJDIR)/test/%.o)
# Dependency files for incremental builds
SRCDEPS = $(SRCOBJS:.o=.d)
TESTDEPS = $(TESTOBJS:.o=.d)
# Targets
TEST_TARGET = $(BUILDDIR)/grf_test
LIBRARY_TARGET = $(BUILDDIR)/libgrf.a
# Default target
.PHONY: all clean test debug release library help directories test_on_mac
.DEFAULT_GOAL := help
# ======================================================================================
# Build targets
# ======================================================================================
all: directories library test ## Build library and tests
library: directories $(LIBRARY_TARGET) ## Build static library only
test: directories $(TEST_TARGET) ## Build and run tests
fast: CXXFLAGS = $(CXXFLAGS_FAST)
fast: directories $(LIBRARY_TARGET) ## Build optimized version quickly
debug: CXXFLAGS = $(CXXFLAGS_DEBUG)
debug: directories $(TEST_TARGET) ## Build debug version
release: CXXFLAGS += -DNDEBUG
release: directories $(LIBRARY_TARGET) ## Build optimized release version
# Quick rebuild - only recompile changed files (no clean)
quick: $(TEST_TARGET) ## Quick incremental build and test
# Ultra-fast build - minimal dependencies checking
ultra-fast: 
	@echo "Ultra-fast build mode..."
	$(MAKE) -j$(shell nproc 2>/dev/null || echo 4) CXXFLAGS="$(CXXFLAGS_FAST)" $(TEST_TARGET)
# ======================================================================================
# Directory creation
# ======================================================================================
directories:
	@echo "Creating build directories..."
	@mkdir -p $(OBJDIR)/src/analysis 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/commons 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/forest 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/prediction 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/prediction/collector 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/relabeling 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/sampling 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/splitting 2>/dev/null || true
	@mkdir -p $(OBJDIR)/src/tree 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/analysis 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/commons 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/forest 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/prediction 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/relabeling 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/sampling 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/splitting 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/tree 2>/dev/null || true
	@mkdir -p $(OBJDIR)/test/utilities 2>/dev/null || true
# ======================================================================================
# Library target
# ======================================================================================
$(LIBRARY_TARGET): $(SRCOBJS) | directories
	@echo "Creating static library..."
	ar rcs $@ $^
	@echo "Library created: $@"
# ======================================================================================
# Test target
# ======================================================================================
$(TEST_TARGET): $(SRCOBJS) $(TESTOBJS) | directories
	@echo "Linking test executable..."
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Test executable created: $@"
# ======================================================================================
# Object file compilation
# ======================================================================================
# Source object files with dependency generation
$(OBJDIR)/src/%.o: $(SRCDIR)/%.cpp | directories
	@echo "Compiling source: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@
# Test object files with dependency generation
$(OBJDIR)/test/%.o: $(TESTDIR)/%.cpp | directories
	@echo "Compiling test: $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@
# Include dependency files for proper incremental builds
-include $(SRCDEPS)
-include $(TESTDEPS)
# ======================================================================================
# Test execution targets
# ======================================================================================
run-tests: $(TEST_TARGET) ## Build and run all tests
	@echo "Running tests..."
	./$(TEST_TARGET)
run-tests-verbose: $(TEST_TARGET) ## Build and run tests with verbose output
	@echo "Running tests with verbose output..."
	./$(TEST_TARGET) --success --durations yes
run-tests-specific: $(TEST_TARGET) ## Build and run specific test (use TEST=pattern)
	@echo "Running specific tests matching: $(TEST)"
	./$(TEST_TARGET) "$(TEST)"
run-tests-junit: $(TEST_TARGET) ## Build and run tests with JUnit XML output
	@echo "Running tests with JUnit output..."
	./$(TEST_TARGET) -r junit -o $(BUILDDIR)/test_results.xml
run-tests-list: $(TEST_TARGET) ## List all available tests
	@echo "Listing all available tests..."
	./$(TEST_TARGET) --list-tests
run-tests-tags: $(TEST_TARGET) ## List all available test tags
	@echo "Listing all available test tags..."
	./$(TEST_TARGET) --list-tags
run-tests-compact: $(TEST_TARGET) ## Run tests with compact output
	@echo "Running tests with compact output..."
	./$(TEST_TARGET) -r compact
run-tests-xml: $(TEST_TARGET) ## Run tests with XML output
	@echo "Running tests with XML output..."
	./$(TEST_TARGET) -r xml -o $(BUILDDIR)/test_results.xml

test_on_mac: $(TEST_TARGET) ## Run tests optimized for macOS (excludes problematic NaN tests)
	@echo "Running macOS-optimized test suite..."
	@echo "Note: Excluding NaN tests due to performance issues on macOS"
	./$(TEST_TARGET) '~[NaN]' -r compact

# ======================================================================================
# Analysis and debugging targets
# ======================================================================================
analyze: ## Run static analysis (requires cppcheck)
	@echo "Running static analysis..."
	cppcheck --enable=all --std=c++11 $(INCLUDES) $(SRCDIR)/
valgrind: $(TEST_TARGET) ## Run tests under valgrind (Linux only)
	@echo "Running tests under valgrind..."
	valgrind --leak-check=full --track-origins=yes ./$(TEST_TARGET)
gdb: $(TEST_TARGET) ## Run tests in GDB debugger
	@echo "Starting GDB debugger..."
	gdb ./$(TEST_TARGET)
# ======================================================================================
# Code coverage targets
# ======================================================================================
coverage: CXXFLAGS += -fprofile-arcs -ftest-coverage
coverage: LDFLAGS += -lgcov
coverage: clean $(TEST_TARGET) ## Build with coverage and run tests
	@echo "Running tests for coverage..."
	./$(TEST_TARGET)
	@echo "Generating coverage report..."
	gcov $(SRCFILES)
	@echo "Coverage files generated (.gcov)"
coverage-html: coverage ## Generate HTML coverage report (requires lcov)
	@echo "Generating HTML coverage report..."
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory $(BUILDDIR)/coverage
	@echo "HTML coverage report generated in $(BUILDDIR)/coverage"
# ======================================================================================
# Documentation targets
# ======================================================================================
docs: ## Generate documentation (requires doxygen)
	@echo "Generating documentation..."
	doxygen Doxyfile
	@echo "Documentation generated in docs/"
# ======================================================================================
# Benchmarking targets
# ======================================================================================
benchmark: CXXFLAGS += -DBENCHMARK -O3
benchmark: $(TEST_TARGET) ## Build optimized version for benchmarking
	@echo "Running benchmarks..."
	./$(TEST_TARGET) "[benchmark]"
# ======================================================================================
# Installation targets
# ======================================================================================
PREFIX ?= /usr/local
INCDIR = $(PREFIX)/include/grf
LIBDIR = $(PREFIX)/lib
install: $(LIBRARY_TARGET) ## Install library and headers
	@echo "Installing library to $(LIBDIR)..."
	install -d $(LIBDIR)
	install -m 644 $(LIBRARY_TARGET) $(LIBDIR)/
	@echo "Installing headers to $(INCDIR)..."
	install -d $(INCDIR)
	find $(SRCDIR) -name "*.h" -exec install -m 644 {} $(INCDIR)/ \;
	@echo "Installation complete"
uninstall: ## Uninstall library and headers
	@echo "Uninstalling..."
	rm -f $(LIBDIR)/libgrf.a
	rm -rf $(INCDIR)
	@echo "Uninstall complete"
# ======================================================================================
# Utility targets
# ======================================================================================
clean: ## Clean build artifacts
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILDDIR)
	rm -f *.gcov *.gcda *.gcno coverage.info
	@echo "Clean complete"
clean-all: clean ## Clean everything including documentation
	rm -rf docs/
distclean: clean-all ## Clean everything for distribution
info: ## Show build information
	@echo "GRF C++ Library Build Information"
	@echo "=================================="
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Includes: $(INCLUDES)"
	@echo "Source files: $(words $(SRCFILES))"
	@echo "Test files: $(words $(TESTFILES))"
	@echo "Build directory: $(BUILDDIR)"
	@echo "Library target: $(LIBRARY_TARGET)"
	@echo "Test target: $(TEST_TARGET)"
list-tests: $(TEST_TARGET) ## List available tests
	@echo "Available tests:"
	./$(TEST_TARGET) --list-tests
# ======================================================================================
# Development targets
# ======================================================================================
format: ## Format code (requires clang-format)
	@echo "Formatting source code..."
	find $(SRCDIR) $(TESTDIR) -name "*.cpp" -o -name "*.h" | xargs clang-format -i
	@echo "Code formatted"
lint: ## Lint code (requires cpplint)
	@echo "Linting code..."
	find $(SRCDIR) $(TESTDIR) -name "*.cpp" -o -name "*.h" | xargs cpplint
	@echo "Linting complete"
check-deps: ## Check for required dependencies
	@echo "Checking dependencies..."
	@which $(CXX) > /dev/null || echo "ERROR: $(CXX) not found"
	@which ar > /dev/null || echo "ERROR: ar not found"
	@echo "Dependencies check complete"
# ======================================================================================
# Help target
# ======================================================================================
help: ## Show this help message
	@echo "GRF C++ Library Makefile"
	@echo "========================"
	@echo ""
	@echo "Available targets:"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@echo ""
	@echo "Fast build options:"
	@echo "  make -j             # Use all CPU cores"
	@echo "  make -j8            # Use 8 cores specifically"
	@echo "  make quick          # Quick incremental build"
	@echo "  make fast           # Optimized fast build"
	@echo "  make ultra-fast     # Maximum speed build"
	@echo ""
	@echo "Examples:"
	@echo "  make -j test                     # Parallel build and test"
	@echo "  make quick                       # Quick incremental build"
	@echo "  make fast                        # Fast optimized build"
	@echo "  make run-tests-specific TEST='*Regression*'  # Run specific tests"
	@echo "  make test_on_mac                 # Run tests optimized for macOS"
	@echo "  make coverage                    # Build with coverage"
	@echo "  make install PREFIX=/opt/grf     # Install to custom location"
	@echo ""
	@echo "Build optimizations:"
	@echo "  - Automatic dependency tracking for incremental builds"
	@echo "  - ccache support (auto-detected)"
	@echo "  - Parallel compilation with -j flag"
	@echo "  - Fast math and link-time optimization available"
# ======================================================================================
# Platform-specific adjustments
# ======================================================================================
# Windows (MinGW/MSYS2)
ifeq ($(OS),Windows_NT)
    TEST_TARGET = $(BUILDDIR)/grf_test.exe
    LIBRARY_TARGET = $(BUILDDIR)/grf.lib
    CXXFLAGS += -D_WIN32_WINNT=0x0601
endif
# macOS
ifeq ($(shell uname),Darwin)
    CXXFLAGS += -stdlib=libc++
endif
# ======================================================================================
# Dependencies tracking
# ======================================================================================
-include $(SRCOBJS:.o=.d)
-include $(TESTOBJS:.o=.d)
%.d: %.cpp
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM -MT '$(@:.d=.o)' $< > $@
