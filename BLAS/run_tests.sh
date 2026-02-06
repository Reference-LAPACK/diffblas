#!/bin/bash

# Test script for running all differentiated BLAS function tests
# This script tests the following modes: forward (d), vector forward (dv), reverse (b), vector reverse (bv)

# Don't use set -e here - we want to always print the summary even if there are errors
# Individual commands will handle their own error checking

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Counters for Forward Mode tests (scalar)
FWD_SCALAR_TOTAL=0
FWD_SCALAR_MACHINE_PRECISION=0
FWD_SCALAR_ACCEPTABLE=0
FWD_SCALAR_OUTSIDE_TOLERANCE=0
FWD_SCALAR_EXECUTION_FAILED=0
FWD_SCALAR_SKIPPED=0

# Counters for Forward Mode tests (vector)
FWD_VECTOR_TOTAL=0
FWD_VECTOR_MACHINE_PRECISION=0
FWD_VECTOR_ACCEPTABLE=0
FWD_VECTOR_OUTSIDE_TOLERANCE=0
FWD_VECTOR_EXECUTION_FAILED=0
FWD_VECTOR_SKIPPED=0

# Legacy combined counters for backward compatibility
FWD_TOTAL=0
FWD_MACHINE_PRECISION=0
FWD_ACCEPTABLE=0
FWD_OUTSIDE_TOLERANCE=0
FWD_EXECUTION_FAILED=0
FWD_SKIPPED=0

# Counters for Reverse Mode tests (scalar)
REV_SCALAR_TOTAL=0
REV_SCALAR_MACHINE_PRECISION=0
REV_SCALAR_ACCEPTABLE=0
REV_SCALAR_OUTSIDE_TOLERANCE=0
REV_SCALAR_EXECUTION_FAILED=0
REV_SCALAR_SKIPPED=0

# Counters for Reverse Mode tests (vector)
REV_VECTOR_TOTAL=0
REV_VECTOR_MACHINE_PRECISION=0
REV_VECTOR_ACCEPTABLE=0
REV_VECTOR_OUTSIDE_TOLERANCE=0
REV_VECTOR_EXECUTION_FAILED=0
REV_VECTOR_SKIPPED=0

# Legacy combined counters for backward compatibility
REV_TOTAL=0
REV_MACHINE_PRECISION=0
REV_ACCEPTABLE=0
REV_OUTSIDE_TOLERANCE=0
REV_EXECUTION_FAILED=0
REV_SKIPPED=0

# Overall counters
TOTAL_TESTS=0
TAPENADE_FAILED_TESTS=0
UNABLE_TO_PROCESS_TESTS=0

# Arrays to store results by mode (scalar forward)
FWD_SCALAR_MACHINE_PRECISION_LIST=()
FWD_SCALAR_ACCEPTABLE_LIST=()
FWD_SCALAR_OUTSIDE_TOLERANCE_LIST=()
FWD_SCALAR_EXECUTION_FAILED_LIST=()
FWD_SCALAR_SKIPPED_LIST=()

# Arrays to store results by mode (vector forward)
FWD_VECTOR_MACHINE_PRECISION_LIST=()
FWD_VECTOR_ACCEPTABLE_LIST=()
FWD_VECTOR_OUTSIDE_TOLERANCE_LIST=()
FWD_VECTOR_EXECUTION_FAILED_LIST=()
FWD_VECTOR_SKIPPED_LIST=()

# Arrays to store results by mode (scalar reverse)
REV_SCALAR_MACHINE_PRECISION_LIST=()
REV_SCALAR_ACCEPTABLE_LIST=()
REV_SCALAR_OUTSIDE_TOLERANCE_LIST=()
REV_SCALAR_EXECUTION_FAILED_LIST=()
REV_SCALAR_SKIPPED_LIST=()

# Arrays to store results by mode (vector reverse)
REV_VECTOR_MACHINE_PRECISION_LIST=()
REV_VECTOR_ACCEPTABLE_LIST=()
REV_VECTOR_OUTSIDE_TOLERANCE_LIST=()
REV_VECTOR_EXECUTION_FAILED_LIST=()
REV_VECTOR_SKIPPED_LIST=()

# Legacy combined arrays for backward compatibility
FWD_MACHINE_PRECISION_LIST=()
FWD_ACCEPTABLE_LIST=()
FWD_OUTSIDE_TOLERANCE_LIST=()
FWD_EXECUTION_FAILED_LIST=()
FWD_SKIPPED_LIST=()

REV_MACHINE_PRECISION_LIST=()
REV_ACCEPTABLE_LIST=()
REV_OUTSIDE_TOLERANCE_LIST=()
REV_EXECUTION_FAILED_LIST=()
REV_SKIPPED_LIST=()

TAPENADE_FAILED_LIST=()
UNABLE_TO_PROCESS_LIST=()

# Function to print colored status
print_status() {
    local status=$1
    local message=$2
    case $status in
        "MACHINE_PRECISION")
            echo -e "${GREEN}[MACHINE_PRECISION]${NC} $message"
            ;;
        "ACCEPTABLE")
            echo -e "${GREEN}[ACCEPTABLE]${NC} $message"
            ;;
        "OUTSIDE_TOLERANCE")
            echo -e "${YELLOW}[OUTSIDE_TOLERANCE]${NC} $message"
            ;;
        "EXECUTION_FAILED")
            echo -e "${RED}[EXECUTION_FAILED]${NC} $message"
            ;;
        "SKIPPED")
            echo -e "${CYAN}[SKIPPED]${NC} $message"
            ;;
        "TAPENADE_FAILED")
            echo -e "${MAGENTA}[TAPENADE_FAILED]${NC} $message"
            ;;
        "UNABLE_TO_PROCESS")
            echo -e "${RED}[UNABLE_TO_PROCESS]${NC} $message"
            ;;
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        *)
            echo -e "[$status] $message"
            ;;
    esac
}

# Timeout for test execution (in seconds)
TEST_TIMEOUT=10

# Function to safely run a test with timeout and signal handling
safe_run_test() {
    local test_executable=$1
    local output_file=$2
    
    # Use timeout to prevent hanging tests
    timeout ${TEST_TIMEOUT}s ./"$test_executable" > "$output_file" 2>&1
    local exit_code=$?
    
    # Check if timeout killed the process (exit code 124)
    if [ $exit_code -eq 124 ]; then
        echo "Test timed out after ${TEST_TIMEOUT}s" >> "$output_file"
        return 124
    fi
    
    # Check if the test crashed (empty output file usually indicates a crash)
    if [ ! -s "$output_file" ]; then
        echo "Test crashed or produced no output" >> "$output_file"
        return 1
    fi
    
    return $exit_code
}

# Function to run a single test (forward or reverse mode)
run_single_test() {
    local test_executable=$1
    local test_name=$2
    local mode=$3  # "FWD", "FWD_VEC", "REV", or "REV_VEC"
    local output_file="test_${mode}_output.log"
    
    # Determine if this is a forward or reverse mode for counter purposes
    local is_forward=false
    local is_forward_scalar=false
    local is_forward_vector=false
    local is_reverse_scalar=false
    local is_reverse_vector=false
    if [[ "$mode" == "FWD"* ]]; then
        is_forward=true
        if [[ "$mode" == "FWD" ]]; then
            is_forward_scalar=true
        elif [[ "$mode" == "FWD_VEC" ]]; then
            is_forward_vector=true
        fi
    elif [[ "$mode" == "REV"* ]]; then
        if [[ "$mode" == "REV" ]]; then
            is_reverse_scalar=true
        elif [[ "$mode" == "REV_VEC" ]]; then
            is_reverse_vector=true
        fi
    fi
    
    if [ ! -f "$test_executable" ]; then
        if [ "$is_forward_scalar" = "true" ]; then
            FWD_SCALAR_SKIPPED=$((FWD_SCALAR_SKIPPED + 1))
            FWD_SCALAR_SKIPPED_LIST+=("$test_name")
            FWD_SKIPPED=$((FWD_SKIPPED + 1))
            FWD_SKIPPED_LIST+=("$test_name")
        elif [ "$is_forward_vector" = "true" ]; then
            FWD_VECTOR_SKIPPED=$((FWD_VECTOR_SKIPPED + 1))
            FWD_VECTOR_SKIPPED_LIST+=("$test_name")
            FWD_SKIPPED=$((FWD_SKIPPED + 1))
            FWD_SKIPPED_LIST+=("$test_name")
        elif [ "$is_reverse_scalar" = "true" ]; then
            REV_SCALAR_SKIPPED=$((REV_SCALAR_SKIPPED + 1))
            REV_SCALAR_SKIPPED_LIST+=("$test_name")
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        elif [ "$is_reverse_vector" = "true" ]; then
            REV_VECTOR_SKIPPED=$((REV_VECTOR_SKIPPED + 1))
            REV_VECTOR_SKIPPED_LIST+=("$test_name")
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        else
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        fi
        print_status "SKIPPED" "$test_name ($mode): Test executable not found"
        return
    fi
    
    if [ ! -x "$test_executable" ]; then
        if [ "$is_forward_scalar" = "true" ]; then
            FWD_SCALAR_SKIPPED=$((FWD_SCALAR_SKIPPED + 1))
            FWD_SCALAR_SKIPPED_LIST+=("$test_name")
            FWD_SKIPPED=$((FWD_SKIPPED + 1))
            FWD_SKIPPED_LIST+=("$test_name")
        elif [ "$is_forward_vector" = "true" ]; then
            FWD_VECTOR_SKIPPED=$((FWD_VECTOR_SKIPPED + 1))
            FWD_VECTOR_SKIPPED_LIST+=("$test_name")
            FWD_SKIPPED=$((FWD_SKIPPED + 1))
            FWD_SKIPPED_LIST+=("$test_name")
        elif [ "$is_reverse_scalar" = "true" ]; then
            REV_SCALAR_SKIPPED=$((REV_SCALAR_SKIPPED + 1))
            REV_SCALAR_SKIPPED_LIST+=("$test_name")
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        elif [ "$is_reverse_vector" = "true" ]; then
            REV_VECTOR_SKIPPED=$((REV_VECTOR_SKIPPED + 1))
            REV_VECTOR_SKIPPED_LIST+=("$test_name")
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        else
            REV_SKIPPED=$((REV_SKIPPED + 1))
            REV_SKIPPED_LIST+=("$test_name")
        fi
        print_status "SKIPPED" "$test_name ($mode): Test executable exists but is not executable"
        return
    fi
    
    # Run the test safely
    safe_run_test "$test_executable" "$output_file"
    local exit_code=$?
    
    # Check for execution failure patterns
    local has_execution_failures=false
    if grep -q "Segmentation fault\\|Aborted\\|Floating point exception\\|Test timed out\\|had an illegal value\\|error while loading shared libraries\\|cannot open shared object file" "$output_file" 2>/dev/null; then
        has_execution_failures=true
    fi
    
    # Check for derivative tolerance patterns (for both forward and reverse mode)
    local has_machine_precision=false
    local has_acceptable=false
    local has_outside_tolerance=false
    
    if grep -q "FAIL: Large errors detected" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    elif grep -q "PASS: Derivatives are accurate to machine precision" "$output_file" 2>/dev/null; then
        has_machine_precision=true
    elif grep -q "PASS: Vector derivatives are accurate to machine precision" "$output_file" 2>/dev/null; then
        has_machine_precision=true
    elif grep -q "PASS: Derivatives are reasonably accurate" "$output_file" 2>/dev/null; then
        has_acceptable=true
    elif grep -q "PASS: Vector derivatives are reasonably accurate" "$output_file" 2>/dev/null; then
        has_acceptable=true
    elif grep -q "WARNING: Derivatives may have significant errors" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    elif grep -q "WARNING: Vector derivatives may have significant errors" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    fi
    
    # Determine test result category and update counters
    if [ $exit_code -eq 0 ] && [ "$has_execution_failures" = false ]; then
        if [ "$has_machine_precision" = true ]; then
            if [ "$is_forward_scalar" = "true" ]; then
                FWD_SCALAR_MACHINE_PRECISION=$((FWD_SCALAR_MACHINE_PRECISION + 1))
                FWD_SCALAR_MACHINE_PRECISION_LIST+=("$test_name")
                FWD_MACHINE_PRECISION=$((FWD_MACHINE_PRECISION + 1))
                FWD_MACHINE_PRECISION_LIST+=("$test_name")
            elif [ "$is_forward_vector" = "true" ]; then
                FWD_VECTOR_MACHINE_PRECISION=$((FWD_VECTOR_MACHINE_PRECISION + 1))
                FWD_VECTOR_MACHINE_PRECISION_LIST+=("$test_name")
                FWD_MACHINE_PRECISION=$((FWD_MACHINE_PRECISION + 1))
                FWD_MACHINE_PRECISION_LIST+=("$test_name")
            elif [ "$is_reverse_scalar" = "true" ]; then
                REV_SCALAR_MACHINE_PRECISION=$((REV_SCALAR_MACHINE_PRECISION + 1))
                REV_SCALAR_MACHINE_PRECISION_LIST+=("$test_name")
                REV_MACHINE_PRECISION=$((REV_MACHINE_PRECISION + 1))
                REV_MACHINE_PRECISION_LIST+=("$test_name")
            elif [ "$is_reverse_vector" = "true" ]; then
                REV_VECTOR_MACHINE_PRECISION=$((REV_VECTOR_MACHINE_PRECISION + 1))
                REV_VECTOR_MACHINE_PRECISION_LIST+=("$test_name")
                REV_MACHINE_PRECISION=$((REV_MACHINE_PRECISION + 1))
                REV_MACHINE_PRECISION_LIST+=("$test_name")
            else
                REV_MACHINE_PRECISION=$((REV_MACHINE_PRECISION + 1))
                REV_MACHINE_PRECISION_LIST+=("$test_name")
            fi
            print_status "MACHINE_PRECISION" "$test_name ($mode): Derivatives match to machine precision"
        elif [ "$has_acceptable" = true ]; then
            if [ "$is_forward_scalar" = "true" ]; then
                FWD_SCALAR_ACCEPTABLE=$((FWD_SCALAR_ACCEPTABLE + 1))
                FWD_SCALAR_ACCEPTABLE_LIST+=("$test_name")
                FWD_ACCEPTABLE=$((FWD_ACCEPTABLE + 1))
                FWD_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_forward_vector" = "true" ]; then
                FWD_VECTOR_ACCEPTABLE=$((FWD_VECTOR_ACCEPTABLE + 1))
                FWD_VECTOR_ACCEPTABLE_LIST+=("$test_name")
                FWD_ACCEPTABLE=$((FWD_ACCEPTABLE + 1))
                FWD_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_reverse_scalar" = "true" ]; then
                REV_SCALAR_ACCEPTABLE=$((REV_SCALAR_ACCEPTABLE + 1))
                REV_SCALAR_ACCEPTABLE_LIST+=("$test_name")
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_reverse_vector" = "true" ]; then
                REV_VECTOR_ACCEPTABLE=$((REV_VECTOR_ACCEPTABLE + 1))
                REV_VECTOR_ACCEPTABLE_LIST+=("$test_name")
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            else
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            fi
            print_status "ACCEPTABLE" "$test_name ($mode): Derivatives are acceptable"
        elif [ "$has_outside_tolerance" = true ]; then
            if [ "$is_forward_scalar" = "true" ]; then
                FWD_SCALAR_OUTSIDE_TOLERANCE=$((FWD_SCALAR_OUTSIDE_TOLERANCE + 1))
                FWD_SCALAR_OUTSIDE_TOLERANCE_LIST+=("$test_name")
                FWD_OUTSIDE_TOLERANCE=$((FWD_OUTSIDE_TOLERANCE + 1))
                FWD_OUTSIDE_TOLERANCE_LIST+=("$test_name")
            elif [ "$is_forward_vector" = "true" ]; then
                FWD_VECTOR_OUTSIDE_TOLERANCE=$((FWD_VECTOR_OUTSIDE_TOLERANCE + 1))
                FWD_VECTOR_OUTSIDE_TOLERANCE_LIST+=("$test_name")
                FWD_OUTSIDE_TOLERANCE=$((FWD_OUTSIDE_TOLERANCE + 1))
                FWD_OUTSIDE_TOLERANCE_LIST+=("$test_name")
            elif [ "$is_reverse_scalar" = "true" ]; then
                REV_SCALAR_OUTSIDE_TOLERANCE=$((REV_SCALAR_OUTSIDE_TOLERANCE + 1))
                REV_SCALAR_OUTSIDE_TOLERANCE_LIST+=("$test_name")
                REV_OUTSIDE_TOLERANCE=$((REV_OUTSIDE_TOLERANCE + 1))
                REV_OUTSIDE_TOLERANCE_LIST+=("$test_name")
            elif [ "$is_reverse_vector" = "true" ]; then
                REV_VECTOR_OUTSIDE_TOLERANCE=$((REV_VECTOR_OUTSIDE_TOLERANCE + 1))
                REV_VECTOR_OUTSIDE_TOLERANCE_LIST+=("$test_name")
                REV_OUTSIDE_TOLERANCE=$((REV_OUTSIDE_TOLERANCE + 1))
                REV_OUTSIDE_TOLERANCE_LIST+=("$test_name")
            else
                REV_OUTSIDE_TOLERANCE=$((REV_OUTSIDE_TOLERANCE + 1))
                REV_OUTSIDE_TOLERANCE_LIST+=("$test_name")
            fi
            print_status "OUTSIDE_TOLERANCE" "$test_name ($mode): Code runs but derivatives outside acceptable tolerance"
        else
            # Test completed but no clear derivative status - treat as acceptable
            if [ "$is_forward_scalar" = "true" ]; then
                FWD_SCALAR_ACCEPTABLE=$((FWD_SCALAR_ACCEPTABLE + 1))
                FWD_SCALAR_ACCEPTABLE_LIST+=("$test_name")
                FWD_ACCEPTABLE=$((FWD_ACCEPTABLE + 1))
                FWD_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_forward_vector" = "true" ]; then
                FWD_VECTOR_ACCEPTABLE=$((FWD_VECTOR_ACCEPTABLE + 1))
                FWD_VECTOR_ACCEPTABLE_LIST+=("$test_name")
                FWD_ACCEPTABLE=$((FWD_ACCEPTABLE + 1))
                FWD_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_reverse_scalar" = "true" ]; then
                REV_SCALAR_ACCEPTABLE=$((REV_SCALAR_ACCEPTABLE + 1))
                REV_SCALAR_ACCEPTABLE_LIST+=("$test_name")
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            elif [ "$is_reverse_vector" = "true" ]; then
                REV_VECTOR_ACCEPTABLE=$((REV_VECTOR_ACCEPTABLE + 1))
                REV_VECTOR_ACCEPTABLE_LIST+=("$test_name")
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            else
                REV_ACCEPTABLE=$((REV_ACCEPTABLE + 1))
                REV_ACCEPTABLE_LIST+=("$test_name")
            fi
            print_status "ACCEPTABLE" "$test_name ($mode): Test completed successfully"
        fi
        echo "  Last line of output:"
        tail -1 "$output_file" | sed 's/^/    /'
    elif [ "$has_execution_failures" = true ]; then
        if [ "$is_forward_scalar" = "true" ]; then
            FWD_SCALAR_EXECUTION_FAILED=$((FWD_SCALAR_EXECUTION_FAILED + 1))
            FWD_SCALAR_EXECUTION_FAILED_LIST+=("$test_name")
            FWD_EXECUTION_FAILED=$((FWD_EXECUTION_FAILED + 1))
            FWD_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_forward_vector" = "true" ]; then
            FWD_VECTOR_EXECUTION_FAILED=$((FWD_VECTOR_EXECUTION_FAILED + 1))
            FWD_VECTOR_EXECUTION_FAILED_LIST+=("$test_name")
            FWD_EXECUTION_FAILED=$((FWD_EXECUTION_FAILED + 1))
            FWD_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_reverse_scalar" = "true" ]; then
            REV_SCALAR_EXECUTION_FAILED=$((REV_SCALAR_EXECUTION_FAILED + 1))
            REV_SCALAR_EXECUTION_FAILED_LIST+=("$test_name")
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_reverse_vector" = "true" ]; then
            REV_VECTOR_EXECUTION_FAILED=$((REV_VECTOR_EXECUTION_FAILED + 1))
            REV_VECTOR_EXECUTION_FAILED_LIST+=("$test_name")
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        else
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        fi
        print_status "EXECUTION_FAILED" "$test_name ($mode): Code fails to complete execution"
        echo "  Error output:"
        grep -E "Segmentation fault|Aborted|Floating point exception|Test timed out|had an illegal value|error while loading shared libraries|cannot open shared object file" "$output_file" | head -3 | sed 's/^/    /'
    else
        # Test completed with non-zero exit code but no clear failure patterns
        if [ "$is_forward_scalar" = "true" ]; then
            FWD_SCALAR_EXECUTION_FAILED=$((FWD_SCALAR_EXECUTION_FAILED + 1))
            FWD_SCALAR_EXECUTION_FAILED_LIST+=("$test_name")
            FWD_EXECUTION_FAILED=$((FWD_EXECUTION_FAILED + 1))
            FWD_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_forward_vector" = "true" ]; then
            FWD_VECTOR_EXECUTION_FAILED=$((FWD_VECTOR_EXECUTION_FAILED + 1))
            FWD_VECTOR_EXECUTION_FAILED_LIST+=("$test_name")
            FWD_EXECUTION_FAILED=$((FWD_EXECUTION_FAILED + 1))
            FWD_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_reverse_scalar" = "true" ]; then
            REV_SCALAR_EXECUTION_FAILED=$((REV_SCALAR_EXECUTION_FAILED + 1))
            REV_SCALAR_EXECUTION_FAILED_LIST+=("$test_name")
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        elif [ "$is_reverse_vector" = "true" ]; then
            REV_VECTOR_EXECUTION_FAILED=$((REV_VECTOR_EXECUTION_FAILED + 1))
            REV_VECTOR_EXECUTION_FAILED_LIST+=("$test_name")
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        else
            REV_EXECUTION_FAILED=$((REV_EXECUTION_FAILED + 1))
            REV_EXECUTION_FAILED_LIST+=("$test_name")
        fi
        print_status "EXECUTION_FAILED" "$test_name ($mode): Test failed with exit code $exit_code"
        echo "  Last line of output:"
        tail -1 "$output_file" | sed 's/^/    /'
    fi
}

# Configuration: Which modes to test
RUN_D=true
RUN_DV=true
RUN_B=true
RUN_BV=true
FLAT_MODE=true

# Function to run test for a single function (flat mode)
run_test_for_func() {
    local funcname=$1
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    print_status "INFO" "Testing $funcname"
    
    # Check if this function is marked as "Unable to process"
    if [ -f ".unable_to_process_${funcname}" ]; then
        UNABLE_TO_PROCESS_TESTS=$((UNABLE_TO_PROCESS_TESTS + 1))
        UNABLE_TO_PROCESS_LIST+=("$funcname")
        print_status "UNABLE_TO_PROCESS" "$funcname: Insufficient parameter documentation"
        echo ""
        return
    fi
    
    # Check if Tapenade failed for this function (flat mode: files in src/ dir)
    # Tapenade preserves the source file extension (.f or .f90), so check both
    local has_any_mode=false
    if [ "$RUN_D" = "true" ] && ([ -f "src/${funcname}_d.f" ] || [ -f "src/${funcname}_d.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_DV" = "true" ] && ([ -f "src/${funcname}_dv.f" ] || [ -f "src/${funcname}_dv.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_B" = "true" ] && ([ -f "src/${funcname}_b.f" ] || [ -f "src/${funcname}_b.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_BV" = "true" ] && ([ -f "src/${funcname}_bv.f" ] || [ -f "src/${funcname}_bv.f90" ]); then
        has_any_mode=true
    fi
    
    if [ "$has_any_mode" = "false" ]; then
        TAPENADE_FAILED_TESTS=$((TAPENADE_FAILED_TESTS + 1))
        TAPENADE_FAILED_LIST+=("$funcname")
        print_status "TAPENADE_FAILED" "$funcname: Tapenade fails to differentiate the code"
        echo ""
        return
    fi
    
    # Run scalar forward mode test (flat mode: test_funcname in build/ dir)
    if [ "$RUN_D" = "true" ]; then
        FWD_SCALAR_TOTAL=$((FWD_SCALAR_TOTAL + 1))
        FWD_TOTAL=$((FWD_TOTAL + 1))
        run_single_test "build/test_$funcname" "$funcname" "FWD"
    fi
    
    # Run vector forward mode test  
    if [ "$RUN_DV" = "true" ]; then
        FWD_VECTOR_TOTAL=$((FWD_VECTOR_TOTAL + 1))
        FWD_TOTAL=$((FWD_TOTAL + 1))
        run_single_test "build/test_${funcname}_vector_forward" "$funcname" "FWD_VEC"
    fi
    
    # Run scalar reverse mode test
    if [ "$RUN_B" = "true" ]; then
        REV_SCALAR_TOTAL=$((REV_SCALAR_TOTAL + 1))
        REV_TOTAL=$((REV_TOTAL + 1))
        run_single_test "build/test_${funcname}_reverse" "$funcname" "REV"
    fi
    
    # Run vector reverse mode test
    if [ "$RUN_BV" = "true" ]; then
        REV_VECTOR_TOTAL=$((REV_VECTOR_TOTAL + 1))
        REV_TOTAL=$((REV_TOTAL + 1))
        run_single_test "build/test_${funcname}_vector_reverse" "$funcname" "REV_VEC"
    fi
    
    echo ""
}

# Function to run test in a directory (nested mode)
run_test_in_dir() {
    local subdir=$1
    local dirname=$(basename "$subdir")
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    print_status "INFO" "Testing $dirname in $subdir"
    
    # Check if this function is marked as "Unable to process"
    if [ -f "$subdir/.unable_to_process" ]; then
        UNABLE_TO_PROCESS_TESTS=$((UNABLE_TO_PROCESS_TESTS + 1))
        UNABLE_TO_PROCESS_LIST+=("$dirname")
        print_status "UNABLE_TO_PROCESS" "$dirname: Insufficient parameter documentation"
        echo ""
        return
    fi
    
    # Check if Tapenade failed for this function (look in mode subdirectories)
    # Tapenade preserves the source file extension (.f or .f90), so check both
    local has_any_mode=false
    if [ "$RUN_D" = "true" ] && ([ -f "$subdir/d/${dirname}_d.f" ] || [ -f "$subdir/d/${dirname}_d.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_DV" = "true" ] && ([ -f "$subdir/dv/${dirname}_dv.f" ] || [ -f "$subdir/dv/${dirname}_dv.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_B" = "true" ] && ([ -f "$subdir/b/${dirname}_b.f" ] || [ -f "$subdir/b/${dirname}_b.f90" ]); then
        has_any_mode=true
    fi
    if [ "$RUN_BV" = "true" ] && ([ -f "$subdir/bv/${dirname}_bv.f" ] || [ -f "$subdir/bv/${dirname}_bv.f90" ]); then
        has_any_mode=true
    fi
    
    if [ "$has_any_mode" = "false" ]; then
        TAPENADE_FAILED_TESTS=$((TAPENADE_FAILED_TESTS + 1))
        TAPENADE_FAILED_LIST+=("$dirname")
        print_status "TAPENADE_FAILED" "$dirname: Tapenade fails to differentiate the code"
        echo ""
        return
    fi
    
    cd "$subdir"
    
    # Run scalar forward mode test
    if [ "$RUN_D" = "true" ]; then
        FWD_SCALAR_TOTAL=$((FWD_SCALAR_TOTAL + 1))
        FWD_TOTAL=$((FWD_TOTAL + 1))
        run_single_test "d/test_$dirname" "$dirname" "FWD"
    fi
    
    # Run vector forward mode test  
    if [ "$RUN_DV" = "true" ]; then
        FWD_VECTOR_TOTAL=$((FWD_VECTOR_TOTAL + 1))
        FWD_TOTAL=$((FWD_TOTAL + 1))
        run_single_test "dv/test_${dirname}_vector_forward" "$dirname" "FWD_VEC"
    fi
    
    # Run scalar reverse mode test
    if [ "$RUN_B" = "true" ]; then
        REV_SCALAR_TOTAL=$((REV_SCALAR_TOTAL + 1))
        REV_TOTAL=$((REV_TOTAL + 1))
        run_single_test "b/test_${dirname}_reverse" "$dirname" "REV"
    fi
    
    # Run vector reverse mode test
    if [ "$RUN_BV" = "true" ]; then
        REV_VECTOR_TOTAL=$((REV_VECTOR_TOTAL + 1))
        REV_TOTAL=$((REV_TOTAL + 1))
        run_single_test "bv/test_${dirname}_vector_reverse" "$dirname" "REV_VEC"
    fi
    
    cd - > /dev/null
    echo ""
}

# Main execution
main() {
    echo "=========================================="
    echo "Running differentiated BLAS function tests"
    echo "=========================================="
    echo "Working directory: $SCRIPT_DIR"
    echo "Mode: $([ "$FLAT_MODE" = "true" ] && echo "Flat" || echo "Nested")"
    echo ""

    if [ "$FLAT_MODE" = "true" ]; then
        # Flat mode: find functions by scanning build/ for test executables or src/ for differentiated files
        funcs=()
        
        # First try to find test executables in build/
        if [ -d "build" ]; then
            for testexe in $(ls build/test_* 2>/dev/null | grep -v '\.' | sort); do
                basename=$(basename "$testexe")
                # Extract function name: test_funcname -> funcname, test_funcname_reverse -> funcname
                funcname=$(echo "$basename" | sed -E 's/^test_//; s/_(reverse|vector_forward|vector_reverse)$//')
                if [ -n "$funcname" ] && [[ ! " ${funcs[*]} " =~ " ${funcname} " ]]; then
                    funcs+=("$funcname")
                fi
            done
        fi
        
        # If no test executables, try to find differentiated sources in src/
        if [ ${#funcs[@]} -eq 0 ] && [ -d "src" ]; then
            for srcfile in $(ls src/*_d.f src/*_d.f90 src/*_b.f src/*_b.f90 2>/dev/null | sort); do
                basename=$(basename "$srcfile")
                # Extract function name: funcname_d.f -> funcname
                funcname=$(echo "$basename" | sed -E 's/_(d|b|dv|bv)\.(f|f90)$//')
                if [ -n "$funcname" ] && [[ ! " ${funcs[*]} " =~ " ${funcname} " ]]; then
                    funcs+=("$funcname")
                fi
            done
        fi

        if [ ${#funcs[@]} -eq 0 ]; then
            print_status "WARN" "No test executables in build/ or differentiated files in src/ found."
            exit 0
        fi

        print_status "INFO" "Found ${#funcs[@]} functions to test: ${funcs[*]}"
        echo ""

        for func in "${funcs[@]}"; do
            run_test_for_func "$func"
        done
    else
        # Nested mode: find all subdirectories with Makefiles
        subdirs=()
        for dir in $(find . -maxdepth 1 -type d | sort); do
            # Skip current directory (.) and the script's directory
            if [[ "$dir" != "." ]] && [[ "$dir" != "$SCRIPT_DIR" ]] && [[ "$(basename "$dir")" != "$(basename "$SCRIPT_DIR")" ]]; then
                # Only include directories that contain a Makefile
                if [ -f "$dir/Makefile" ]; then
                    subdirs+=("$dir")
                fi
            fi
        done

        if [ ${#subdirs[@]} -eq 0 ]; then
            print_status "WARN" "No subdirectories found to test."
            exit 0
        fi

        print_status "INFO" "Found ${#subdirs[@]} subdirectories to test"
        echo ""

        for subdir in "${subdirs[@]}"; do
            run_test_in_dir "$subdir"
        done
    fi

    # Print comprehensive summary
    echo "=========================================="
    echo "COMPREHENSIVE TEST SUMMARY"
    echo "=========================================="
    echo -e "Total functions tested: ${CYAN}$TOTAL_TESTS${NC}"
    echo -e "Tapenade Failed: ${MAGENTA}$TAPENADE_FAILED_TESTS${NC}"
    echo -e "Unable to Process: ${RED}$UNABLE_TO_PROCESS_TESTS${NC}"
    echo ""
    
    # Forward Mode Summary (only if forward mode tests were run)
    if [ "$RUN_D" = "true" ] || [ "$RUN_DV" = "true" ]; then
        echo "=========================================="
        echo "FORWARD MODE RESULTS"
        echo "=========================================="
        echo -e "Total forward mode tests: ${CYAN}$FWD_TOTAL${NC}"
        echo ""
        
        # Scalar Forward Mode Results
        if [ "$RUN_D" = "true" ]; then
            echo -e "${BLUE}Scalar Forward Mode:${NC}"
            echo -e "  Total: ${CYAN}$FWD_SCALAR_TOTAL${NC}"
            echo -e "  Machine Precision: ${GREEN}$FWD_SCALAR_MACHINE_PRECISION${NC}"
            echo -e "  Acceptable: ${GREEN}$FWD_SCALAR_ACCEPTABLE${NC}"
            echo -e "  Outside Tolerance: ${YELLOW}$FWD_SCALAR_OUTSIDE_TOLERANCE${NC}"
            echo -e "  Execution Failed: ${RED}$FWD_SCALAR_EXECUTION_FAILED${NC}"
            echo -e "  Skipped: ${CYAN}$FWD_SCALAR_SKIPPED${NC}"
            echo ""
            
            if [ ${#FWD_SCALAR_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}FWD Scalar Machine Precision:${NC} ${FWD_SCALAR_MACHINE_PRECISION_LIST[*]}"
            fi
            if [ ${#FWD_SCALAR_ACCEPTABLE_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}FWD Scalar Acceptable:${NC} ${FWD_SCALAR_ACCEPTABLE_LIST[*]}"
            fi
            if [ ${#FWD_SCALAR_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
                echo -e "${YELLOW}FWD Scalar Outside Tolerance:${NC} ${FWD_SCALAR_OUTSIDE_TOLERANCE_LIST[*]}"
            fi
            if [ ${#FWD_SCALAR_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
                echo -e "${RED}FWD Scalar Execution Failed:${NC} ${FWD_SCALAR_EXECUTION_FAILED_LIST[*]}"
            fi
            if [ ${#FWD_SCALAR_SKIPPED_LIST[@]} -gt 0 ]; then
                echo -e "${CYAN}FWD Scalar Skipped:${NC} ${FWD_SCALAR_SKIPPED_LIST[*]}"
            fi
            echo ""
        fi
        
        # Vector Forward Mode Results
        if [ "$RUN_DV" = "true" ]; then
            echo -e "${BLUE}Vector Forward Mode:${NC}"
            echo -e "  Total: ${CYAN}$FWD_VECTOR_TOTAL${NC}"
            echo -e "  Machine Precision: ${GREEN}$FWD_VECTOR_MACHINE_PRECISION${NC}"
            echo -e "  Acceptable: ${GREEN}$FWD_VECTOR_ACCEPTABLE${NC}"
            echo -e "  Outside Tolerance: ${YELLOW}$FWD_VECTOR_OUTSIDE_TOLERANCE${NC}"
            echo -e "  Execution Failed: ${RED}$FWD_VECTOR_EXECUTION_FAILED${NC}"
            echo -e "  Skipped: ${CYAN}$FWD_VECTOR_SKIPPED${NC}"
            echo ""
            
            if [ ${#FWD_VECTOR_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}FWD Vector Machine Precision:${NC} ${FWD_VECTOR_MACHINE_PRECISION_LIST[*]}"
            fi
            if [ ${#FWD_VECTOR_ACCEPTABLE_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}FWD Vector Acceptable:${NC} ${FWD_VECTOR_ACCEPTABLE_LIST[*]}"
            fi
            if [ ${#FWD_VECTOR_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
                echo -e "${YELLOW}FWD Vector Outside Tolerance:${NC} ${FWD_VECTOR_OUTSIDE_TOLERANCE_LIST[*]}"
            fi
            if [ ${#FWD_VECTOR_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
                echo -e "${RED}FWD Vector Execution Failed:${NC} ${FWD_VECTOR_EXECUTION_FAILED_LIST[*]}"
            fi
            if [ ${#FWD_VECTOR_SKIPPED_LIST[@]} -gt 0 ]; then
                echo -e "${CYAN}FWD Vector Skipped:${NC} ${FWD_VECTOR_SKIPPED_LIST[*]}"
            fi
            echo ""
        fi
        
        # Combined Forward Mode Results (for backward compatibility)
        echo -e "${BLUE}Combined Forward Mode:${NC}"
        echo -e "  Machine Precision: ${GREEN}$FWD_MACHINE_PRECISION${NC}"
        echo -e "  Acceptable: ${GREEN}$FWD_ACCEPTABLE${NC}"
        echo -e "  Outside Tolerance: ${YELLOW}$FWD_OUTSIDE_TOLERANCE${NC}"
        echo -e "  Execution Failed: ${RED}$FWD_EXECUTION_FAILED${NC}"
        echo -e "  Skipped: ${CYAN}$FWD_SKIPPED${NC}"
        echo ""
        
        if [ ${#FWD_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
            echo -e "${GREEN}FWD Machine Precision:${NC} ${FWD_MACHINE_PRECISION_LIST[*]}"
        fi
        if [ ${#FWD_ACCEPTABLE_LIST[@]} -gt 0 ]; then
            echo -e "${GREEN}FWD Acceptable:${NC} ${FWD_ACCEPTABLE_LIST[*]}"
        fi
        if [ ${#FWD_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
            echo -e "${YELLOW}FWD Outside Tolerance:${NC} ${FWD_OUTSIDE_TOLERANCE_LIST[*]}"
        fi
        if [ ${#FWD_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
            echo -e "${RED}FWD Execution Failed:${NC} ${FWD_EXECUTION_FAILED_LIST[*]}"
        fi
        if [ ${#FWD_SKIPPED_LIST[@]} -gt 0 ]; then
            echo -e "${CYAN}FWD Skipped:${NC} ${FWD_SKIPPED_LIST[*]}"
        fi
    fi
    
    echo ""
    
    # Reverse Mode Summary (only if reverse mode tests were run)
    if [ "$RUN_B" = "true" ] || [ "$RUN_BV" = "true" ]; then
        echo "=========================================="
        echo "REVERSE MODE RESULTS"
        echo "=========================================="
        echo -e "Total reverse mode tests: ${CYAN}$REV_TOTAL${NC}"
        echo ""
        
        # Scalar Reverse Mode Results
        if [ "$RUN_B" = "true" ]; then
            echo -e "${BLUE}Scalar Reverse Mode:${NC}"
            echo -e "  Total: ${CYAN}$REV_SCALAR_TOTAL${NC}"
            echo -e "  Machine Precision: ${GREEN}$REV_SCALAR_MACHINE_PRECISION${NC}"
            echo -e "  Acceptable: ${GREEN}$REV_SCALAR_ACCEPTABLE${NC}"
            echo -e "  Outside Tolerance: ${YELLOW}$REV_SCALAR_OUTSIDE_TOLERANCE${NC}"
            echo -e "  Execution Failed: ${RED}$REV_SCALAR_EXECUTION_FAILED${NC}"
            echo -e "  Skipped: ${CYAN}$REV_SCALAR_SKIPPED${NC}"
            echo ""
            
            if [ ${#REV_SCALAR_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}REV Scalar Machine Precision:${NC} ${REV_SCALAR_MACHINE_PRECISION_LIST[*]}"
            fi
            if [ ${#REV_SCALAR_ACCEPTABLE_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}REV Scalar Acceptable:${NC} ${REV_SCALAR_ACCEPTABLE_LIST[*]}"
            fi
            if [ ${#REV_SCALAR_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
                echo -e "${YELLOW}REV Scalar Outside Tolerance:${NC} ${REV_SCALAR_OUTSIDE_TOLERANCE_LIST[*]}"
            fi
            if [ ${#REV_SCALAR_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
                echo -e "${RED}REV Scalar Execution Failed:${NC} ${REV_SCALAR_EXECUTION_FAILED_LIST[*]}"
            fi
            if [ ${#REV_SCALAR_SKIPPED_LIST[@]} -gt 0 ]; then
                echo -e "${CYAN}REV Scalar Skipped:${NC} ${REV_SCALAR_SKIPPED_LIST[*]}"
            fi
            echo ""
        fi
        
        # Vector Reverse Mode Results
        if [ "$RUN_BV" = "true" ]; then
            echo -e "${BLUE}Vector Reverse Mode:${NC}"
            echo -e "  Total: ${CYAN}$REV_VECTOR_TOTAL${NC}"
            echo -e "  Machine Precision: ${GREEN}$REV_VECTOR_MACHINE_PRECISION${NC}"
            echo -e "  Acceptable: ${GREEN}$REV_VECTOR_ACCEPTABLE${NC}"
            echo -e "  Outside Tolerance: ${YELLOW}$REV_VECTOR_OUTSIDE_TOLERANCE${NC}"
            echo -e "  Execution Failed: ${RED}$REV_VECTOR_EXECUTION_FAILED${NC}"
            echo -e "  Skipped: ${CYAN}$REV_VECTOR_SKIPPED${NC}"
            echo ""
            
            if [ ${#REV_VECTOR_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}REV Vector Machine Precision:${NC} ${REV_VECTOR_MACHINE_PRECISION_LIST[*]}"
            fi
            if [ ${#REV_VECTOR_ACCEPTABLE_LIST[@]} -gt 0 ]; then
                echo -e "${GREEN}REV Vector Acceptable:${NC} ${REV_VECTOR_ACCEPTABLE_LIST[*]}"
            fi
            if [ ${#REV_VECTOR_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
                echo -e "${YELLOW}REV Vector Outside Tolerance:${NC} ${REV_VECTOR_OUTSIDE_TOLERANCE_LIST[*]}"
            fi
            if [ ${#REV_VECTOR_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
                echo -e "${RED}REV Vector Execution Failed:${NC} ${REV_VECTOR_EXECUTION_FAILED_LIST[*]}"
            fi
            if [ ${#REV_VECTOR_SKIPPED_LIST[@]} -gt 0 ]; then
                echo -e "${CYAN}REV Vector Skipped:${NC} ${REV_VECTOR_SKIPPED_LIST[*]}"
            fi
            echo ""
        fi
        
        # Combined Reverse Mode Results (for backward compatibility)
        echo -e "${BLUE}Combined Reverse Mode:${NC}"
        echo -e "  Machine Precision: ${GREEN}$REV_MACHINE_PRECISION${NC}"
        echo -e "  Acceptable: ${GREEN}$REV_ACCEPTABLE${NC}"
        echo -e "  Outside Tolerance: ${YELLOW}$REV_OUTSIDE_TOLERANCE${NC}"
        echo -e "  Execution Failed: ${RED}$REV_EXECUTION_FAILED${NC}"
        echo -e "  Skipped: ${CYAN}$REV_SKIPPED${NC}"
        echo ""
        
        if [ ${#REV_MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
            echo -e "${GREEN}REV Machine Precision:${NC} ${REV_MACHINE_PRECISION_LIST[*]}"
        fi
        if [ ${#REV_ACCEPTABLE_LIST[@]} -gt 0 ]; then
            echo -e "${GREEN}REV Acceptable:${NC} ${REV_ACCEPTABLE_LIST[*]}"
        fi
        if [ ${#REV_OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
            echo -e "${YELLOW}REV Outside Tolerance:${NC} ${REV_OUTSIDE_TOLERANCE_LIST[*]}"
        fi
        if [ ${#REV_EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
            echo -e "${RED}REV Execution Failed:${NC} ${REV_EXECUTION_FAILED_LIST[*]}"
        fi
        if [ ${#REV_SKIPPED_LIST[@]} -gt 0 ]; then
            echo -e "${CYAN}REV Skipped:${NC} ${REV_SKIPPED_LIST[*]}"
        fi
    fi
    
    echo ""
    
    # Tapenade failures
    if [ ${#TAPENADE_FAILED_LIST[@]} -gt 0 ]; then
        echo -e "${MAGENTA}Tapenade Failed:${NC} ${TAPENADE_FAILED_LIST[*]}"
        echo ""
    fi
    
    # Unable to process
    if [ $UNABLE_TO_PROCESS_TESTS -gt 0 ]; then
        echo -e "${RED}Unable to Process:${NC} $UNABLE_TO_PROCESS_TESTS"
        echo ""
    fi
    if [ ${#UNABLE_TO_PROCESS_LIST[@]} -gt 0 ]; then
        echo -e "${RED}Unable to Process:${NC} ${UNABLE_TO_PROCESS_LIST[*]}"
        echo ""
    fi
    
    # Calculate overall success rate
    local fwd_success=$((FWD_MACHINE_PRECISION + FWD_ACCEPTABLE))
    local fwd_scalar_success=$((FWD_SCALAR_MACHINE_PRECISION + FWD_SCALAR_ACCEPTABLE))
    local fwd_vector_success=$((FWD_VECTOR_MACHINE_PRECISION + FWD_VECTOR_ACCEPTABLE))
    local rev_success=$((REV_MACHINE_PRECISION + REV_ACCEPTABLE))
    local rev_scalar_success=$((REV_SCALAR_MACHINE_PRECISION + REV_SCALAR_ACCEPTABLE))
    local rev_vector_success=$((REV_VECTOR_MACHINE_PRECISION + REV_VECTOR_ACCEPTABLE))
    local fwd_executed=$((FWD_TOTAL - FWD_SKIPPED))
    local fwd_scalar_executed=$((FWD_SCALAR_TOTAL - FWD_SCALAR_SKIPPED))
    local fwd_vector_executed=$((FWD_VECTOR_TOTAL - FWD_VECTOR_SKIPPED))
    local rev_executed=$((REV_TOTAL - REV_SKIPPED))
    local rev_scalar_executed=$((REV_SCALAR_TOTAL - REV_SCALAR_SKIPPED))
    local rev_vector_executed=$((REV_VECTOR_TOTAL - REV_VECTOR_SKIPPED))
    local total_failures=$((FWD_EXECUTION_FAILED + REV_EXECUTION_FAILED))
    local total_tolerance_issues=$((FWD_OUTSIDE_TOLERANCE + REV_OUTSIDE_TOLERANCE))
    
    echo "=========================================="
    echo "OVERALL RESULTS"
    echo "=========================================="
    if [ "$RUN_D" = "true" ] || [ "$RUN_DV" = "true" ]; then
        echo -e "Forward Mode (Combined): ${fwd_success}/${fwd_executed} successful"
        if [ "$RUN_D" = "true" ]; then
            echo -e "  Scalar Forward Mode: ${fwd_scalar_success}/${fwd_scalar_executed} successful"
        fi
        if [ "$RUN_DV" = "true" ]; then
            echo -e "  Vector Forward Mode: ${fwd_vector_success}/${fwd_vector_executed} successful"
        fi
    fi
    if [ "$RUN_B" = "true" ] || [ "$RUN_BV" = "true" ]; then
        echo -e "Reverse Mode (Combined): ${rev_success}/${rev_executed} successful"
        if [ "$RUN_B" = "true" ]; then
            echo -e "  Scalar Reverse Mode: ${rev_scalar_success}/${rev_scalar_executed} successful"
        fi
        if [ "$RUN_BV" = "true" ]; then
            echo -e "  Vector Reverse Mode: ${rev_vector_success}/${rev_vector_executed} successful"
        fi
    fi
    echo ""
    
    if [ $total_failures -eq 0 ] && [ $total_tolerance_issues -eq 0 ]; then
        echo -e "${GREEN}Overall result: ALL TESTS PASSED${NC}"
        return 0
    elif [ $total_failures -eq 0 ]; then
        echo -e "${YELLOW}Overall result: TESTS COMPLETED WITH SOME TOLERANCE ISSUES${NC}"
        return 0
    else
        echo -e "${RED}Overall result: SOME TESTS FAILED EXECUTION${NC}"
        return 1
    fi
}

# Handle command line arguments
case "${1:-}" in
    -h|--help)
        echo "Usage: $(basename "$0") [options]"
        echo ""
        echo "Options:"
        echo "  -h, --help     Show this help message"
        echo "  -v, --verbose  Show more detailed output"
        echo ""
        echo "This script runs tests in all subdirectories of the current directory."
        echo "Each subdirectory should contain a test executable named test_<function>."
        exit 0
        ;;
    -v|--verbose)
        set -x  # Enable debug mode
        shift
        ;;
    *)
        # No arguments or unknown arguments, run main
        ;;
esac

# Run main function and capture exit code
main "$@"
exit_code=$?

# Always exit with the captured exit code (ensures summary is printed)
exit $exit_code
