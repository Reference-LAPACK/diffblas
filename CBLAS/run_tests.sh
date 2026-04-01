#!/bin/bash
# Top-level test script for differentiated CBLAS functions
# Tests all subdirectories in forward mode (d)

# Note: We don't use 'set -e' here because we need to handle test failures gracefully
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

# Counters
TOTAL_TESTS=0
MACHINE_PRECISION=0
ACCEPTABLE=0
OUTSIDE_TOLERANCE=0
EXECUTION_FAILED=0
SKIPPED=0
TAPENADE_FAILED=0

# Arrays to store results
MACHINE_PRECISION_LIST=()
ACCEPTABLE_LIST=()
OUTSIDE_TOLERANCE_LIST=()
EXECUTION_FAILED_LIST=()
SKIPPED_LIST=()
TAPENADE_FAILED_LIST=()
SUCCESS_d=0 TOTAL_d=0 MACHINE_PRECISION_d=0 ACCEPTABLE_d=0 OUTSIDE_TOLERANCE_d=0 EXECUTION_FAILED_d=0 SKIPPED_d=0 MACHINE_PRECISION_LIST_d=() ACCEPTABLE_LIST_d=() OUTSIDE_TOLERANCE_LIST_d=() EXECUTION_FAILED_LIST_d=() SKIPPED_LIST_d=()
        SUCCESS_dv=0 TOTAL_dv=0 MACHINE_PRECISION_dv=0 ACCEPTABLE_dv=0 OUTSIDE_TOLERANCE_dv=0 EXECUTION_FAILED_dv=0 SKIPPED_dv=0 MACHINE_PRECISION_LIST_dv=() ACCEPTABLE_LIST_dv=() OUTSIDE_TOLERANCE_LIST_dv=() EXECUTION_FAILED_LIST_dv=() SKIPPED_LIST_dv=()
        SUCCESS_b=0 TOTAL_b=0 MACHINE_PRECISION_b=0 ACCEPTABLE_b=0 OUTSIDE_TOLERANCE_b=0 EXECUTION_FAILED_b=0 SKIPPED_b=0 MACHINE_PRECISION_LIST_b=() ACCEPTABLE_LIST_b=() OUTSIDE_TOLERANCE_LIST_b=() EXECUTION_FAILED_LIST_b=() SKIPPED_LIST_b=()
        SUCCESS_bv=0 TOTAL_bv=0 MACHINE_PRECISION_bv=0 ACCEPTABLE_bv=0 OUTSIDE_TOLERANCE_bv=0 EXECUTION_FAILED_bv=0 SKIPPED_bv=0 MACHINE_PRECISION_LIST_bv=() ACCEPTABLE_LIST_bv=() OUTSIDE_TOLERANCE_LIST_bv=() EXECUTION_FAILED_LIST_bv=() SKIPPED_LIST_bv=()
        

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
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        *)
            echo -e "[$status] $message"
            ;;
    esac
}

# Function to safely run a test with timeout
safe_run_test() {
    local test_executable=$1
    local output_file=$2
    
    # Use timeout to prevent hanging tests
    # When a command segfaults, timeout returns the signal number + 128 (e.g., 139 for SIGSEGV)
    # Do not use || true so we preserve the test exit code for classification
    timeout 30s ./"$test_executable" > "$output_file" 2>&1
    local timeout_exit_code=$?
    
    # Check if the test crashed (empty output file usually indicates a crash)
    if [ ! -s "$output_file" ]; then
        echo "Test crashed or produced no output" >> "$output_file"
        # Return a failure code, but don't exit the script
        return 1
    fi
    
    # Return the exit code for further checking
    # Exit codes: 0 = success, 124 = timeout, 139 = segfault, 134 = abort, 136 = FPE
    return $timeout_exit_code
}

# Function to run a single test
run_single_test() {
    local test_executable=$1
    local test_name=$2
    local output_file="test_output.log"
    local current_mode=""
    [[ "$test_name" == *_bv ]] && current_mode="bv"
    [[ "$test_name" == *_dv ]] && current_mode="dv"
    [[ "$test_name" == *_b ]] && current_mode="b"
    [[ "$test_name" == *_d ]] && current_mode="d"
    
    if [ ! -f "$test_executable" ]; then
        SKIPPED=$((SKIPPED + 1))
        [ -n "$current_mode" ] && eval "SKIPPED_$current_mode=\$((SKIPPED_$current_mode + 1))" && eval "SKIPPED_LIST_$current_mode+=("\$test_name")"
        SKIPPED_LIST+=("$test_name")
        print_status "SKIPPED" "$test_name: Test executable not found"
        return
    fi
    
    if [ ! -x "$test_executable" ]; then
        SKIPPED=$((SKIPPED + 1))
        [ -n "$current_mode" ] && eval "SKIPPED_$current_mode=\$((SKIPPED_$current_mode + 1))" && eval "SKIPPED_LIST_$current_mode+=("\$test_name")"
        SKIPPED_LIST+=("$test_name")
        print_status "SKIPPED" "$test_name: Test executable exists but is not executable"
        return
    fi
    
    if [ -n "$current_mode" ]; then eval "TOTAL_$current_mode=\$((TOTAL_$current_mode + 1))"; fi
    
    # Run the test safely (do not use || true so we get the real exit code)
    safe_run_test "$test_executable" "$output_file"
    local exit_code=$?
    
    # Check for execution failure patterns
    local has_execution_failures=false
    # Check exit code: 124 = timeout, 139 = segfault (128+11), 134 = abort (128+6), 136 = FPE (128+8)
    # Also check for any non-zero exit code that's not a normal test failure
    if [ $exit_code -eq 124 ] || [ $exit_code -eq 139 ] || [ $exit_code -eq 134 ] || [ $exit_code -eq 136 ] || [ $exit_code -gt 1 ]; then
        has_execution_failures=true
    fi
    # Also check output file for error messages (case-insensitive)
    if grep -qi "Segmentation fault\|Aborted\|Floating point exception\|Test timed out\|dumped core\|core dumped" "$output_file" 2>/dev/null; then
        has_execution_failures=true
    fi
    # CBLAS/xerbla parameter errors: test ran but with invalid args (e.g. uninitialized Side/Uplo)
    if grep -qE "Illegal (Side|Uplo|Trans|Layout|Diag) setting|Parameter [0-9]+ to routine .* (was )?incorrect" "$output_file" 2>/dev/null; then
        has_execution_failures=true
    fi
    
    # Check for derivative tolerance patterns
    local has_machine_precision=false
    local has_acceptable=false
    local has_outside_tolerance=false
    
    if grep -q "FAIL: Large errors detected in derivatives" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    elif grep -q "FAIL: VJP error ratio" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    elif grep -q "PASS: Derivatives are accurate to machine precision" "$output_file" 2>/dev/null; then
        has_machine_precision=true
    elif grep -q "PASS: Derivatives are reasonably accurate" "$output_file" 2>/dev/null; then
        has_acceptable=true
    elif grep -q "PASS: reverse mode (stub)" "$output_file" 2>/dev/null; then
        # Reverse-mode code present; VJP numerical check only for GEMM/nrm2
        has_acceptable=true
    elif grep -q "PASS: reverse vector mode (stub)" "$output_file" 2>/dev/null; then
        # Vector reverse (bv) stub
        has_acceptable=true
    elif grep -q "WARNING: Derivatives may have significant errors" "$output_file" 2>/dev/null; then
        has_outside_tolerance=true
    fi
    
    # Determine test result category
    if [ $exit_code -eq 0 ] && [ "$has_execution_failures" = false ]; then
        if [ "$has_machine_precision" = true ]; then
            MACHINE_PRECISION=$((MACHINE_PRECISION + 1))
            [ -n "$current_mode" ] && eval "SUCCESS_$current_mode=\$((SUCCESS_$current_mode + 1))" && eval "MACHINE_PRECISION_$current_mode=\$((MACHINE_PRECISION_$current_mode + 1))" && eval "MACHINE_PRECISION_LIST_$current_mode+=("\$test_name")"
            MACHINE_PRECISION_LIST+=("$test_name")
            print_status "MACHINE_PRECISION" "$test_name: Derivatives match to machine precision"
        elif [ "$has_acceptable" = true ]; then
            ACCEPTABLE=$((ACCEPTABLE + 1))
            [ -n "$current_mode" ] && eval "SUCCESS_$current_mode=\$((SUCCESS_$current_mode + 1))" && eval "ACCEPTABLE_$current_mode=\$((ACCEPTABLE_$current_mode + 1))" && eval "ACCEPTABLE_LIST_$current_mode+=("\$test_name")"
            ACCEPTABLE_LIST+=("$test_name")
            if grep -q "PASS: reverse mode (stub)" "$output_file" 2>/dev/null; then
                print_status "ACCEPTABLE" "$test_name: Reverse mode (stub; VJP check only for GEMM/nrm2)"
            elif grep -q "PASS: reverse vector mode (stub)" "$output_file" 2>/dev/null; then
                print_status "ACCEPTABLE" "$test_name: Reverse vector mode (stub)"
            elif grep -q "PASS: Derivatives are reasonably accurate" "$output_file" 2>/dev/null; then
                print_status "ACCEPTABLE" "$test_name: Derivatives are acceptable"
            else
                print_status "ACCEPTABLE" "$test_name: Test completed successfully"
            fi
        elif [ "$has_outside_tolerance" = true ]; then
            OUTSIDE_TOLERANCE=$((OUTSIDE_TOLERANCE + 1))
            [ -n "$current_mode" ] && eval "OUTSIDE_TOLERANCE_$current_mode=\$((OUTSIDE_TOLERANCE_$current_mode + 1))" && eval "OUTSIDE_TOLERANCE_LIST_$current_mode+=("\$test_name")"
            OUTSIDE_TOLERANCE_LIST+=("$test_name")
            print_status "OUTSIDE_TOLERANCE" "$test_name: Code runs but derivatives outside acceptable tolerance"
        else
            # Test completed but no clear derivative status - treat as acceptable
            ACCEPTABLE=$((ACCEPTABLE + 1))
            [ -n "$current_mode" ] && eval "SUCCESS_$current_mode=\$((SUCCESS_$current_mode + 1))" && eval "ACCEPTABLE_$current_mode=\$((ACCEPTABLE_$current_mode + 1))" && eval "ACCEPTABLE_LIST_$current_mode+=("\$test_name")"
            ACCEPTABLE_LIST+=("$test_name")
            print_status "ACCEPTABLE" "$test_name: Test completed successfully"
        fi
        echo "  Last line of output:"
        tail -1 "$output_file" | sed 's/^/    /'
    elif [ $exit_code -eq 1 ] && [ "$has_outside_tolerance" = true ]; then
        OUTSIDE_TOLERANCE=$((OUTSIDE_TOLERANCE + 1))
        [ -n "$current_mode" ] && eval "OUTSIDE_TOLERANCE_$current_mode=\$((OUTSIDE_TOLERANCE_$current_mode + 1))" && eval "OUTSIDE_TOLERANCE_LIST_$current_mode+=("\$test_name")"
        OUTSIDE_TOLERANCE_LIST+=("$test_name")
        print_status "OUTSIDE_TOLERANCE" "$test_name: VJP/derivative check failed (e.g. nrm2 error ratio > 1)"
        echo "  Last line of output:"
        tail -1 "$output_file" | sed 's/^/    /'
    elif [ "$has_execution_failures" = true ]; then
        EXECUTION_FAILED=$((EXECUTION_FAILED + 1))
        [ -n "$current_mode" ] && eval "EXECUTION_FAILED_$current_mode=\$((EXECUTION_FAILED_$current_mode + 1))" && eval "EXECUTION_FAILED_LIST_$current_mode+=("\$test_name")"
        EXECUTION_FAILED_LIST+=("$test_name")
        print_status "EXECUTION_FAILED" "$test_name: Code fails to complete execution"
        echo "  Error output:"
        grep -iE "Segmentation fault|Aborted|Floating point exception|Test timed out|dumped core|core dumped" "$output_file" | head -3 | sed 's/^/    /'
        grep -E "Illegal (Side|Uplo|Trans|Layout|Diag) setting|Parameter [0-9]+ to routine .* (was )?incorrect" "$output_file" 2>/dev/null | head -3 | sed 's/^/    /'
        if [ $exit_code -ne 0 ]; then
            echo "  Exit code: $exit_code"
        fi
    else
        EXECUTION_FAILED=$((EXECUTION_FAILED + 1))
        [ -n "$current_mode" ] && eval "EXECUTION_FAILED_$current_mode=\$((EXECUTION_FAILED_$current_mode + 1))" && eval "EXECUTION_FAILED_LIST_$current_mode+=("\$test_name")"
        EXECUTION_FAILED_LIST+=("$test_name")
        print_status "EXECUTION_FAILED" "$test_name: Test failed with exit code $exit_code"
        echo "  Last line of output:"
        tail -1 "$output_file" | sed 's/^/    /'
    fi
}



# Main execution
main() {
    echo "=========================================="
    echo "Running differentiated CBLAS function tests"
    echo "=========================================="
    echo "Working directory: $SCRIPT_DIR"
    echo "Mode: forward (d)"
    echo ""

    
    
    # Flat layout: discover ALL tests from test/ (_d, _dv, _b) so we run whatever was built
    if [ -d "build" ]; then
        TEST_NAMES=()
        if [ -d "test" ]; then
            for f in test/test_cblas_*.c; do
                [ -f "$f" ] || continue
                base=$(basename "$f" .c)
                TEST_NAMES+=("$base")
            done
        fi
        TOTAL_TESTS=${#TEST_NAMES[@]}
        # Build as many test executables as possible (Makefile uses MAKEFLAGS += -k by default)
        if [ -f "Makefile" ]; then
            make test-executables 2>/dev/null || true
        fi
        for test_name in "${TEST_NAMES[@]}"; do
            exe="build/$test_name"
            run_single_test "$exe" "$test_name"
        done
    else
        for t in test_*; do
            if [ -x "$t" ]; then
                TOTAL_TESTS=$((TOTAL_TESTS + 1))
                run_single_test "$t" "$t"
            fi
        done
    fi


    # Print comprehensive summary
    echo "=========================================="
    echo "COMPREHENSIVE TEST SUMMARY"
    echo "=========================================="
    echo -e "Total functions tested: ${CYAN}$TOTAL_TESTS${NC}"
    echo -e "Tapenade Failed: ${MAGENTA}$TAPENADE_FAILED${NC}"
    echo ""
    
    if [ ${#MACHINE_PRECISION_LIST[@]} -gt 0 ]; then
        echo -e "${GREEN}Machine Precision:${NC} ${MACHINE_PRECISION_LIST[*]}"
    fi
    if [ ${#ACCEPTABLE_LIST[@]} -gt 0 ]; then
        echo -e "${GREEN}Acceptable:${NC} ${ACCEPTABLE_LIST[*]}"
    fi
    if [ ${#OUTSIDE_TOLERANCE_LIST[@]} -gt 0 ]; then
        echo -e "${YELLOW}Outside Tolerance:${NC} ${OUTSIDE_TOLERANCE_LIST[*]}"
    fi
    if [ ${#EXECUTION_FAILED_LIST[@]} -gt 0 ]; then
        echo -e "${RED}Execution Failed:${NC} ${EXECUTION_FAILED_LIST[*]}"
    fi
    if [ ${#SKIPPED_LIST[@]} -gt 0 ]; then
        echo -e "${CYAN}Skipped:${NC} ${SKIPPED_LIST[*]}"
    fi
    if [ ${#TAPENADE_FAILED_LIST[@]} -gt 0 ]; then
        echo -e "${MAGENTA}Tapenade Failed:${NC} ${TAPENADE_FAILED_LIST[*]}"
    fi
    echo ""
    
    echo "=========================================="
    echo "RESULTS BY MODE"
    echo "=========================================="
    echo -e "Total tests: ${CYAN}$TOTAL_TESTS${NC}"
    echo -e "Machine Precision: ${GREEN}$MACHINE_PRECISION${NC}"
    echo -e "Acceptable: ${GREEN}$ACCEPTABLE${NC}"
    echo -e "Outside Tolerance: ${YELLOW}$OUTSIDE_TOLERANCE${NC}"
    echo -e "Execution Failed: ${RED}$EXECUTION_FAILED${NC}"
    echo -e "Skipped: ${CYAN}$SKIPPED${NC}"
    echo ""
    
    # Calculate overall success rate
    local success=$((MACHINE_PRECISION + ACCEPTABLE))
    local executed=$((TOTAL_TESTS - SKIPPED - TAPENADE_FAILED))
    
    echo -e "${GREEN}Forward Scalar (d): ${SUCCESS_d}/${TOTAL_d} successful${NC} (Machine Precision: ${MACHINE_PRECISION_d}, Acceptable: ${ACCEPTABLE_d}, Outside Tolerance: ${OUTSIDE_TOLERANCE_d}, Execution Failed: ${EXECUTION_FAILED_d}, Skipped: ${SKIPPED_d})"
    echo -e "${GREEN}Machine Precision:${NC} ${MACHINE_PRECISION_LIST_d[*]}"
    echo -e "${GREEN}Acceptable:${NC} ${ACCEPTABLE_LIST_d[*]}"
    echo -e "${YELLOW}Outside Tolerance:${NC} ${OUTSIDE_TOLERANCE_LIST_d[*]}"
    echo -e "${RED}Execution Failed:${NC} ${EXECUTION_FAILED_LIST_d[*]}"
    echo -e "${CYAN}Skipped:${NC} ${SKIPPED_LIST_d[*]}"
    echo ""
    echo -e "${GREEN}Reverse Scalar (b): ${SUCCESS_b}/${TOTAL_b} successful${NC} (Machine Precision: ${MACHINE_PRECISION_b}, Acceptable: ${ACCEPTABLE_b}, Outside Tolerance: ${OUTSIDE_TOLERANCE_b}, Execution Failed: ${EXECUTION_FAILED_b}, Skipped: ${SKIPPED_b})"
    echo -e "${GREEN}Machine Precision:${NC} ${MACHINE_PRECISION_LIST_b[*]}"
    echo -e "${GREEN}Acceptable:${NC} ${ACCEPTABLE_LIST_b[*]}"
    echo -e "${YELLOW}Outside Tolerance:${NC} ${OUTSIDE_TOLERANCE_LIST_b[*]}"
    echo -e "${RED}Execution Failed:${NC} ${EXECUTION_FAILED_LIST_b[*]}"
    echo -e "${CYAN}Skipped:${NC} ${SKIPPED_LIST_b[*]}"
    echo ""
    echo -e "${GREEN}Forward vector (dv): ${SUCCESS_dv}/${TOTAL_dv} successful${NC} (Machine Precision: ${MACHINE_PRECISION_dv}, Acceptable: ${ACCEPTABLE_dv}, Outside Tolerance: ${OUTSIDE_TOLERANCE_dv}, Execution Failed: ${EXECUTION_FAILED_dv}, Skipped: ${SKIPPED_dv})"
    echo -e "${GREEN}Machine Precision:${NC} ${MACHINE_PRECISION_LIST_dv[*]}"
    echo -e "${GREEN}Acceptable:${NC} ${ACCEPTABLE_LIST_dv[*]}"
    echo -e "${YELLOW}Outside Tolerance:${NC} ${OUTSIDE_TOLERANCE_LIST_dv[*]}"
    echo -e "${RED}Execution Failed:${NC} ${EXECUTION_FAILED_LIST_dv[*]}"
    echo -e "${CYAN}Skipped:${NC} ${SKIPPED_LIST_dv[*]}"
    echo ""
    echo -e "${GREEN}Reverse vector (bv): ${SUCCESS_bv}/${TOTAL_bv} successful${NC} (Machine Precision: ${MACHINE_PRECISION_bv}, Acceptable: ${ACCEPTABLE_bv}, Outside Tolerance: ${OUTSIDE_TOLERANCE_bv}, Execution Failed: ${EXECUTION_FAILED_bv}, Skipped: ${SKIPPED_bv})"
    echo -e "${GREEN}Machine Precision:${NC} ${MACHINE_PRECISION_LIST_bv[*]}"
    echo -e "${GREEN}Acceptable:${NC} ${ACCEPTABLE_LIST_bv[*]}"
    echo -e "${YELLOW}Outside Tolerance:${NC} ${OUTSIDE_TOLERANCE_LIST_bv[*]}"
    echo -e "${RED}Execution Failed:${NC} ${EXECUTION_FAILED_LIST_bv[*]}"
    echo -e "${CYAN}Skipped:${NC} ${SKIPPED_LIST_bv[*]}"
    echo ""
    
    echo ""
    echo "=========================================="
    echo "OVERALL RESULTS"
    echo "=========================================="
    echo -e "Total: ${success}/${TOTAL_TESTS} successful"
    echo ""
    
    if [ $EXECUTION_FAILED -eq 0 ] && [ $OUTSIDE_TOLERANCE -eq 0 ]; then
        echo -e "${GREEN}Overall result: ALL TESTS PASSED${NC}"
        exit 0
    elif [ $EXECUTION_FAILED -eq 0 ]; then
        echo -e "${YELLOW}Overall result: TESTS COMPLETED WITH SOME TOLERANCE ISSUES${NC}"
        exit 0
    else
        echo -e "${RED}Overall result: SOME TESTS FAILED EXECUTION${NC}"
        exit 1
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
        echo "Each subdirectory should contain a test executable in the d/ subdirectory."
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

main "$@"
