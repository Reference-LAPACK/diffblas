#!/bin/bash

# Test script for Python interfaces (original and differentiated)
# This script integrates with the existing test infrastructure

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_TEST_SCRIPT="$SCRIPT_DIR/test_python_interfaces.py"
REPORT_FILE="$SCRIPT_DIR/python_interface_test_report.json"

# Colors for output
RED='[0;31m'
GREEN='[0;32m'
YELLOW='[1;33m'
BLUE='[0;34m'
MAGENTA='[0;35m'
CYAN='[0;36m'
NC='[0m' # No Color

# Counters for reporting
TOTAL_FUNCTIONS=0
ORIGINAL_SUCCESS=0
ORIGINAL_FAILED=0
DIFFERENTIATED_SUCCESS=0
DIFFERENTIATED_FAILED=0
PYTHON_ERRORS=0

# Arrays to store results
ORIGINAL_SUCCESS_LIST=()
ORIGINAL_FAILED_LIST=()
DIFFERENTIATED_SUCCESS_LIST=()
DIFFERENTIATED_FAILED_LIST=()
PYTHON_ERROR_LIST=()

# Function to print colored status
print_status() {
    local status=$1
    local message=$2
    case $status in
        "SUCCESS")
            echo -e "${GREEN}[SUCCESS]${NC} $message"
            ;;
        "WARNING")
            echo -e "${YELLOW}[WARNING]${NC} $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} $message"
            ;;
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        "SKIP")
            echo -e "${CYAN}[SKIP]${NC} $message"
            ;;
        "FAIL")
            echo -e "${MAGENTA}[FAIL]${NC} $message"
            ;;
        *)
            echo -e "[$status] $message"
            ;;
    esac
}

# Function to check if Python is available
check_python() {
    if ! command -v python3 &> /dev/null; then
        print_status "ERROR" "Python3 is not available"
        return 1
    fi
    
    # Check if required Python packages are available
    python3 -c "import numpy, ctypes" 2>/dev/null || {
        print_status "ERROR" "Required Python packages (numpy, ctypes) are not available"
        return 1
    }
    
    return 0
}

# Function to check if Python test script exists
check_test_script() {
    if [ ! -f "$PYTHON_TEST_SCRIPT" ]; then
        print_status "ERROR" "Python test script not found: $PYTHON_TEST_SCRIPT"
        return 1
    fi
    
    if [ ! -x "$PYTHON_TEST_SCRIPT" ]; then
        print_status "WARNING" "Python test script is not executable, making it executable"
        chmod +x "$PYTHON_TEST_SCRIPT"
    fi
    
    return 0
}

# Function to run Python interface tests
run_python_tests() {
    local function_filter="$1"
    local verbose_flag="$2"
    
    print_status "INFO" "Running Python interface tests..."
    
    # Build command
    local cmd="python3 $PYTHON_TEST_SCRIPT --base-dir $SCRIPT_DIR --output $REPORT_FILE"
    
    if [ -n "$function_filter" ]; then
        cmd="$cmd --function $function_filter"
    fi
    
    if [ "$verbose_flag" = "true" ]; then
        cmd="$cmd --verbose"
    fi
    
    # Run the tests
    print_status "INFO" "Executing: $cmd"
    
    if eval "$cmd"; then
        print_status "SUCCESS" "Python interface tests completed successfully"
        return 0
    else
        local exit_code=$?
        print_status "ERROR" "Python interface tests failed with exit code $exit_code"
        return $exit_code
    fi
}

# Function to parse JSON report and extract statistics
parse_report() {
    if [ ! -f "$REPORT_FILE" ]; then
        print_status "ERROR" "Report file not found: $REPORT_FILE"
        return 1
    fi
    
    # Extract statistics using Python (more reliable than shell parsing)
    local stats=$(python3 -c "
import json
import sys
try:
    with open('$REPORT_FILE', 'r') as f:
        report = json.load(f)
    summary = report['summary']
    print(f\"{summary['total_functions']},{summary['original_success']},{summary['original_failed']},{summary['differentiated_success']},{summary['differentiated_failed']}\")
except Exception as e:
    print(f\"0,0,0,0,0\")
    sys.exit(1)
")
    
    if [ $? -ne 0 ]; then
        print_status "ERROR" "Failed to parse report file"
        return 1
    fi
    
    # Parse the statistics
    IFS=',' read -r TOTAL_FUNCTIONS ORIGINAL_SUCCESS ORIGINAL_FAILED DIFFERENTIATED_SUCCESS DIFFERENTIATED_FAILED <<< "$stats"
    
    return 0
}

# Function to print detailed results
print_detailed_results() {
    if [ ! -f "$REPORT_FILE" ]; then
        print_status "WARNING" "No detailed report available"
        return
    fi
    
    print_status "INFO" "Detailed results from $REPORT_FILE:"
    
    # Use Python to format the results nicely
    python3 -c "
import json
import sys

try:
    with open('$REPORT_FILE', 'r') as f:
        report = json.load(f)
    
    print('\n' + '='*60)
    print('DETAILED PYTHON INTERFACE TEST RESULTS')
    print('='*60)
    
    for func_name, results in report['results'].items():
        print(f'\n{func_name.upper()}:')
        
        # Original interface
        orig = results['original']
        if orig['status'] == 'SUCCESS':
            print(f'  Original: SUCCESS ({orig[\"execution_time\"]:.3f}s)')
        else:
            print(f'  Original: {orig[\"status\"]} - {orig.get(\"error\", \"Unknown error\")}')
        
        # Differentiated interface
        diff = results['differentiated']
        if diff['status'] == 'SUCCESS':
            print(f'  Differentiated: SUCCESS ({diff[\"execution_time\"]:.3f}s)')
        else:
            print(f'  Differentiated: {diff[\"status\"]} - {diff.get(\"error\", \"Unknown error\")}')
    
    print('\n' + '='*60)
    
except Exception as e:
    print(f'Error reading report: {e}')
    sys.exit(1)
"
}

# Function to print summary
print_summary() {
    echo ""
    echo "=========================================="
    echo "PYTHON INTERFACE TEST SUMMARY"
    echo "=========================================="
    echo -e "Total functions tested: ${CYAN}$TOTAL_FUNCTIONS${NC}"
    echo ""
    echo -e "Original interfaces: ${GREEN}$ORIGINAL_SUCCESS${NC} success, ${RED}$ORIGINAL_FAILED${NC} failed"
    echo -e "Differentiated interfaces: ${GREEN}$DIFFERENTIATED_SUCCESS${NC} success, ${RED}$DIFFERENTIATED_FAILED${NC} failed"
    echo ""
    
    local total_tests=$((ORIGINAL_SUCCESS + ORIGINAL_FAILED + DIFFERENTIATED_SUCCESS + DIFFERENTIATED_FAILED))
    local total_success=$((ORIGINAL_SUCCESS + DIFFERENTIATED_SUCCESS))
    local total_failed=$((ORIGINAL_FAILED + DIFFERENTIATED_FAILED))
    
    echo -e "Overall: ${GREEN}$total_success${NC} success, ${RED}$total_failed${NC} failed out of $total_tests tests"
    
    if [ $total_failed -eq 0 ]; then
        echo -e "${GREEN}Overall result: ALL PYTHON INTERFACE TESTS PASSED${NC}"
        return 0
    else
        echo -e "${YELLOW}Overall result: SOME PYTHON INTERFACE TESTS FAILED${NC}"
        return 1
    fi
}

# Main execution function
main() {
    local function_filter=""
    local verbose=false
    local detailed=false
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                echo "Usage: $(basename "$0") [options]"
                echo ""
                echo "Options:"
                echo "  -h, --help          Show this help message"
                echo "  -v, --verbose       Show verbose output"
                echo "  -d, --detailed      Show detailed results"
                echo "  -f, --function FUNC Test specific function only"
                echo ""
                echo "This script tests Python interfaces for BLAS functions."
                echo "It runs both original and differentiated interface tests."
                exit 0
                ;;
            -v|--verbose)
                verbose=true
                shift
                ;;
            -d|--detailed)
                detailed=true
                shift
                ;;
            -f|--function)
                function_filter="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done
    
    echo "=========================================="
    echo "Python Interface Test Suite"
    echo "=========================================="
    echo "Working directory: $SCRIPT_DIR"
    echo ""
    
    # Check prerequisites
    if ! check_python; then
        exit 1
    fi
    
    if ! check_test_script; then
        exit 1
    fi
    
    # Run the tests
    if ! run_python_tests "$function_filter" "$verbose"; then
        print_status "ERROR" "Python interface tests failed"
        exit 1
    fi
    
    # Parse the report
    if ! parse_report; then
        print_status "ERROR" "Failed to parse test results"
        exit 1
    fi
    
    # Print detailed results if requested
    if [ "$detailed" = "true" ]; then
        print_detailed_results
    fi
    
    # Print summary
    print_summary
}

# Run main function with all arguments
main "$@"
