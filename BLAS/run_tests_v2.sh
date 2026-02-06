#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TEST_TIMEOUT=10

# ======================================
# Parse mode option
# ======================================

case "$1" in
    -d)
        MODE="forward"
        SUFFIX=""
        ;;
    -dv)
        MODE="vector-forward"
        SUFFIX="_vector_forward"
        ;;
    -b)
        MODE="reverse"
        SUFFIX="_reverse"
        ;;
    -bv)
        MODE="vector-reverse"
        SUFFIX="_vector_reverse"
        ;;
    *)
        echo "Usage: $0 {-d|-dv|-b|-bv}"
        exit 1
        ;;
esac

LOGFILE="results-${MODE}.log"
> "$LOGFILE"

# ---------- Colors ----------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BLUE='\033[0;34m'
NC='\033[0m'

# ---------- Print ----------
print_status() {
    local s=$1
    local msg=$2

    case $s in
        MACHINE_PRECISION) echo -e "${GREEN}[MACHINE_PRECISION]${NC} $msg" ;;
        ACCEPTABLE)        echo -e "${GREEN}[ACCEPTABLE]${NC} $msg" ;;
        OUTSIDE_TOLERANCE) echo -e "${YELLOW}[OUTSIDE_TOLERANCE]${NC} $msg" ;;
        EXECUTION_FAILED)  echo -e "${RED}[EXECUTION_FAILED]${NC} $msg" ;;
        SKIPPED)           echo -e "${CYAN}[SKIPPED]${NC} $msg" ;;
        INFO)              echo -e "${BLUE}[INFO]${NC} $msg" ;;
        *) echo "[$s] $msg" ;;
    esac
}

# ---------- Julia log ----------
write_log() {
    echo "$1:$MODE:$2" >> "$LOGFILE"
}

# ---------- Run with timeout ----------
safe_run() {
    local exe=$1
    local out=$2

    timeout ${TEST_TIMEOUT}s ./"$exe" > "$out" 2>&1
    code=$?

    if [ $code -eq 124 ]; then
        echo "Test timed out" >> "$out"
        return 124
    fi

    if [ ! -s "$out" ]; then
        echo "Test crashed" >> "$out"
        return 1
    fi

    return $code
}

# ---------- Detect result ----------
detect_status() {
    local file=$1

    if grep -qiE "segmentation fault|aborted|floating point exception|timed out|illegal value|cannot open shared object" "$file"; then
        echo "EXECUTION_FAILED"
        return
    fi

    if grep -qi "FAIL: Large errors detected" "$file"; then
        echo "OUTSIDE_TOLERANCE"
        return
    fi

    if grep -qi "WARNING: .*significant errors" "$file"; then
        echo "OUTSIDE_TOLERANCE"
        return
    fi

    if grep -qi "PASS: .*machine precision" "$file"; then
        echo "MACHINE_PRECISION"
        return
    fi

    if grep -qi "PASS: .*reasonably accurate" "$file"; then
        echo "ACCEPTABLE"
        return
    fi

    echo "ACCEPTABLE"
}

# ==========================================
# Main
# ==========================================

echo "================================="
echo "Running mode: $MODE"
echo "Logfile: $LOGFILE"
echo "================================="
echo ""

if [ ! -d build ]; then
    echo "No build/ directory found"
    exit 1
fi

# Collect routines
routines=()

for exe in build/test_*; do
    [ -e "$exe" ] || continue

    base=$(basename "$exe")

    name=$(echo "$base" | sed -E "s/^test_//; s/${SUFFIX}$//")

    if [[ ! " ${routines[*]} " =~ " ${name} " ]]; then
        routines+=("$name")
    fi
done

if [ ${#routines[@]} -eq 0 ]; then
    echo "No tests found."
    exit 0
fi

print_status "INFO" "Found ${#routines[@]} routines: ${routines[*]}"
echo ""

# Run only selected mode

for r in "${routines[@]}"; do

    exe="build/test_${r}${SUFFIX}"

    if [ "$SUFFIX" = "" ]; then
        exe="build/test_$r"
    fi

    if [ ! -f "$exe" ] || [ ! -x "$exe" ]; then
        print_status "SKIPPED" "$r"
        write_log "$r" "SKIPPED"
        continue
    fi

    out="output_${r}_${MODE}.log"

    safe_run "$exe" "$out"
    status=$(detect_status "$out")

    print_status "$status" "$r"
    write_log "$r" "$status"
done

echo ""
echo "================================="
echo "Log written to: $LOGFILE"
echo "================================="
