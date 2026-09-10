#!/bin/bash
set -eu

# =====================================================================
# FUNCTION: apply_cice_overrides
# 
# Description: 
#   A bash utility to apply CICE set_nml options directly 
#   to an ice_in namelist. Modifies the file in-place, including 
#   Fortran is not case-sensitive.
#
# Example usage:
#   ./CICE_override.sh ice_in < set_nml.bgc
#   echo "kice = 2" | ./CICE_override.sh ice_in
#   printf "kice = 2\ndt = 1800.0\n" | ./CICE_override.sh ice_in
# =====================================================================

apply_cice_overrides() {
    local ice_in="$1"

    if [[ -z "$ice_in" || ! -f "$ice_in" ]]; then
        echo "Error: Missing or invalid ice_in file." >&2
        return 1
    fi

    # Create a safe temp file to store edited ice_in
    local tmp_file
    tmp_file=$(mktemp)

    # Read overrides from standard input (stdin)
    while read -r line; do
        # Skip empty lines and lines starting with '#'
        [[ -z "$line" || "$line" == "#"* ]] && continue
        
        # Extract the override key and value
        local key
        key=$(echo "$line" | cut -d'=' -f1 | awk '{print $1}')
        
        local value
        value=$(echo "$line" | cut -d'=' -f2- | sed 's/^[[:space:]]*//')
        
        [[ -z "$key" || -z "$value" ]] && continue
        
        # Find how key is capitalized in ice_in using grep -i
        local exact_match
        exact_match=$(grep -i -m 1 "^[[:space:]]*${key}[[:space:]]*=" "$ice_in" || true)

        # Error if override key isn't in ice_in
        if [[ -z "$exact_match" ]]; then
            echo "ERROR: Override key '${key}' not found in '${ice_in}' " >&2
            rm -f "$tmp_file"
            return 1
        fi
        
        # Extract the exact capitalization (e.g., "kice" becomes "KICE")
        local exact_key
        exact_key=$(echo "$exact_match" | cut -d'=' -f1 | awk '{print $1}')
            
        # Overwrite key,value using standard, case-sensitive sed
        sed -E "s/^[[:space:]]*${exact_key}[[:space:]]*=.*/  ${exact_key} = ${value}/" "$ice_in" > "$tmp_file"
            
        # Overwrite ice_in
        cat "$tmp_file" > "$ice_in"

    done

    # Clean up
    rm -f "$tmp_file"
}

# =====================================================================
# UNIT TESTS
# =====================================================================
run_tests() {
    echo "Running unit tests for apply_cice_overrides..."
    local fails=0
    local pass=0

    assert_match() {
        local expected="$1"
        local file="$2"
        local test_name="$3"
        
        if grep -qF "$expected" "$file"; then
            echo "[PASS] $test_name"
            pass=$((pass + 1))
        else
            echo "[FAIL] $test_name (Expected to find: '$expected')"
            fails=$((fails + 1))
        fi
    }

    local mock_ice_in
    mock_ice_in=$(mktemp)
    
    reset_mock_ice_in() {
        # Grouped commands to avoid SC2129 multiple redirects
        {
            echo "&setup_nml"
            echo "  days_per_year = 365"
            echo "  use_restart   = .true."
            echo "  KICE = 1"
            echo "  DT = 3600.0"
            echo "/"
        } > "$mock_ice_in"
    }

    # ---------------------------------------------------------
    # TEST 1: File Redirect
    # ---------------------------------------------------------
    reset_mock_ice_in
    
    local mock_set_nml
    mock_set_nml=$(mktemp)
    
    {
        echo "days_per_year = 360"
        echo "use_restart = .false."
    } > "$mock_set_nml"

    apply_cice_overrides "$mock_ice_in" < "$mock_set_nml"
    
    assert_match "  days_per_year = 360" "$mock_ice_in" "Test 1: File Redirect (Integer)"
    assert_match "  use_restart = .false." "$mock_ice_in" "Test 1: File Redirect (Logical)"
    rm -f "$mock_set_nml"

    # ---------------------------------------------------------
    # TEST 2: Piped Multi-line (Verifying Case-Insensitivity)
    # ---------------------------------------------------------
    reset_mock_ice_in
    
    printf "dt = 1800.0\nkice = 2\n" | apply_cice_overrides "$mock_ice_in"

    assert_match "  DT = 1800.0" "$mock_ice_in" "Test 2: Piped Multi-line (Float)"
    assert_match "  KICE = 2" "$mock_ice_in" "Test 2: Case-insensitive extraction (KICE)"

    # ---------------------------------------------------------
    # TEST 3: Piped Single Line
    # ---------------------------------------------------------
    reset_mock_ice_in
    
    echo "days_per_year = 366" | apply_cice_overrides "$mock_ice_in"
    
    assert_match "  days_per_year = 366" "$mock_ice_in" "Test 3: Piped Single Line"

    # ---------------------------------------------------------
    # TEST 4: Fail on Missing/Misspelled Key
    # ---------------------------------------------------------
    reset_mock_ice_in

    # We expect this to fail, so we capture the return code
    if echo "bad_spelling = 99" | apply_cice_overrides "$mock_ice_in" 2>/dev/null; then
        echo "[FAIL] Test 4: Missing Key (Expected script to abort)"
        fails=$((fails + 1))
    else
        echo "[PASS] Test 4: Missing Key aborted successfully"
        pass=$((pass + 1))
    fi

    rm -f "$mock_ice_in"
    echo "--------------------------------"
    echo "Tests completed: $pass passed, $fails failed."
    return $fails
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # If an argument is passed (e.g., "ice_in"), apply the overrides
    if [[ $# -gt 0 ]]; then
        apply_cice_overrides "$1" || exit 1
    # If no arguments are passed, run the unit tests
    else
        run_tests
    fi
fi

