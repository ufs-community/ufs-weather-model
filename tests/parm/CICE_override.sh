#!/bin/bash
set -eu

# =====================================================================
# FUNCTION: apply_namelist_overrides
# 
# Description: 
#   A unified utility to override individual keys AND/OR replace entire 
#   Fortran namelist blocks in a target file. 
#
# Features:
#   - Modifies the file in-place while preserving permissions/symlinks.
#   - Safely handles Fortran case-insensitivity.
#   - Fails fast and breaks execution if keys or blocks are missing.
#   - Reads from standard input (stdin).
#   - Strictly preserves chronological sequence (interleaving safe).
#
# Example usage:
#   ./namelist_override.sh ice_in < set_nml.bgc
#   echo "dt = 1800.0" | ./namelist_override.sh ice_in
# =====================================================================

apply_namelist_overrides() {
    local target_file="$1"

    if [[ -z "$target_file" || ! -f "$target_file" ]]; then
        echo "FATAL ERROR: Missing or invalid target file: '${target_file}'" >&2
        return 1
    fi

    # Create a temporary working directory
    local tmp_dir
    tmp_dir=$(mktemp -d) || { echo "FATAL ERROR: Failed to create temp directory." >&2; return 1; }
    
    # 1. Read all input from stdin into a master file
    cat > "${tmp_dir}/master_input"

    # 2. Smart Parser: Split input into sequentially numbered chunks
    awk -v outdir="$tmp_dir" '
    BEGIN { in_block = 0; seq = 0; last_type = "" }
    { line_lower = tolower($0) }
    
    # Skip empty lines or pure comments (# or !) outside blocks
    in_block == 0 && (line_lower ~ /^[ \t]*$/ || line_lower ~ /^[ \t]*[#!]/) { next }
    
    # Match the start of a namelist group (e.g., &grid_nml)
    line_lower ~ /^[ \t]*&[a-z0-9_]+/ && in_block == 0 {
        match(line_lower, /&[a-z0-9_]+/)
        group = substr(line_lower, RSTART+1, RLENGTH-1)
        seq++
        last_type = "block"
        file = sprintf("%s/%04d_block_%s.txt", outdir, seq, group)
        in_block = 1
    }
    
    # Inside a block, dump line to the specific block file
    in_block == 1 {
        print $0 > file
    }
    
    # Match the end of a namelist group (/)
    line_lower ~ /^[ \t]*\/[ \t\r\n]*$/ && in_block == 1 {
        close(file)
        in_block = 0
        next
    }
    
    # If outside a block and contains an "=" sign, it is an individual key override
    in_block == 0 && $0 ~ /=/ {
        # Group contiguous key overrides into the same numbered file
        if (last_type != "key") {
            seq++
            last_type = "key"
            file = sprintf("%s/%04d_keys.txt", outdir, seq)
        }
        print $0 >> file
    }
    ' "${tmp_dir}/master_input"

    # 3. Unified Execution Engine: Process sequentially sorted files
    for override_file in "${tmp_dir}"/*.txt; do
        # Skip if glob didn't expand
        [[ -e "$override_file" ]] || continue
        
        local filename
        filename=$(basename "$override_file")
        
        # ---------------------------------------------------------
        # HANDLE BLOCK REPLACEMENT
        # ---------------------------------------------------------
        if [[ "$filename" == *_block_* ]]; then
            local group_name
            group_name=${filename#*_block_}
            group_name=${group_name%.txt}
            
            local exact_match
            exact_match=$(awk -v grp="$group_name" '
                BEGIN { gl = tolower(grp) } 
                { ll = tolower($0); if (ll ~ "^[ \t]*&" gl "([ \t\r\n]|$)") { print $0; exit } }
            ' "$target_file" || true)

            if [[ -z "$exact_match" ]]; then
                echo "FATAL ERROR: Namelist block '&${group_name}' not found in '${target_file}'." >&2
                rm -rf "$tmp_dir"
                return 1
            fi

            local tmp_out
            tmp_out=$(mktemp)
            
            awk -v repl_file="$override_file" -v group="$group_name" '
            BEGIN { in_block = 0; group_lower = tolower(group) }
            { line_lower = tolower($0) }
            in_block == 0 && line_lower ~ "^[ \t]*&" group_lower "([ \t\r\n]|$)" {
                in_block = 1
                while ((getline repl_line < repl_file) > 0) { print repl_line }
                next
            }
            in_block == 1 && line_lower ~ "^[ \t]*/[ \t\r\n]*$" { in_block = 0; next }
            in_block == 1 { next }
            in_block == 0 { print }
            ' "$target_file" > "$tmp_out"
            
            cat "$tmp_out" > "$target_file"
            rm -f "$tmp_out"

        # ---------------------------------------------------------
        # HANDLE INDIVIDUAL KEYS
        # ---------------------------------------------------------
        elif [[ "$filename" == *_keys.txt ]]; then
            local tmp_key_out
            tmp_key_out=$(mktemp)

            while read -r line; do
                local key
                key=$(echo "$line" | cut -d'=' -f1 | awk '{print $1}')
                local value
                value=$(echo "$line" | cut -d'=' -f2- | sed 's/^[[:space:]]*//')
                
                [[ -z "$key" || -z "$value" ]] && continue
                
                local exact_match
                exact_match=$(grep -i -m 1 "^[[:space:]]*${key}[[:space:]]*=" "$target_file" || true)
                
                if [[ -z "$exact_match" ]]; then
                    echo "FATAL ERROR: Override key '${key}' not found in '${target_file}'. Misspelled?" >&2
                    rm -f "$tmp_key_out"
                    rm -rf "$tmp_dir"
                    return 1
                fi
                
                local exact_key
                exact_key=$(echo "$exact_match" | cut -d'=' -f1 | awk '{print $1}')
                
                # Standard POSIX sed used here for maximum cross-platform compatibility
                sed "s/^[[:space:]]*${exact_key}[[:space:]]*=.*/  ${exact_key} = ${value}/" "$target_file" > "$tmp_key_out"
                cat "$tmp_key_out" > "$target_file"
                
            done < "$override_file"
            
            rm -f "$tmp_key_out"
        fi
    done

    # Clean up the working directory safely
    rm -rf "$tmp_dir"
}

# =====================================================================
# UNIT TESTS
# =====================================================================
run_tests() {
    echo "Running unit tests for apply_namelist_overrides..."
    local fails=0
    local pass=0

    local mock_target
    mock_target=$(mktemp)
    
    reset_mock_target() {
        {
            echo "&setup_nml"
            echo "  days_per_year = 365"
            echo "  dt = 3600.0"
            echo "/"
            echo "&grid_nml"
            echo "  grid_format = 'nc'"
            echo "  ncat = 1"
            echo "/"
        } > "$mock_target"
    }

    local mock_repl
    mock_repl=$(mktemp)

    # ---------------------------------------------------------
    # TEST 1: Unified - Block + Individual Keys
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "dt = 1800.0"
        echo "&grid_nml"
        echo "  grid_format = 'bin'"
        echo "  ncat = 5"
        echo "/"
        echo "days_per_year = 360"
    } > "$mock_repl"

    apply_namelist_overrides "$mock_target" < "$mock_repl"
    
    if grep -qF "ncat = 5" "$mock_target" && grep -qF "dt = 1800.0" "$mock_target"; then
        echo "[PASS] Test 1: Unified Overrides (Blocks + Keys) successful"
        pass=$((pass + 1))
    else
        echo "[FAIL] Test 1: Unified Overrides Failed"
        fails=$((fails + 1))
    fi

    # ---------------------------------------------------------
    # TEST 2: Fail Fast (Missing Key)
    # ---------------------------------------------------------
    reset_mock_target
    if echo "bad_key = 99" | apply_namelist_overrides "$mock_target" 2>/dev/null; then
        echo "[FAIL] Test 2: Missing Key (Expected script to abort)"
        fails=$((fails + 1))
    else
        echo "[PASS] Test 2: Missing Key aborted successfully"
        pass=$((pass + 1))
    fi

    # ---------------------------------------------------------
    # TEST 3: Fail Fast (Missing Block)
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "&fake_nml"
        echo "  ncat = 5"
        echo "/"
    } > "$mock_repl"
    
    if apply_namelist_overrides "$mock_target" < "$mock_repl" 2>/dev/null; then
        echo "[FAIL] Test 3: Missing Block (Expected script to abort)"
        fails=$((fails + 1))
    else
        echo "[PASS] Test 3: Missing Block aborted successfully"
        pass=$((pass + 1))
    fi

    # ---------------------------------------------------------
    # TEST 4: Interleaving - Block Overridden by Key
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "&grid_nml"
        echo "  grid_format = 'bin'"
        echo "  ncat = 5"
        echo "/"
        echo "ncat = 99"
    } > "$mock_repl"
    
    apply_namelist_overrides "$mock_target" < "$mock_repl"
    
    if grep -qF "ncat = 99" "$mock_target" && ! grep -qF "ncat = 5" "$mock_target"; then
        echo "[PASS] Test 4: Interleaving (Block -> Key) successful"
        pass=$((pass + 1))
    else
        echo "[FAIL] Test 4: Interleaving (Block -> Key) failed"
        fails=$((fails + 1))
    fi

    # ---------------------------------------------------------
    # TEST 5: Interleaving - Key Overridden by Block
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "dt = 111.0"
        echo "&setup_nml"
        echo "  days_per_year = 365"
        echo "  dt = 222.0"
        echo "/"
    } > "$mock_repl"
    
    apply_namelist_overrides "$mock_target" < "$mock_repl"
    
    if grep -qF "dt = 222.0" "$mock_target" && ! grep -qF "dt = 111.0" "$mock_target"; then
        echo "[PASS] Test 5: Interleaving (Key -> Block) successful"
        pass=$((pass + 1))
    else
        echo "[FAIL] Test 5: Interleaving (Key -> Block) failed"
        fails=$((fails + 1))
    fi

    # ---------------------------------------------------------
    # TEST 6: Interleaving - The "Ping-Pong" Sequence
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "dt = 111.0"
        echo "&setup_nml"
        echo "  dt = 222.0"
        echo "/"
        echo "dt = 333.0"
        echo "&setup_nml"
        echo "  dt = 444.0"
        echo "/"
    } > "$mock_repl"
    
    apply_namelist_overrides "$mock_target" < "$mock_repl"
    
    if grep -qF "dt = 444.0" "$mock_target" && ! grep -qF "333.0" "$mock_target" && ! grep -qF "222.0" "$mock_target"; then
        echo "[PASS] Test 6: Interleaving (Ping-Pong Chronology) successful"
        pass=$((pass + 1))
    else
        echo "[FAIL] Test 6: Interleaving (Ping-Pong Chronology) failed"
        fails=$((fails + 1))
    fi

    # ---------------------------------------------------------
    # TEST 7: Chaotic Interspersed Inputs & Spacing
    # ---------------------------------------------------------
    reset_mock_target
    {
        echo "   dT = 720.0   "
        echo "&GRID_nml"
        echo "  grid_format = 'bin'"
        echo "  ncat = 5"
        echo "/"
        echo "  DAYS_per_YEAR= 366"
    } > "$mock_repl"

    apply_namelist_overrides "$mock_target" < "$mock_repl"
    
    if grep -qF "dt = 720.0" "$mock_target" && grep -qF "days_per_year = 366" "$mock_target" && grep -qF "ncat = 5" "$mock_target"; then
        echo "[PASS] Test 7: Chaotic interspersed inputs handled correctly"
        pass=$((pass + 1))
    else
        echo "[FAIL] Test 7: Chaotic interspersed inputs failed"
        fails=$((fails + 1))
    fi

    rm -f "$mock_target" "$mock_repl"
    echo "---------------------------------------------------------"
    if [[ $fails -eq 0 ]]; then
        echo "✅ All tests completed successfully: $pass passed, $fails failed."
    else
        echo "❌ Tests completed with errors: $pass passed, $fails failed."
    fi
    return $fails
}

# =====================================================================
# SCRIPT EXECUTION ROUTING
# =====================================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -gt 0 ]]; then
        apply_namelist_overrides "$1" || exit 1
    else
        run_tests
    fi
fi
