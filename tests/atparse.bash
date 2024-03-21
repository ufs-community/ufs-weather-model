#! /usr/bin/env bash
function atparse {
    # Ensure "set -x" is inactive, but remember if it was active so we can reset it later.
    local __set_x
    [ -o xtrace ] && __set_x='set -x' || __set_x='set +x'
    set +x
    # Use __ in names to avoid clashing with variables in {var} blocks.
    local __text __before __after __during __had_eoln
    # Allow setting variables on the atparse command line rather than the environment.
    # They will be local variables in this function.
    for __text in "$@" ; do
        if [[ $__text =~ ^([a-zA-Z][a-zA-Z0-9_]*)=(.*)$ ]] ; then
            eval "local ${BASH_REMATCH[1]}"
            eval "${BASH_REMATCH[1]}="'"${BASH_REMATCH[2]}"'
        else
            echo "ERROR: Ignoring invalid argument $__text\n" 1>&2
        fi
    done
    while [[ 1 == 1 ]] ; do
        # Read the next line of text. This will "fail" if no more text
        # is left OR if the last line lacks an end-of-line character.
        __had_eoln=YES
        read -d '' -r __text

        # Stop when "read" reports it is done ($? -ne 0) AND the text is
        # non-empty (! -n "$__text"). This ensures we read the final line
        # even if it lacks an end-of-line character.
        if [[ $? -ne 0 ]] ; then
            if [[ -n "$__text" ]] ; then
                # Text remained, but it had no end-of-line.
                __had_eoln=NO
            else
                break
            fi
        fi
        # Search for strings like @[varname] or @['string'] or @[@]
        while [[ "$__text" =~ ^([^@]*)(@\[[a-zA-Z_][a-zA-Z_0-9]*\]|@\[\'[^\']*\'\]|@\[@\]|@)(.*) ]] ; do
            __before="${BASH_REMATCH[1]}"
            __during="${BASH_REMATCH[2]}"
            __after="${BASH_REMATCH[3]}"
            printf %s "$__before"
            # @['string'] inserts string
            if [[ "$__during" =~ ^@\[\'(.*)\'\]$ ]] ; then
                printf %s "${BASH_REMATCH[1]}"
            # @[@] inserts @
            elif [[ "$__during" == '@[@]' ]] ; then
                printf @
            # @[varname] inserts $varname
            elif [[ "$__during" =~ ^@\[([a-zA-Z_][a-zA-Z_0-9]*)\] ]] ; then
                eval 'printf %s "$'"${BASH_REMATCH[1]}"'"'
            # Unrecognized sequences are inserted verbatim.
            else
                printf '%s' "$__during"
            fi
            # Continue until we run out of text in this line.
            if [[ "$__after" == "$__text" ]] ; then
                break
            fi
            __text="$__after"
        done
        # Print the corrected text
        if [[ "$__had_eoln" == YES ]] ; then
            printf '%s\n' "$__text"
        else
            printf '%s' "$__text"
        fi
    done
    eval "$__set_x"
}

function test_atparse {
    # Note that these cannot be local since they will be invisible
    # to atparse:
    testvar='[testvar]'
    var1='[var1]'
    var2='[var2]'
    var4='[var4]'
    ( cat<<\EOF ; echo -n "line with no end-of-line character [var4] = @[var4]" ) | atparse var3='**'
Nothing special here. = @['Nothing special here.']
[testvar] = @[testvar]
[var1] [var2] = @[var1] @[var2]
** = @[var3]
[var4] == @[var4]
@ = @[@] = @['@']
-n
 eval "export PE$c=\${PE$c:-0}" = @[' eval "export PE$c=\${PE$c:-0}"']
EOF
    echo " ... this text should be on the same line as the line with no end-of-line character"
    echo "After block, \$var3 = \"$var3\" should be empty"
}
