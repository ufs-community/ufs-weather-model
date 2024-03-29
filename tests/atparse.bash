#! /usr/bin/env bash
function atparse {
    # Usage:
    #     source atparse.bash # defines the "atparse" function; only do this once
    #     atparse [ var1=value1 [ var2=value2 [...] ] ] < input_file > output_file
    # This function filters text from stdin to stdout. It scans for text sequences like:
    #    @[varname]
    # And replaces them with the value of the corresponding ${varname} variable.
    # You can provide variables that are not set in bash by providing them on the command line.
    # If set -u is enabled, it will exit the process when a variable is empty or undefined via set -u.

    # Use __ in names to avoid clashing with variables in {var} blocks.
    local __text # current line of text being parsed, or the current command-line argument being parsed
    local __before # all text before the next @[...] option
    local __after # all text after the next @[...] option
    local __during # the contents of the @[...] option, including the @[ and ]
    local __set_x=":" # will be "set -x" if the calling script had that option enabled
    local __set_u=":" # will be "set -u" if the calling script had that option enabled
    local __set_e=":" # will be "set -e" if the calling script had that option enabled
    local __abort_on_undefined=NO # YES = script should abort if a variable is undefined, NO otherwise

    # Ensure "set -x -e -u" are all inactive, but remember if they
    # were active so we can reset them later.
    if [[ -o xtrace ]] ; then
        __set_x="set -x"
    fi
    if [[ -o errexit ]] ; then
        __set_e="set -e"
    fi
    if [[ -o nounset ]] ; then
        __set_u="set -u"
        __abort_on_undefined=YES
    fi
    set +eux

    # Allow setting variables on the atparse command line rather than the environment.
    # They will be local variables in this function.
    for __text in "$@" ; do
        if [[ $__text =~ ^([a-zA-Z][a-zA-Z0-9_]*)=(.*)$ ]] ; then
            eval "local ${BASH_REMATCH[1]}"
            eval "${BASH_REMATCH[1]}="'"${BASH_REMATCH[2]}"'
        else
            echo "ERROR: Ignoring invalid argument $__text" 1>&2
        fi
    done

    # Loop over all lines of text.
    while [[ 1 == 1 ]] ; do
        # Read the next line of text. This will "fail" if no more text
        # is left OR if the last line lacks an end-of-line character.
        read -d '' -r __text

        # Stop when "read" reports it is done ($? -ne 0) AND the text is
        # non-empty (! -n "$__text"). This ensures we read the final line
        # even if it lacks an end-of-line character.
        if [[ $? -ne 0 ]] ; then
            if [[ -n "$__text" ]] ; then
                # Text remained, but it had no end-of-line.
                :
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
                # Flag unknown variables at this step only.
                if [[ ${__abort_on_undefined} == YES ]] ; then
                    set -u
                fi
                eval 'printf %s "$'"${BASH_REMATCH[1]}"'"'
                if [[ ${__abort_on_undefined} == YES ]] ; then
                    set +u
                fi
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
        printf '%s\n' "$__text"
    done

    # Restore the calling script's shell options.
    eval "$__set_x"
    eval "$__set_u"
    eval "$__set_e"
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
@[undefined_variable_that_should_exit_script_if_set_minus_u_is_used]
-n
 eval "export PE$c=\${PE$c:-0}" = @[' eval "export PE$c=\${PE$c:-0}"']
EOF
    echo " ... this text should be on the same line as the line with no end-of-line character"
    echo "After block, \$var3 = \"$var3\" should be empty"
}
