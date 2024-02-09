#! /usr/bin/env bash
function atparse {
    local __set_x
    [ -o xtrace ] && __set_x='set -x' || __set_x='set +x'
    set +x
    # Use __ in names to avoid clashing with variables in {var} blocks.
    local __text __before __after __during
    for __text in "$@" ; do
        if [[ $__text =~ ^([a-zA-Z][a-zA-Z0-9_]*)=(.*)$ ]] ; then
            eval "local ${BASH_REMATCH[1]}"
            eval "${BASH_REMATCH[1]}="'"${BASH_REMATCH[2]}"'
        else
            echo "ERROR: Ignoring invalid argument $__text\n" 1>&2
        fi
    done
    while IFS= read -r __text ; do
        while [[ "$__text" =~ ^([^@]*)(@\[[a-zA-Z_][a-zA-Z_0-9]*\]|@\[\'[^\']*\'\]|@\[@\]|@)(.*) ]] ; do
            __before="${BASH_REMATCH[1]}"
            __during="${BASH_REMATCH[2]}"
            __after="${BASH_REMATCH[3]}"
#            printf 'PARSE[%s|%s|%s]\n' "$__before" "$__during" "$__after"
            printf %s "$__before"
            if [[ "$__during" =~ ^@\[\'(.*)\'\]$ ]] ; then
                printf %s "${BASH_REMATCH[1]}"
            elif [[ "$__during" == '@[@]' ]] ; then
                printf @
            elif [[ "$__during" =~ ^@\[([a-zA-Z_][a-zA-Z_0-9]*)\] ]] ; then
                eval 'printf %s "$'"${BASH_REMATCH[1]}"'"'
            else
                printf '%s' "$__during"
            fi
            if [[ "$__after" == "$__text" ]] ; then
                break
            fi
            __text="$__after"
        done
        printf '%s\n' "$__text"
    done
    eval "$__set_x"
}

function test_atparse {
    # Note that these cannot be local since they will be invisible
    # to atparse:
    testvar='[testvar]'
    var1='[var1]'
    var2='[var2]'
    cat<<\EOF | atparse var3='**'
Nothing special here. = @['Nothing special here.']
[testvar] = @[testvar]
[var1] [var2] = @[var1] @[var2]
** = @[var3]
@ = @[@] = @['@']
-n
 eval "export PE$c=\${PE$c:-0}" = @[' eval "export PE$c=\${PE$c:-0}"']
EOF
    echo "After block, \$var3 = \"$var3\" should be empty"
}
