#!/bin/bash

usage() {
  set +x
  echo
  echo "Usage: $0 -a | -l | -s | -t | -w"
  echo
  echo "  -a  clean all"
  echo "  -l  clean links only"
  echo "  -s  clean shell scripts only"
  echo "  -t  clean tests-dev directory"
  echo "  -w  clean without trace -- i.e. return tests-dev to original state by removing new files from a run and/or -s sync command"
  echo
  exit 1
}

[[ $# -eq 0 ]] && usage

CLEAN_ALL=false
CLEAN_LINKS=false
CLEAN_SHELL_SCRPT=false
CLEAN_TESTS_DIR=false
CLEAN_WITHOUT_TRACE=false

while getopts "alstw" opt; do
  case ${opt} in
    a)
        CLEAN_ALL=true
	;;
    l)
        CLEAN_LINKS=true
	;;
    s)
        CLEAN_SHELL_SCRPT=true
	;;
    t)
        CLEAN_TESTS_DIR=true
	;;
    w)
        CLEAN_ALL=true
        CLEAN_LINKS=true
        CLEAN_SHELL_SCRPT=true
        CLEAN_TESTS_DIR=true
        CLEAN_WITHOUT_TRACE=true
	;;
    *)
        # DEFAULT
        CLEAN_ALL=true
        CLEAN_LINKS=true
        CLEAN_SHELL_SCRPT=true
        CLEAN_TESTS_DIR=false
        CLEAN_WITHOUT_TRACE=false
	;;
  esac
done

#################################################
### CLEAN links
#################################################
if [[ "${CLEAN_ALL,,}" == "true" || "${CLEAN_WITHOUT_TRACE,,}" == "true" || "${CLEAN_LINKS,,}" == "true" ]]; then
  echo "Cleaning links:"
  num_links=$(find . -type l -printf '.' | wc -c)
  echo "  == Deleting $num_links symbolic links =="
  find . -type l -delete
fi

#################################################
### CLEAN shell files
#################################################
if [[ "${CLEAN_ALL,,}" == "true" || "${CLEAN_WITHOUT_TRACE,,}" == "true" || "${CLEAN_SHELL_SCRPT,,}" == "true" ]]; then
   echo "Cleaning shell files:"
   files=(detect_machine.sh rt_utils.sh module-setup.sh compile.sh)
   count=0
   for f in "${files[@]}"; do
       if [[ -e "$f" ]]; then
           ((count++))
           rm "$f"
       fi
   done
   echo "  == Found and deleted $count files =="
fi

#################################################
### CLEAN tests directory
#################################################
if [[ "${CLEAN_ALL,,}" == "true" || "${CLEAN_WITHOUT_TRACE,,}" == "true" || "${CLEAN_TESTS_DIR,,}" == "true" ]]; then
   echo "Cleaning tests directory:"

   shopt -s nullglob

   dirs=(tests exp_conf parm diag_table field_table)
   count=0
   for d in "${dirs[@]}"; do
       for f in test_cases/"$d"/*; do
           case ${d} in
             tests)       rd=../tests/tests ;;
             exp_conf)    rd=../tests/fv3_conf ;;
             parm)        rd=../tests/parm ;;
             diag_table)  rd=../tests/parm/diag_table ;;
             field_table) rd=../tests/parm/field_table ;;
           esac

           [[ -f "$f" ]] || continue
           filename=$(basename "$f")
           target="$rd/$filename"

           if [[ -e "$target" ]]; then
               ((count++))
               rm "$target"
           fi
       done
   done

echo "  == Found and deleted $count files =="
fi

#################################################
### CLEAN the lefotovers
#################################################
if [[ "${CLEAN_WITHOUT_TRACE,,}" == "true" ]]; then
shopt -s nullglob
   echo "Cleaning leftovers:"

   files=(test_changes.list ufs_test_temp.yaml rocoto_w* logs/log_* lock rt_utils.sh)
   count=0
   for f in "${files[@]}"; do
       if [[ -e "$f" ]]; then
           ((count++))
           rm -rf "$f"
       fi
   done
   echo "  == Found and deleted $count files =="
fi
