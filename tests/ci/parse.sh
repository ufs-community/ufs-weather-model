#!/bin/bash
set -eu

name_=$(sed -n 1p ci.test)
case_=$(sed -n 2p ci.test)
img_=$(sed -n 3p ci.test)

imd1_=''
[[ $case_ =~ thr || $case_ =~ mpi || $case_ =~ dcp || $case_ =~ rst ]] && imd1_+='std'
[[ $case_ =~ bit ]] && imd1_+=' bit'
[[ $case_ =~ dbg ]] && imd1_+=' dbg'
imd1_=$(echo $imd1_ | sed -e 's/^ *//' -e 's/ *$//')

bld_='{"bld_set":['
for i in $imd1_; do bld_+="\"$i\","; done
bld_=$(echo $bld_ | sed -e 's/,$//')
bld_+=']}'

imd2_=()
test_='{"test_set":['
for i in $case_; do
  test_+="\"$i\","
  [[ $i =~ thr || $i =~ mpi || $i =~ dcp || $i =~ rst ]] && imd2_+=( std )
  [[ $i =~ bit ]] && imd2_+=( bit )
  [[ $i =~ dbg ]] && imd2_+=( dbg )
done

test_=$(echo $test_ | sed -e 's/,$//')
test_+='],"include":['
j=0
for i in $case_; do test_+="{\"test_set\":\"$i\",\"artifact\":\"${imd2_[$j]}\"},"; j=$((j+1)); done
test_=$(echo $test_ | sed -e 's/,$//')
test_+=']}'

echo $name_ $bld_ $test_ $img_
