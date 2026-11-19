#!/bin/bash
set -eux

cat ufs_trace_*.trace > all.traces
sed -i '$ s/.$//' all.traces
echo '[' > out.trace
cat all.traces >> out.trace
echo ']' >> out.trace
