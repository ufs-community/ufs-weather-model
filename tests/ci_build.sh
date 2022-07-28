#!/bin/bash

export RT_MACHINE=linux

export RT_COMPILER=gnu

/bin/bash ../modulefiles/ufs_linux.gnu && MACHINE_ID=linux && ./opnReqTest -n control -c thr -z
