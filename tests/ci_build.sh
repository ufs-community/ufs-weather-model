#!/bin/bash

export RT_MACHINE=linux

export RT_COMPILER=gnu

../opnReqTest -n $test_name -c $build_case -z 
