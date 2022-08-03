#!/bin/bash
set -eu -x

usage_and_exit() {
  echo
  echo "Note: main purpose of this script is to build Docker base image to run opnReqTest script"
  echo
  echo "Usage: $0 -b <build-image>"
  echo "  -b specify Docker image name to build"
  echo
  exit 2
}

IMG_NAME=ci-test-weather
TEST_NAME=""
RUN_CASE=""

while getopts :b: opt; do
  case $opt in
    b)
      IMG_NAME=$OPTARG
      ;;
  esac
done

#CONTAINER_NAME="${CNT_NAME}"
#OLD="$(docker ps --all --quiet --filter=name="$CONTAINER_NAME")"
#if [ -n "$OLD" ]; then
#   docker stop $OLD && docker rm $OLD
#fi

docker build --build-arg test_name=$TEST_NAME \
             --build-arg build_case=$RUN_CASE \
             --no-cache \
             --compress \
            -f Dockerfile -t ${IMG_NAME} ../..

#docker create --name "${CNT_NAME}" "${IMG_NAME}"

#docker logs --details --timestamps "${TEST_NAME}_${TEST_CASE}"
#exit $(docker inspect "${TEST_NAME}_${TEST_CASE}" --format='{{.State.ExitCode}}')
