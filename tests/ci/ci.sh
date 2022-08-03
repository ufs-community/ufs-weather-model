#!/bin/bash
set -eu -x

usage_and_exit() {
  echo
  echo "Note: main purpose of this script is to build Docker base image and containers to run opnReqTest script"
  echo
  echo "Usage: $0 -b <build-image> | -c <create-container>"
  echo "  -b specify Docker image name to build"
  echo "  -c specify Docker container name to create"
  echo
  exit 2
}

IMG_NAME=my-image
CNT_NAME=ci-test-weather
TEST_NAME=""
BUILD_CASE=""

while getopts :b:c: opt; do
  case $opt in
    b)
      IMG_NAME=$OPTARG
      ;;
    c)
      CNT_NAME=$OPTARG
      ;;
  esac
done

CONTAINER_NAME="${CNT_NAME}"
OLD="$(docker ps --all --quiet --filter=name="$CONTAINER_NAME")"
if [ -n "$OLD" ]; then
   docker stop $OLD && docker rm $OLD
fi

docker build --build-arg test_name=$TEST_NAME \
             --build-arg build_case=$BUILD_CASE \
             --no-cache \
             --compress \
            -f Dockerfile -t ${IMG_NAME} ../..

docker create --name "${CNT_NAME}_${BUILD_CASE}" "${IMG_NAME}"

#docker logs --details --timestamps "${TEST_NAME}_${TEST_CASE}"
#exit $(docker inspect "${TEST_NAME}_${TEST_CASE}" --format='{{.State.ExitCode}}')
