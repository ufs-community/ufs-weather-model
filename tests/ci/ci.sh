#!/bin/bash
set -eu -x

check_memory_usage() {
  dirName=/sys/fs/cgroup/memory/docker/$1
  # /sys/fs/cgroup/memory/actions_job/${containerID}
  set +x
  while [ -d $dirName ] ; do
    awk '/(^cache |^rss |^shmem )/' $dirName/memory.stat | cut -f2 -d' ' | paste -s -d,
    sleep 1
  done
  set -x
}

usage_and_exit() {
  echo
  echo "Note: main purpose of this script is to interface between CI automation and opnReqTest script"
  echo "and therefore, direct invocation via CLI may result in unexpected behavior"
  echo
  echo "Usage: $0 -b <build-case> | -r <test-case>"
  echo "  -b specify cases to build: comma-separated list of any combination of std,big,dbg"
  echo "  -r specify tests to run: comma-separated list of any combination of std,thr,mpi,dcp,rst,bit,dbg"
  echo
  exit 2
}

IMG_NAME=ci-test-weather
BUILD="false"
RUN="false"
TEST_NAME=""
BUILD_CASE=""
TEST_CASE=""

while getopts :b:r:n: opt; do
  case $opt in
    b)
      BUILD="true"
      BUILD_CASE=$OPTARG
      ;;
    r)
      RUN="true"
      TEST_CASE=$OPTARG
      ;;
    n)
      TEST_NAME=$OPTARG
      ;;
  esac
done

# Read in TEST_NAME if not passed on
echo "test name is ${TEST_NAME}"

if [ $BUILD = "true" ] && [ $RUN = "true" ]; then
  echo "Specify either build (-b) or run (-r) option, not both"
  usage_and_exit
fi

if [ $BUILD = "false" ] && [ $RUN = "false" ]; then
  echo "Specify either build (-b) or run (-r) option"
  usage_and_exit
fi

if [ $BUILD = "true" ]; then

  docker build --build-arg test_name=$TEST_NAME \
               --build-arg build_case=$BUILD_CASE \
               --no-cache \
               --compress \
               -f Dockerfile -t ${IMG_NAME} ../..

  docker create --name "${TEST_NAME}_${BUILD_CASE}" "${IMG_NAME}"
  #docker cp -a "${TEST_NAME}_${BUILD_CASE}:/home/builder/ufs-weather-model/tests/fv3.tar.gz" /home/maintuser/workspace/sandbox/ufs-wm-ci-debug/
  #docker rm "${TEST_NAME}_${BUILD_CASE}"

elif [ $RUN == "true" ]; then

  CONTAINER_NAME="${TEST_NAME}_${TEST_CASE}"
  OLD="$(docker ps --all --quiet --filter=name="$CONTAINER_NAME")"
  if [ -n "$OLD" ]; then
    docker stop $OLD && docker rm $OLD
  fi

  docker pull noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2e

  #docker volume rm -f DataVolume >/dev/null &&
  docker run -d --rm -v DataVolume:/tmp noaaepic/input-data:20220414

  docker create -u builder -e "CI_TEST=true" -e "USER=builder" \
                -e "RT_MACHINE=linux.gnu" -e "RT_COMPILER=gnu" \
                -w "/home/builder/ufs-weather-model/tests" \
                -v DataVolume:/home/builder/data/NEMSfv3gfs \
                --shm-size=512m --name "${TEST_NAME}_${TEST_CASE}" noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2e \
                /bin/bash -c "./opnReqTest -n ${TEST_NAME} -c ${TEST_CASE} -z  && ./opnReqTest -n ${TEST_NAME} -c ${TEST_CASE} -x"

  cd $GITHUB_WORKSPACE
  docker cp . "${TEST_NAME}_${TEST_CASE}:/home/builder/ufs-weather-model"
  docker start "${TEST_NAME}_${TEST_CASE}"

  sleep 3

  docker logs --details --timestamps "${TEST_NAME}_${TEST_CASE}"

  #exit $(docker inspect "${TEST_NAME}_${TEST_CASE}" --format='{{.State.ExitCode}}')
fi
