#!/bin/bash
set -eu

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
  echo "Note: main purpose of this script is to interface between CI automation and utest script"
  echo "and therefore, direct invocation via CLI may result in unexpected behavior"
  echo
  echo "Usage: $0 -b <build-case> | -r <test-case>"
  echo "  -b specify cases to build: comma-separated list of any combination of std,big,dbg"
  echo "  -r specify tests to run: comma-separated list of any combination of std,thr,mpi,dcp,rst,bit,dbg"
  echo
  exit 2
}

IMG_NAME=$(sed -n 3p ci.test)
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
TEST_NAME=${TEST_NAME:-$(sed -n 1p ci.test)}
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

  sed -i -e '/affinity.c/d' ../../CMakeLists.txt

  sudo docker build --build-arg test_name=$TEST_NAME \
                    --build-arg build_case=$BUILD_CASE \
                    --no-cache \
                    --squash --compress \
                    -f Dockerfile -t ${IMG_NAME} ../..
  exit $?

elif [ $RUN == "true" ]; then

  sudo docker run -d --rm -v DataVolume:/tmp minsukjinoaa/fv3-input-data:develop-20200713
  sudo docker run -d -e test_case=${TEST_CASE} -v DataVolume:/home/tester/data/NEMSfv3gfs/develop-20200713 --name my-container ${IMG_NAME}

  echo 'cache,rss,shmem' >memory_stat
  sleep 3
  containerID=$(sudo docker ps -q --no-trunc)
  check_memory_usage $containerID >>memory_stat &

  sudo docker logs -f $containerID
  exit $(sudo docker inspect $containerID --format='{{.State.ExitCode}}')

fi
