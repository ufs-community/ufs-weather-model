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

IMG_NAME=ci-test-weather
BUILD="false"
RUN="false"
TEST_NAME=""
TEST_CASE=""

# Read in TEST_NAME
TEST_NAME=$(sed -n 1p ci.test)

while getopts :c:br opt; do
  case $opt in
    c)
      TEST_CASE=$OPTARG
      ;;
    b)
      BUILD="true"
      ;;
    r)
      RUN="true"
      ;;
  esac
done

if [ $BUILD = "true" ] && [ $RUN = "true" ]; then
  echo "Specify either build (-b) or run (-r) option, not both"
  exit 2
fi

if [ $BUILD = "true" ] && [ Q$TEST_CASE = Q"" ]; then
  echo "Build option (-b) should accompany TEST_CASE option (-c)"
  exit 2
fi

if [ $RUN = "true" ] && [ Q$TEST_CASE != Q"" ]; then
  echo "Run option (-r) should not accompany TEST_CASE option (-c)"
  exit 2
fi

if [ $BUILD = "true" ]; then
  sed -i -e '/affinity.c/d' ../CMakeLists.txt

  docker build --build-arg test_name=$TEST_NAME \
               --build-arg test_case=$TEST_CASE \
               --no-cache \
               --squash --compress \
               -f ../Dockerfile -t ${IMG_NAME} ..
  exit $?

elif [ $RUN == "true" ]; then
  docker run -d ${IMG_NAME}
  echo 'cache,rss,shmem' >memory_stat
  sleep 3
  containerID=$(docker ps -q --no-trunc)
  check_memory_usage $containerID >>memory_stat &
  docker logs -f $containerID
  exit $(docker inspect $containerID --format='{{.State.ExitCode}}')

fi
