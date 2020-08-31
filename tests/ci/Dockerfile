From minsukjinoaa/ci-test-base:ubuntu20.04-v3
#From minsukjinoaa/ci-test-base:centos7-v3

ENV HOME=/home/tester
COPY --chown=tester:tester . $HOME/ufs-weather-model

USER tester
ENV USER=tester
ARG test_name
ARG build_case
ENV test_name=$test_name
ENV build_case=$build_case
ENV test_case=
ENV CI_TEST=true
ENV NEMS_COMPILER=gnu
ENV NEMS_MACHINE=linux.gnu

WORKDIR $HOME/ufs-weather-model/tests
RUN ./utest -n $test_name -c $build_case -z
CMD ./utest -n $test_name -c $test_case -x
