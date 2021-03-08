From noaaemc/ubuntu-hpc:v1.1

ENV HOME=/home/builder
COPY --chown=builder:builder . $HOME/ufs-weather-model

USER builder
ENV USER=builder
ARG test_name
ARG build_case
ENV test_name=$test_name
ENV build_case=$build_case
ENV test_case=
ENV CI_TEST=true
ENV RT_COMPILER=gnu
ENV RT_MACHINE=linux.gnu

WORKDIR $HOME/ufs-weather-model/tests
RUN ./utest -n $test_name -c $build_case -z
CMD ./utest -n $test_name -c $test_case -x
