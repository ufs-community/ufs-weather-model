From noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2e

CMD ["/bin/bash"]

ENV HOME=/home/builder
COPY --chown=builder:builder . $HOME/ufs-weather-model

USER builder
ENV USER=builder
ARG test_name
ARG build_case
ENV test_name=$test_name
ENV build_case=$build_case
ENV test_case=thr
ENV CI_TEST=true
ENV RT_COMPILER=gnu
ENV RT_MACHINE=linux
ENV MACHINE_ID=linux

WORKDIR $HOME/ufs-weather-model/tests
CMD ./opnReqTest -n $test_name -c $build_case
