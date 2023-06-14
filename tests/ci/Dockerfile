From noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2

CMD ["/bin/bash"]

ENV HOME=/home/builder
COPY --chown=builder:builder . $HOME/ufs-weather-model

USER builder
ENV USER=builder
ARG test_name
ARG run_case
ENV test_name=$test_name
ENV run_case=$run_case
ENV CI_TEST=true
ENV RT_COMPILER=gnu
ENV RT_MACHINE=linux
ENV MACHINE_ID=linux
# ACCNR CAN BE ANYTHING; DOCKER DOES NOT NEED THIS
ENV ACCNR=nems 

WORKDIR $HOME/ufs-weather-model/tests
CMD ./opnReqTest -n $test_name -c $run_case
