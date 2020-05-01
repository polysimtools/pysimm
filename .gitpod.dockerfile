FROM gitpod/workspace-full

COPY --from=registry.gitlab.com/mefortunato/docker-images/lammps:7Aug2019 /usr/local/bin/lmp_mpi /usr/local/bin/lmp_mpi

USER root

COPY requirements.txt requirements.txt 

RUN apt-get update && \
    apt-get install -y libopenmpi-dev openmpi-bin && \
    pip3 install -r requirements.txt
    
USER gitpod

ENV LAMMPS_EXEC='lmp_mpi'
