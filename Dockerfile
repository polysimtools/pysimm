FROM python:3.7-buster

COPY --from=registry.gitlab.com/mefortunato/docker-images/lammps:7Aug2019 /usr/local/bin/lmp_mpi /usr/local/bin/lmp_mpi

COPY . /usr/local/pysimm

RUN apt-get update && \
    apt-get install -y libopenmpi-dev openmpi-bin && \
    pip install -r /usr/local/pysimm/requirements.txt
