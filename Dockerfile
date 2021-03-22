FROM python:3.7-buster

COPY --from=registry.gitlab.com/mefortunato/docker-images/lammps:7Aug2019 /usr/local/bin/lmp_mpi /usr/local/bin/lmp_mpi

COPY . /usr/local/pysimm


RUN apt-get update && \
    apt-get install -y libopenmpi-dev openmpi-bin && \
    pip install -r /usr/local/pysimm/requirements.txt && \
    pip install -e /usr/local/pysimm

ENV LAMMPS_EXEC="lmp_mpi"

RUN apt-get update && apt-get install -y gfortran

RUN cd /usr/local/lib && \
    git clone https://github.com/SarkisovGroup/PoreBlazer.git && \
    cd PoreBlazer/src/ && \
    make -f Makefile_gfort && \
    cp poreblazer.exe /usr/local/bin && \
    cd ../../ && \
    rm -rf PoreBlazer

RUN TB_NAME=zeopp.v0p3.tar.gz && DIR_NAME=zeopp.v0p3 && \
    cd /usr/local/lib && \
    curl http://www.zeoplusplus.org/zeo++-0.3.tar.gz --output $TB_NAME && \
    tar -xf $TB_NAME && \
    rm $TB_NAME && mv zeo++-0.3 $DIR_NAME && \ 
    cd $DIR_NAME/voro++/src && make && \
    cd ../.. && make && \ 
    cp network /usr/local/bin && \ 
    rm -rf $DIR_NAME

ENV ZEOpp_EXEC="network"