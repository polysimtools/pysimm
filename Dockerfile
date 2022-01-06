FROM debian:buster as build

ARG PACKAGES="molecule extra-molecule kspace misc qeq class2 manybody"
ARG BIN_DIR="/usr/local/bin"

RUN apt-get update && \
    apt-get install -y make git g++ libopenmpi-dev openmpi-bin curl


RUN curl --silent "https://api.github.com/repos/lammps/lammps/releases/latest"| grep '"tag_name":'|sed -E 's/.*"([^"]+)".*/\1/' > /root/lmp_version && \
    git clone -b $(cat /root/lmp_version) https://github.com/lammps/lammps.git && \
    cd lammps/src && \
    for PACKAGE in $PACKAGES; do make yes-$PACKAGE; done && \
    make -j2 mpi && \
    cp lmp_mpi $BIN_DIR && \
    cd ../../ && \
    rm -rf lammps /root/lmp_version


RUN cd /usr/local/lib && \
    git clone -b v1.2.5 https://github.com/MaginnGroup/Cassandra.git && \
    cd Cassandra/Src/ && \
    make clean && make -f Makefile.gfortran.openMP && \
    cp cassandra_gfortran_openMP.exe $BIN_DIR/cs_gfort_omp.exe && \
    cd ../../ && rm -rf Cassandra


RUN cd /usr/local/lib && \
    git clone https://github.com/SarkisovGroup/PoreBlazer.git && \
    cd PoreBlazer/src/ && \
    make -f Makefile_gfort && \
    cp poreblazer.exe $BIN_DIR && \
    cd ../../ && \
    rm -rf PoreBlazer


RUN TB_NAME=zeopp.v0p3.tar.gz && DIR_NAME=zeopp.v0p3 && \
    cd /usr/local/lib && \
    curl http://www.zeoplusplus.org/zeo++-0.3.tar.gz --output $TB_NAME && \
    tar -xf $TB_NAME && \
    rm $TB_NAME && mv zeo++-0.3 $DIR_NAME && \
    cd $DIR_NAME/voro++/src && make && \
    cd ../.. && make && \
    cp network $BIN_DIR && \
    rm -rf $DIR_NAME


FROM python:3.8-buster

COPY --from=build $BIN_DIR/* /usr/local/bin/

RUN apt-get update && \
    apt-get install -y libopenmpi-dev openmpi-bin vim mc

COPY . /usr/local/pysimm
RUN pip install -r /usr/local/pysimm/requirements.txt && \
    pip install -e /usr/local/pysimm

ENV LAMMPS_EXEC="lmp_mpi"
ENV CASSANDRA_EXEC="cs_gfort_omp.exe"
ENV ZEOpp_EXEC="network"
