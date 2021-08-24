FROM thomasrobinson/centos7-netcdff:4.5.3-c4.7.4-gcc-mpich-slurm
## Dockerfile used to create CM4

## Set up spack
RUN . /opt/spack/share/spack/setup-env.sh
## Make the CM4 directory
RUN mkdir -p /opt/CM4
## Build the CM4 from github
RUN git clone --recursive https://github.com/NOAA-GFDL/CM4.git -b 2021.02 \
    && cd CM4/exec \ 
    && make gcc=on HDF_INCLUDE=-I/opt/hdf5/include SH=sh CLUBB=off \
    && cp cm4.x /opt/CM4 \
    && make clean_all
## Add the CM4 executable to the path
ENV PATH=/opt/CM4/:${PATH}
## Add permissions to the CM4
RUN chmod 777 /opt/CM4/cm4.x

