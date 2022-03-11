# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM python:3.8-buster

MAINTAINER Anthony S. Castanza <acastanza@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir /src

# for install log files - check here for log files when debugging
RUN mkdir /logs

# install system dependencies
RUN apt-get update --yes
RUN apt-get install build-essential=12.6 --yes | tee /logs/build-essential_install.log
RUN apt-get install libcurl4-gnutls-dev=7.64.0-4+deb10u2 --yes | tee /logs/libcurl4-gnutls_install.log
# libhdf5-serial-dev is virtual pkg of libhdf5-dev 1.10.4+repack-10
# RUN apt-get install libhdf5-serial-dev=1.10.0-patch1+docs-4 --yes | tee /logs/libhdf5-serial-dev.log
RUN apt-get install libhdf5-dev=1.10.4+repack-10 --yes | tee /logs/libhdf5-dev.log
RUN apt-get install libxml2-dev=2.9.4+dfsg1-7+deb10u2 --yes | tee /logs/libxml2_install.log

# install RUST
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# install python with conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && \
    bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# install R dependencies

# install python dependencies
RUN pip install pandas==1.2.5
RUN pip install scipy==1.7.1
RUN pip install anndata==0.7.6
RUN pip install scvelo==0.2.4
RUN pip install git+https://github.com/k3yavi/vpolo.git@v0.3.0


# copy module files
COPY module/* /build/
RUN chmod a+x /build/alevin.vpolo.py

# display software versions
RUN python --version
RUN pip --version

# default command
CMD ["python --version"]

# build using this:
# docker build --rm -t genepattern/convert_alevin:<tag> .
