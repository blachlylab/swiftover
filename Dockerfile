FROM ubuntu:20.04

LABEL version="1.0"
LABEL description="swiftover -- ultrafast liftover"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y wget build-essential git

# ldc2
RUN apt-get install -y libxml2
RUN wget https://github.com/ldc-developers/ldc/releases/download/v1.24.0/ldc2-1.24.0-linux-x86_64.tar.xz
RUN tar xfvJ ldc2-1.24.0-linux-x86_64.tar.xz
ENV PATH="${PWD}/ldc2-1.24.0-linux-x86_64/bin/:${PATH}"

# htslib
RUN apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
RUN tar xfvj htslib-1.11.tar.bz2
RUN cd htslib-1.11 && ./configure && make && make install
RUN ldconfig

# swiftover
RUN git clone https://github.com/blachlylab/swiftover.git
RUN cd swiftover && dub build -b=release-nobounds
ENV PATH="${PWD}/swiftover:${PATH}"

