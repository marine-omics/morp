FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential unzip wget time locales \
  libbz2-dev zlib1g zlib1g-dev liblzma-dev pkg-config libncurses5-dev \
  python3-pip cpanminus curl r-base


RUN cpanm LWP::Simple

WORKDIR /usr/local/


ARG RSEMVER=1.3.3
RUN wget https://github.com/deweylab/RSEM/archive/refs/tags/v${RSEMVER}.tar.gz &&\
  tar -xf  v${RSEMVER}.tar.gz &&\
  cd RSEM-${RSEMVER} &&\
  make && make install &&\
  rm -rf RSEM-${RSEMVER}

# Bowtie2
RUN wget 'https://github.com/BenLangmead/bowtie2/releases/download/v2.5.0/bowtie2-2.5.0-linux-x86_64.zip' &&\
  unzip bowtie2-2.5.0-linux-x86_64.zip && rm bowtie2-2.5.0-linux-x86_64.zip

ENV PATH="${PATH}:/usr/local/bowtie2-2.5.0-linux-x86_64"


# Fastp
WORKDIR /usr/local/bin/
RUN wget http://opengene.org/fastp/fastp.0.23.1 && \
    mv fastp.0.23.1 fastp &&\
    chmod a+x ./fastp

# MultiQC
RUN pip install multiqc


# Cleanup apt package lists to save space
RUN rm -rf /var/lib/apt/lists/*

