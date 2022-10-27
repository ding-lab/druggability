FROM ubuntu:18.04
MAINTAINER R. Jay Mashl "rmashl@wustl.edu"

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN  apt-get update  &&  apt-get install -y --no-install-recommends \
  ca-certificates \
  git \
  wget \
  && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh \
    && mkdir -p /root/.conda \
    && bash Miniconda3-py37_4.12.0-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py37_4.12.0-Linux-x86_64.sh

RUN activate

RUN pip install pandas==1.0.5

WORKDIR /usr/local
RUN git clone --recurse-submodules https://github.com/ding-lab/druggability.git
