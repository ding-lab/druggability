FROM ubuntu:18.04
MAINTAINER R. Jay Mashl "rmashl@wustl.edu"

RUN  apt-get update  &&  apt-get install -y --no-install-recommends \
  python3 \
  && rm -rf /var/lib/apt/lists/*


COPY  *.py /usr/local/
COPY  druggability_databases  /usr/local/druggability_databases

CMD ["/bin/bash"]
