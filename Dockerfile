FROM continuumio/miniconda

RUN apt-get -y update && \
    apt-get -y install build-essential tcsh && \
    # apt-get -y install blast2 && \
    apt-get -y clean && \
    pip install -U pip

RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# RUN conda install blast-legacy==2.2.22 hmmer==3.1b2
RUN conda install hmmer==3.1b2

ENV APPS /software/apps
RUN mkdir -p ${APPS}
WORKDIR ${APPS}
COPY . ${APPS}/inmembrane
WORKDIR ${APPS}/inmembrane
RUN pip install -U -r requirements.txt && pip install -e .

RUN mkdir /data
WORKDIR /data
ENTRYPOINT ["inmembrane_scan"]
# CMD ["-t"]
# CMD ["--help"]
