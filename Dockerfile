FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install -y eatmydata \
    && eatmydata apt-get install -y build-essential wget bzip2 \
      ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 libz-dev \
      git unzip python3.10 pip \
    && apt-get clean

# make sure we can just issue "python"
RUN ln -s /usr/bin/python3.10 /usr/bin/python

# bowtie
RUN wget -q -O bowtie.zip https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip; \
	unzip bowtie.zip -d /opt/; \
	ln -s /opt/bowtie-1.3.1-linux-x86_64 /opt/bowtie; \
	rm bowtie.zip
ENV PATH $PATH:/opt/bowtie

# copy the source
COPY . /src/gbrs

RUN pip install -U pip; cd /src/gbrs; pip install .; cd
