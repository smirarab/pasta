FROM debian:latest
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apt-get update && apt-get -y upgrade

# set up PASTA
RUN apt-get install -y python3 python3-setuptools default-jre git libgomp1
RUN ln -s $(which python3) /usr/local/bin/python
RUN cd /usr/local/bin
RUN git clone https://github.com/smirarab/pasta.git
RUN git clone https://github.com/smirarab/sate-tools-linux.git
# RUN cd sate-tools-linux && git checkout d78ef029b533e4f4ac13ba1e9cfdac7944b1f70e && cd ..
RUN cd pasta && python3 setup.py develop
RUN rm /usr/local/bin/run_pasta_gui.py
ENV CONTRALIGN_DIR /usr/local/bin/sate-tools-linux

RUN mkdir /data
WORKDIR /data
