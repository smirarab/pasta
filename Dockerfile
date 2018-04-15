FROM debian:latest
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apt-get update && apt-get -y upgrade

# set up PASTA
RUN apt-get install -y python3 python3-setuptools openjdk-8-jre git
RUN cd /usr/local/bin
RUN git clone https://github.com/niemasd/pasta.git
RUN git clone https://github.com/smirarab/sate-tools-linux.git
RUN cd pasta && python3 setup.py develop
