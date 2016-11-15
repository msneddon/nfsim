FROM ubuntu:latest
MAINTAINER Michael Sneddon "michaelsneddon@gmail.com"

RUN apt-get update && \
    apt-get -y install g++ make

COPY ./ /nfsim

RUN cd /nfsim && make

RUN mkdir /work
WORKDIR /work

ENTRYPOINT ["/nfsim/bin/NFsim"]
