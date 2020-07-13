FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      "file" \
      "build-essential" \
      "libboost-all-dev" \
      "liblemon-dev" \
      "libopencv-highgui-dev"\
      "libpng-dev" \
      "libtiff5-dev" \
      "unzip" \
      "wget" \
      "mc"\
      "libopencv-dev"
      
##
# Download sources
##

RUN wget --no-check-certificate "https://cmake.org/files/v3.15/cmake-3.15.5.tar.gz" -O /cmake-3.15.5.tar.gz
RUN tar -xvf /cmake-3.15.5.tar.gz

#ADD ./json-3.7.3.zip /json-3.7.3.zip
RUN wget --no-check-certificate "https://github.com/nlohmann/json/archive/v3.7.3.zip" -O /json-3.7.3.zip
RUN unzip /json-3.7.3.zip

##
# Build CMake
##

RUN cd cmake-3.15.5 && ./bootstrap && make -j 4 && make install

##
# Build NLohmann-JSON
##

RUN mkdir json-3.7.3/build
RUN cd json-3.7.3/build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .. && make -j 4 && make install