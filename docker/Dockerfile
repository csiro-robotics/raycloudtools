FROM alpine:3

RUN apk update

RUN apk add make cmake g++ git proj-dev eigen-dev boost-dev libgeotiff-dev gtest-dev

RUN cd /home && \
    git clone https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    git checkout tags/1.0.7 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make && make install

RUN cd /home && \
    git clone https://github.com/LASzip/LASzip.git && \
    cd LASzip && \
    git checkout tags/2.0.1 && \
    mkdir build && cd build && \ 
    cmake .. && \
    make && make install

RUN cd /home && \
    git clone https://github.com/libLAS/libLAS.git && \
    cd libLAS && \ 
    mkdir build && cd build && \
    cmake .. -DWITH_LASZIP=ON -DWITH_GEOTIFF=OFF && \
    make && make install

RUN cd /home && \
    git clone http://github.com/qhull/qhull.git && \
    cd qhull && \
    git checkout tags/v7.3.2 && \
    mkdir -p build && cd build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true && \
    make && make install

RUN cd /home && \
    git clone https://github.com/csiro-robotics/raycloudtools.git && \
    cd raycloudtools && \
    mkdir build && cd build && \
    cmake .. -DWITH_QHULL=ON -DWITH_LAS=ON -DRAYCLOUD_BUILD_TESTS=ON && \
    make && make install 
    
RUN cd /home && \
    git clone https://github.com/csiro-robotics/treetools.git && \
    cd treetools && \
    mkdir build && cd build && \
    cmake .. && \
    make && make install

WORKDIR /data
