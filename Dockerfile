FROM alpine:3

RUN apk uptate && \
    apk add make cmake g++ git proj-dev eigen-dev boost-dev gtest-dev && \
    cd /usr/share/cmake/Modules && \
    wget https://raw.githubusercontent.com/ufz/geotiff/master/cmake/FindGeoTIFF.cmake

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
    mkdir build && cd build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true && \
    make && make install

RUN cd /home && \
    git clone https://github.com/tiagodc/raycloudtools.git && \
    cd raycloudtools && \
    mkdir build && cd build && \
    cmake .. -DGEOTIFF_INCLUDE_DIR=/usr/include/geotiff -DGEOTIFF_LIBRARIES=/usr/lib/libgeotiff.so -DWITH_QHULL=ON -DWITH_LAS=ON -DWITH_TIFF=ON -DRAYCLOUD_BUILD_TESTS=ON && \
    make && make install

