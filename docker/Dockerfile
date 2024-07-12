FROM frolvlad/alpine-glibc:alpine-3

RUN apk update && apk add --no-cache make cmake g++ git proj-dev eigen-dev boost-dev libgeotiff-dev gtest-dev wget

# Adding FindGeoTIFF.cmake
RUN cd /usr/share/cmake/Modules && \
    wget https://raw.githubusercontent.com/ufz/geotiff/master/cmake/FindGeoTIFF.cmake

# Building libnabo
RUN cd /home && \
    git clone https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    git checkout tags/1.0.7 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

# Building LASzip
RUN cd /home && \
    git clone https://github.com/LASzip/LASzip.git && \
    cd LASzip && \
    git checkout tags/2.0.1 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

# Building libLAS
RUN cd /home && \
    git clone https://github.com/libLAS/libLAS.git && \
    cd libLAS && \
    mkdir build && cd build && \
    cmake .. -DWITH_LASZIP=ON -DWITH_GEOTIFF=OFF -DCMAKE_CXX_STANDARD=11 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

# Building Qhull
RUN cd /home && \
    git clone http://github.com/qhull/qhull.git && \
    cd qhull && \
    git checkout tags/v7.3.2 && \
    mkdir -p build && cd build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

# Building raycloudtools
RUN cd /home && \
    git clone https://github.com/csiro-robotics/raycloudtools.git && \
    cd raycloudtools && \
    mkdir build && cd build && \
    cmake .. -DGEOTIFF_LIBRARIES=/usr/lib/libgeotiff.so -DWITH_QHULL=ON -DWITH_LAS=ON -DWITH_TIFF=ON -DRAYCLOUD_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

# Building treetools
RUN cd /home && \
    git clone https://github.com/csiro-robotics/treetools.git && \
    cd treetools && \
    mkdir build && cd build && \
    cmake .. -DWITH_TIFF=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3" && \
    make -j$(nproc) && make install

WORKDIR /data
