FROM ubuntu:22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set LC_ALL and unset LANGUAGE
ENV LC_ALL=C
ENV LANGUAGE=

WORKDIR /src

COPY . .

# Update, install dependencies and cleanup
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    git \
    build-essential \
    cmake \
    clang-tidy \
    libeigen3-dev \
    libproj-dev \
    libgeotiff-dev \
    libgtest-dev \
    libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /deps

# Create a directory for building packages
RUN mkdir -p /deps && \
    cd /deps

# Clone, build and clean up LASzip (3.x supports LAS 1.4 via laszip_api.h C API)
RUN git clone https://github.com/LASzip/LASzip.git && \
    cd LASzip && \
    git checkout tags/3.4.4 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf LASzip

# Clone, build and clean up Qhull
RUN git clone https://github.com/qhull/qhull.git && \
    cd qhull && \
    git checkout tags/v7.3.2 && \
    cd build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf qhull

# Clone, build and clean up libnabo
RUN git clone https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    git checkout tags/1.0.7 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf libnabo

WORKDIR /src

# Build and install raycloudtools
RUN mkdir /build && \
    cd /build && \
    cmake /src \
        -DGeoTIFF_INCLUDE_DIR=/usr/include/geotiff \
        -DGeoTIFF_LIBRARY=/usr/lib/x86_64-linux-gnu/libgeotiff.so \
        -DPROJ_INCLUDE_DIR=/usr/include/proj \
        -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
        -DWITH_QHULL=ON \
        -DWITH_LAS=ON \
        -DDOUBLE_RAYS=ON \
        -DWITH_TIFF=ON \
        -DRAYCLOUD_BUILD_TESTS=ON \
        -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf raycloudtools

# Clone and build TreeTools
RUN git clone https://github.com/csiro-robotics/treetools.git && \
    cd treetools && \
    mkdir -p build && \
    cd build && \
    cmake .. \
        -DPROJ_INCLUDE_DIR=/usr/include/proj \
        -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
        -DWITH_TIFF=ON \
        -DTREE_BUILD_TESTS=ON \
        -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf treetools
    
# Update ldconfig
RUN ldconfig /usr/local/lib

# Clean up build directory
RUN rm -rf /build

# Set the working directory
WORKDIR /workspace

# Handle entrypoint for interactive and non-interactive
RUN echo '#!/bin/bash\n\
if [ $# -eq 0 ]; then\n\
    exec /bin/bash\n\
else\n\
    exec "$@"\n\
fi' > /entrypoint.sh && chmod +x /entrypoint.sh

# Set the entrypoint to script
ENTRYPOINT ["/entrypoint.sh"]
