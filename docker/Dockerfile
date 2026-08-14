FROM ubuntu:26.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set LC_ALL and unset LANGUAGE
ENV LC_ALL=C
ENV LANGUAGE=

# Update, install dependencies and cleanup
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
        build-essential \
        clang-tidy \
        cmake \
        git \
        kolourpaint \
        libboost-all-dev \
        libeigen3-dev \
        libgeotiff-dev \
        libgtest-dev \
        liblaszip-dev \
        libnetcdf-c++4-dev \
        libnetcdf-dev \
        libproj-dev \
        libqhull-dev \
        libtbb-dev \
        libtiff-dev \
        meshlab && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /deps

# Create a directory for building packages
RUN mkdir -p /deps && \
    cd /deps

# Clone, build and clean up libnabo
RUN git clone https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j "$(nproc)" && \
    make install && \
    cd ../.. && \
    rm -rf libnabo

WORKDIR /src

COPY . .

# Build and install raycloudtools
RUN mkdir /build && \
    cd /build && \
    cmake /src \
        -DCMAKE_BUILD_TYPE=Release \
        -DRAYCLOUD_BUILD_TESTS=ON \
        -DRAYLIB_DOUBLE_RAYS=ON \
        -DRAYLIB_WITH_LAS=ON \
        -DRAYLIB_WITH_QHULL=ON \
        -DRAYLIB_WITH_TIFF=ON && \
    make -j "$(nproc)" && \
    make install && \
    cd ../.. && \
    rm -rf raycloudtools

# Clone and build TreeTools
RUN git clone https://github.com/csiro-robotics/treetools.git && \
    cd treetools && \
    rm --force cmake/FindQhull.cmake cmake/FindTBB.cmake && \
    mkdir -p build && \
    cd build && \
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DTREE_BUILD_TESTS=ON \
        -DWITH_TIFF=ON && \
    make -j "$(nproc)" && \
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
