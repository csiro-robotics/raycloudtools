FROM ubuntu:22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set LC_ALL and unset LANGUAGE
ENV LC_ALL=C
ENV LANGUAGE=

# Update and install dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    git \
    build-essential \
    cmake \
    clang-tidy \
    libeigen3-dev \
    libproj-dev \
    libboost-all-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Copy Riegl libraries if they exist - If using -DWITH_RIEGL put them in a folder 'riegl_libs' in top directory before making container.
COPY --chown=root:root riegl_libs/rdblib-2.4.1-x86_64-linux /opt/rdblib-2.4.1
COPY --chown=root:root riegl_libs/rivlib-2_6_0-x86_64-linux-gcc11 /opt/rivlib-2.6.0

# Set environment variables for RiVLib and rdblib paths
ENV RiVLib_DIR=/opt/rivlib-2.6.0 \
    RiVLib_INCLUDE_DIRS=/opt/rivlib-2.6.0/include \
    rdb_DIR=/opt/rdblib-2.4.1/interface/cpp \
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/rdblib-2.4.1/library

# Clone and build LASzip
ENV LASZIP_VERSION=2.0.1
RUN git clone --branch ${LASZIP_VERSION} --depth 1 https://github.com/LASzip/LASzip.git && \
    cd LASzip && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    sudo make install && \
    sudo cp bin/Release/liblas* /usr/lib/ && \
    cd ../.. && rm -rf LASzip

# Clone and build libLAS
ENV LIBLAS_VERSION=master
RUN git clone --branch ${LIBLAS_VERSION} --depth 1 https://github.com/libLAS/libLAS.git && \
    cd libLAS && \
    mkdir build && \
    cd build && \
    cmake .. -DWITH_LASZIP=ON && \
    make -j$(nproc) && \
    sudo make install && \
    cd ../.. && rm -rf libLAS

# Install Qhull
ENV QHULL_VERSION=v7.3.2
RUN git clone --branch ${QHULL_VERSION} --depth 1 https://github.com/qhull/qhull.git && \
    cd qhull/build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf qhull

# Install libnabo
ENV LIBNABO_VERSION=1.0.7
RUN git clone --branch ${LIBNABO_VERSION} --depth 1 https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf libnabo

# Update ldconfig
RUN ldconfig /usr/local/lib

# Set the working directory
WORKDIR /workspace

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]