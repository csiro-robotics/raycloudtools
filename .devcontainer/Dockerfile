FROM ubuntu:latest

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set LC_ALL and unset LANGUAGE
ENV LC_ALL=C
ENV LANGUAGE=

# Update, install dependencies and cleanup
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    git \
    build-essential \
    gdb \
    cmake \
    cmake-curses-gui \
    clang-tidy \
    libeigen3-dev \
    libproj-dev \
    libgeotiff-dev \
    libboost-all-dev \
    wget \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    libgtest-dev && \
    rm -rf /var/lib/apt/lists/*

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
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf libnabo

# Update ldconfig
RUN ldconfig /usr/local/lib

# Copy build script into image
COPY .devcontainer/commands/compile_install.sh /usr/local/bin/compile_install.sh
RUN chmod +x /usr/local/bin/compile_install.sh

# Set the working directory
WORKDIR /workspaces/raycloudtools

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]