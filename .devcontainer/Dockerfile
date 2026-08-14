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
        cmake-curses-gui \
        doxygen \
        gdb \
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
        meshlab \
        sudo \
        wget && \
    rm -rf /var/lib/apt/lists/*

# Clone, build and clean up libnabo
RUN git clone https://github.com/ethz-asl/libnabo.git && \
    cd libnabo && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make -j "$(nproc)" && \
    make install && \
    cd ../.. && \
    rm -rf libnabo

# Update ldconfig
RUN ldconfig /usr/local/lib

# Set the working directory
WORKDIR /workspaces/raycloudtools

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]

# Make sure ubuntu user exists and allow them to use sudo.
RUN if ! id -u ubuntu; then \
        useradd \
            --create-home \
            --groups=sudo \
            ubuntu; \
    fi && \
    echo "%sudo ALL=(ALL:ALL) NOPASSWD: ALL" >> /etc/sudoers
