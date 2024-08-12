FROM ubuntu:22.04

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
    cmake \
    cmake-curses-gui \
    clang-tidy \
    libeigen3-dev \
    libproj-dev \
    libgeotiff-dev \
    libboost-all-dev \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Copy the conditional copy script
COPY .devcontainer/copy_external_libs.sh /usr/local/bin/copy_external_libs.sh
RUN chmod +x /usr/local/bin/copy_external_libs.sh

# Run the script to handle copying RiVLib and rdblib if they exist
RUN /usr/local/bin/copy_external_libs.sh

# Clone, build and clean up LASzip
RUN git clone https://github.com/LASzip/LASzip.git && \
    cd LASzip && \
    git checkout tags/2.0.1 && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cp bin/Release/liblas* /usr/lib/ && \
    cd ../.. && rm -rf LASzip

# Clone, build and clean up libLAS
RUN git clone https://github.com/libLAS/libLAS.git && \
    cd libLAS && \
    mkdir build && cd build && \
    cmake .. -DWITH_LASZIP=ON -DWITH_GEOTIFF=OFF -DCMAKE_CXX_STANDARD=11 && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && rm -rf libLAS

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

# Build RayCloudTools 
RUN mkdir build && \
    cd build && \
    cmake .. \
        -DGeoTIFF_INCLUDE_DIR=/usr/include/geotiff \
        -DGeoTIFF_LIBRARY=/usr/lib/x86_64-linux-gnu/libgeotiff.so \
        -DPROJ_INCLUDE_DIR=/usr/include/proj \
        -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
        -DWITH_LAS=ON \
        -DWITH_RIEGL=ON \
        -DDOUBLE_RAYS=ON \
        -DWITH_TIFF=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make -j"$(nproc)" && \
    make install && \
    cd ../.. && \
    ldconfig /usr/local/lib

# Clone and build TreeTools
RUN git clone https://github.com/tim-devereux/treetools.git  && \
    cd treetools  && \
    mkdir build  && \
    cd build  && \
    cmake ..  \
        -DPROJ_INCLUDE_DIR=/usr/include/proj \
        -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
        -DWITH_TIFF=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    make -j"$(nproc)"  && \
    make install  && \
    cd ../..

# Update ldconfig
RUN ldconfig /usr/local/lib

# Set the working directory
WORKDIR /workspaces/raycloudtools

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]
