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

# Create a script to handle library setup
RUN echo '#!/bin/bash\n\
setup_riegl_lib() {\n\
    local lib_name=$1\n\
    local lib_version=$2\n\
    local lib_path="/opt/${lib_name}-${lib_version}"\n\
    if [ -d "$lib_path" ]; then\n\
        echo "${lib_name} found. Setting up environment variables."\n\
        case $lib_name in\n\
            "rdblib")\n\
                export rdb_DIR="${lib_path}/interface/cpp"\n\
                export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${lib_path}/library"\n\
                ;;\n\
            "rivlib")\n\
                export RiVLib_DIR="${lib_path}"\n\
                export RiVLib_INCLUDE_DIRS="${lib_path}/include"\n\
                ;;\n\
        esac\n\
    else\n\
        echo "${lib_name} not found. Skipping setup."\n\
    fi\n\
}\n\
setup_riegl_lib "rdblib" "2.4.1"\n\
setup_riegl_lib "rivlib" "2.6.0"\n\
exec "$@"' > /usr/local/bin/entrypoint.sh && chmod +x /usr/local/bin/entrypoint.sh

# Copy Riegl libraries if they exist
COPY --chown=root:root riegl_libs/rdblib-2.4.1-x86_64-linux /opt/rdblib-2.4.1
COPY --chown=root:root riegl_libs/rivlib-2_6_0-x86_64-linux-gcc11 /opt/rivlib-2.6.0

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

# Use the entrypoint script
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# Default command
CMD ["/bin/bash"]