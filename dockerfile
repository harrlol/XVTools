FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# base R and python
RUN apt-get update && apt-get install -y --no-install-recommends \
    tar netbase ca-certificates coreutils \
    build-essential wget curl git  tzdata locales \
    gfortran cmake liblapack-dev libblas-dev \
    python3 python3-pip python3-venv \
    r-base libcurl4-openssl-dev libssl-dev libxml2-dev \
    libhdf5-dev zlib1g-dev libbz2-dev liblzma-dev \
    libpng-dev libcairo2-dev \
    libzstd-dev \
 && rm -rf /var/lib/apt/lists/*

# fetch, build, install
RUN wget https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.2.tar.gz -O JAGS-4.3.2.tar.gz \
    && tar -xzf JAGS-4.3.2.tar.gz \
    && cd JAGS-4.3.2 \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf JAGS-4.3.2*

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")' \
 && Rscript -e 'BiocManager::install("infercnv", ask=FALSE, update=FALSE)'

# RUN Rscript -e 'install.packages(c("optparse","data.table","ggplot2","reshape2","fastcluster","foreach","doParallel","futile.logger","RColorBrewer","gplots","png","Cairo"), repos="https://cloud.r-project.org")'

# python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --upgrade pip \
 && if [ -f /tmp/requirements.txt ]; then pip install -r /tmp/requirements.txt; fi \
 || pip install anndata scanpy pandas numpy scipy h5py tqdm google-cloud-storage google-auth


WORKDIR /app
COPY . /app
ENV R_HELPER=/app/R/src/infercnv_helper.R

CMD ["/bin/bash"]
