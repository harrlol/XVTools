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

# ---- Install latest R (4.4.x) from CRAN ----
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl gnupg ca-certificates \
    && mkdir -p /etc/apt/keyrings \
    && curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | gpg --dearmor -o /etc/apt/keyrings/cran-archive-keyring.gpg \
    && echo "deb [signed-by=/etc/apt/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
        > /etc/apt/sources.list.d/cran.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev \
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

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) \
      install.packages("BiocManager", repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(version = "3.21")' \
    && Rscript -e 'BiocManager::install("infercnv", ask = FALSE, update = FALSE)'

# RUN Rscript -e 'install.packages(c("optparse","data.table","ggplot2","reshape2","fastcluster","foreach","doParallel","futile.logger","RColorBrewer","gplots","png","Cairo"), repos="https://cloud.r-project.org")'

# python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --upgrade pip \
 && if [ -f /tmp/requirements.txt ]; then pip install -r /tmp/requirements.txt; fi \
 || pip install anndata scanpy pandas numpy scipy h5py tqdm google-cloud-storage google-auth

WORKDIR /app
COPY . /app

RUN pip install -e .

CMD ["/bin/bash"]
