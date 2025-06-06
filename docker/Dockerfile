FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL image.author.name "M.Laurent"
LABEL image.author.email "maxime.laurent@uclouvain.be"

# Install dependencies, python and R
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    wget \
    unzip \
    libssl-dev \
    libcurl4-openssl-dev \
    unixodbc-dev \
    libxml2-dev \
    libmariadb-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \ 
    libjpeg-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install plink 1.9
RUN wget -c http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip -d /usr/local/bin/ && \
    rm plink_linux_x86_64_20231211.zip && \
    chmod +x /usr/local/bin/plink

# Install FastLMM via pip
RUN pip3 install fastlmm
RUN chmod +x /usr/local/lib/python3.10/dist-packages/fastlmm/association/Fastlmm_autoselect/fastlmmc

RUN R -e 'install.packages("xml2", dependencies = T, INSTALL_opts = c("--no-lock"))'
RUN R -e 'install.packages("tidyverse")'

# Add FastLMM directory to the PATH
ENV PATH="/usr/local/lib/python3.10/dist-packages/fastlmm/association/Fastlmm_autoselect/:${PATH}"