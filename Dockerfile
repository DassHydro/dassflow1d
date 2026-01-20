FROM ubuntu:22.04 as stage0
MAINTAINER Kevin Larnier "kevin.larnier@hydro-matters.fr"

## Stage 1
FROM stage0 as stage1
# Install packages
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
RUN mkdir -p /usr/share/man/man1
RUN apt-get update
RUN apt-get install -y build-essential gfortran g++ zlib1g-dev libblas-dev
RUN apt-get install -y python3-pip cython3
RUN apt search jre
RUN apt-get install -y openjdk-21-jre-headless

# Setup Python environment
RUN pip3 install numpy==1.26.4
RUN pip3 install f90wrap==0.2.13
RUN pip3 install geopandas
RUN pip3 install matplotlib
RUN pip3 install netcdf4
RUN pip3 install pandas
RUN pip3 install scipy
RUN pip3 install tqdm

## Stage 2
FROM stage1 as stage2
# Create directories
RUN mkdir -p /mnt/rundir
RUN mkdir -p /app/DassFlow1D_AQ
RUN mkdir -p /app/DassFlow1D_AQ/build
RUN mkdir -p /app/DassFlow1D_AQ/build/tap

# Copy DassFlow-1D sources
COPY ./src /app/DassFlow1D_AQ/src
COPY ./libs /app/DassFlow1D_AQ/libs
COPY ./tools /app/DassFlow1D_AQ/tools
COPY ./* /app/DassFlow1D_AQ/

# Copy and setup TAPENADE
ENV TAPENADE_HOME /app/DassFlow1D_AQ/libs/tapenade_3.16
ENV PATH $TAPENADE_HOME/bin:$PATH

# Compile DassFlow-1D and setup CLI
WORKDIR /app/DassFlow1D_AQ
RUN make ADJOINT=1 OPTIM=1 print_params
RUN make cleanlibs alllibs
RUN make ADJOINT=1 clean
RUN make ADJOINT=1 generate_adjoint
RUN make ADJOINT=1 OPTIM=1
COPY ./src/run_dassflow1d.py /app/DassFlow1D_AQ/bin/dassflow1d_cli.py

## Stage 3
FROM stage2 as stage3

# Set WORKDIR and ENV
WORKDIR /mnt/rundir
ENV PYTHONPATH /app/DassFlow1D_AQ/build/api

# Entry point
LABEL version="2.1"
LABEL description="DassFlow-1D v2.1"
LABEL maintainer="Kevin Larnier (kevin.larnier@hydro-matters.fr)"
ENTRYPOINT  ["python3", "/app/DassFlow1D_AQ/bin/dassflow1d_cli.py" ]
