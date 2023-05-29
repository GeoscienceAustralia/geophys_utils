# alexip/geophys_utils:0.1.0 contains geophys_utils v0.1.0 and all dependencies
# GDAL Docker image sourced from https://github.com/OSGeo/gdal/

FROM ghcr.io/osgeo/gdal:alpine-normal-3.6.4

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# System dependencies for geophys_utils 0.1.0 / GDAL 3.6.4 / Python 3.10.11-r0
RUN apk add --update --no-cache \
    'build-base' \
    'gcc' \
    'g++' \
    'gfortran' \
    'make' \
    'cmake' \
    'musl-dev' \
    'python3-dev' \
    'git' \
    'openblas-dev' \
    'lapack-dev' \
    'py3-pip' \
    'python3-dev=3.10.11-r0' \
    'hdf5=1.12.2-r0' \
    'hdf5-dev=1.12.2-r0' \
    'netcdf=4.8.1-r2' \
    'netcdf-dev=4.8.1-r2' \
    'geos=3.10.3-r0' \
    'geos-dev=3.10.3-r0'

# Take a complete copy of the project directory into /geophys_utils
RUN mkdir /geophys_utils
WORKDIR /geophys_utils
COPY . /geophys_utils

# Install geophys_utils package from disk
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --upgrade setuptools wheel && \
    pip install --no-cache-dir -e .
