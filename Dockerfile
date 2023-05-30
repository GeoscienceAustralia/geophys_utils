# alexip/geophys_utils:0.1.0 contains geophys_utils v0.1.0 and all dependencies
# GDAL Docker image sourced from https://github.com/OSGeo/gdal/

FROM ghcr.io/osgeo/gdal:alpine-normal-3.5.3

ENV PYTHONUNBUFFERED=1

# System dependencies for geophys_utils 0.1.0 / GDAL 3.5.3 / Python 3.9
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
    'py3-scipy' \
    'py3-pip' \
    'python3-dev' \
    'hdf5' \
    'hdf5-dev' \
    'netcdf' \
    'netcdf-dev' \
    'netcdf-utils' \
    'geos' \
    'geos-dev'

# Take a complete copy of the project directory into /geophys_utils
RUN mkdir /geophys_utils
WORKDIR /geophys_utils
COPY .. /geophys_utils

# Install geophys_utils package from disk
# need to force scikit-image==0.19.3 due to broken scipy==1.9.1 dependency
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --upgrade setuptools wheel && \
    pip install scikit-image==0.19.3  && \
    pip install --no-cache-dir -e .
