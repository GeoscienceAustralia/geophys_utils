FROM python:3.9-alpine

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# System dependencies for geophys_utils 0.1.0 / GDAL 3.6.4 / Python 3.9
RUN apk add --update --no-cache \
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
    'hdf5-dev=1.14.0-r0' \
    'netcdf-dev=4.9.2-r2' \
    'gdal=3.6.4-r4' \
    'gdal-dev=3.6.4-r4' \
    'proj=9.2.0-r0' \
    'proj-dev=9.2.0-r0' \
    'proj-util=9.2.0-r0' \
    'proj-data=1.13-r0' \
    'geos=3.11.2-r0' \
    'geos-dev=3.11.2-r0'

# Take a complete copy of the project directory into /geophys_utils
RUN mkdir /geophys_utils
WORKDIR /geophys_utils
COPY . /geophys_utils

# Install geophys_utils from disk
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -e .

#RUN adduser -D user
#USER user
