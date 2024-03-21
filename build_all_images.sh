#!/bin/bash
# build_all_images.sh - A script to build all Docker images defined in the Dockerfiles in the docker subdirectory
# Run this script from the project root (i.e. the parent directory of the "docker" subdirectory)
for dockerfile in $(ls docker/Dockerfile_*)
do
  tag="$(echo "${dockerfile}" | cut -d_ -f2-)"
  image="$(basename "$(pwd)"):$tag"
  docker build --progress=plain -t "${image}" -f "${dockerfile}" .
done
