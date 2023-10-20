# This Dockerfile exists only for Binder deployment
# It uses a base image created using repo2docker in BinderHub to speed up the build process
FROM alexip/conda-geophys_utils_base:0.1.0_gdal3.7.0_mcpy3.10_r2d

USER root

# Take a complete copy of the project directory into /geophys_utils (could do a git pull)
COPY . /home/jovyan/geophys_utils

# install geophys_utils and set up scripts
RUN chown -R jovyan:jovyan /home/jovyan/geophys_utils && \
    chmod -R u+rwX /home/jovyan/geophys_utils && \
    /srv/conda/envs/notebook/bin/pip install -e /home/jovyan/geophys_utils && \
    for script in $(ls /home/jovyan/geophys_utils/bin/* | grep -v '\.bat$'); \
        do chmod 755 "${script}"; \
        ln -s "/home/jovyan/geophys_utils/bin/${script}" /usr/local/sbin/$(basename "${script}"); \
    done

USER jovyan