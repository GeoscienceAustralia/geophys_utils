# geophys\_utils - Utilities for Accessing netCDF Encoded Geophysics Data
The geophys_utils Python library is a collection of utilities for discovering, accessing and using geophysics data via 
web services or from netCDF files. This library is maintained by Geoscience Australia with ongoing contributions from
the original Author, Alex Ip.

An overview of the netCDF encodings for grid, point and line geophysical data can be viewed at 
<https://doi.org/10.1080/22020586.2019.12073191>

Feedback, suggestions or contributions will be gratefully accepted.

# Docker deployment
There is a docker directory in the project root which contains Dockerfiles which will allow you to create fresh Docker 
images from the current source code using the command
```docker build --progress=plain -t <image_name>:<image_tag> -f docker/Dockerfile_<version> .``` 
in the project root directory. Alternatively, you can run the ```build_all_images.sh``` script.
Note that, depending on the OS, Python, and GDAL version, some builds can initially take longer due to the extensive 
dependencies for large packages such as NetCDF and SciPy.

Alternatively, you can pull a pre-made Ubuntu image with the latest versions using the following commands:
```bash
docker pull alexip/geophys_utils:latest
# Run unit tests
docker run -it --name gu-test alexip/geophys_utils:latest python -m geophys_utils.test
# Start a shell in the container
docker run -it --name gu-test alexip/geophys_utils:latest /bin/bash
```

## License
The content of this repository is licensed for use under the 
[Apache 2.0 License](http://www.apache.org/licenses/LICENSE-2.0). 
See the [license deed](https://github.com/GeoscienceAustralia/geophys_utils/blob/master/LICENSE) for full details.

## Contacts
**Geoscience Australia**
*Publisher*
<https://www.ga.gov.au>
<clientservices@ga.gov.au>

**Alex Ip**
*Principal Author*
<mailto://alex@trentham.net.au>
<https://orcid.org/0000-0001-8937-8904>

**Andrew Turner**
*Contributing Author*
<https://orcid.org/0000-0001-5085-8783>
