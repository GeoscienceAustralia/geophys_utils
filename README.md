# geophys\_utils - Utilities for Accessing netCDF Encoded Geophysics Data
The geophys_utils Python library is a collection of utilities for discovering, accessing and using geophysics data via 
web services or from netCDF files.

Details on the netCDF encoding for point and line data can be viewed at 
<https://docs.google.com/document/d/1C-SsT1vOcAaPT_4jdY1S_NjUjbk-WbUjb1FCw7uPrxw/edit?usp=sharing>

Feedback, docker/Dockerfile_0.1.0_gdal3.7.0_py3.10.6_ubuntu22.04.2suggestions or contributions will be gratefully 
accepted.

# Docker deployment
There is a Dockerfile in the project root directory which will allow you to create a fresh Docker image from the 
current source code using the command
```docker build --progress=plain -t <image_name>:<image_tag> -f docker/Dockerfile_<version> .``` 
in the project root directory. 
Note that, depending on the OS, Python, and GDAL version, some builds can take longer due to the extensive 
dependencies for large packages such as NetCDF and SciPy.

Alternatively, you can pull a pre-made Ubuntu image with the latest versins using the following commands:
```
docker pull alexip/geophys_utils:latest
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
