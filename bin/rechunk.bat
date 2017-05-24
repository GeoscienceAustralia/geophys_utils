@ECHO OFF
::===============================================================================
::    Copyright 2017 Geoscience Australia
:: 
::    Licensed under the Apache License, Version 2.0 (the "License");
::    you may not use this file except in compliance with the License.
::    You may obtain a copy of the License at
:: 
::        http://www.apache.org/licenses/LICENSE-2.0
:: 
::    Unless required by applicable law or agreed to in writing, software
::    distributed under the License is distributed on an "AS IS" BASIS,
::    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
::    See the License for the specific language governing permissions and
::    limitations under the License.
::===============================================================================
:: Batch file to invoke _netcdf_utils Python script to re-chunk netCDF file in MS-Windows
:: Written by Written by Alex Ip 2/3/2017
:: Example invocation: rechunk infile.nc outfile.nc --chunking=8192,8192

python -m geophys_utils._netcdf_utils --copy  %*
