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
:: Batch script to invoke netcdf2aseg_gdf_converter Python script in MS-Windows
:: Written by Alex Ip 27/6/2018
:: Example invocation: aseg2nc <netcdf_input_path> <aseg_gdf_dat_output_path>

python -m geophys_utils.nc2aseg %*
