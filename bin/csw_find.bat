@ECHO OFF
:: Batch file to invoke _csw_utils Python script in MS-Windows
:: Written by Written by Alex Ip & Andrew Turner 27/2/2017
:: Example invocation: csw_find -k "NCI, AU, grid, potassium" -b 148.996,-35.48,149.399,-35.124

python -m geophys_utils._csw_utils  %*
