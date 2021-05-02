import SoapySDR

# soapysdr is not on pip, it is a python binding to a c++ library

soapy1 = SoapySDR.Device_enumerate()[-1]
