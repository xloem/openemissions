import numpy as np
import os

import SoapySDR
from SoapySDR import *

from util import simulation

# soapysdr is not on pip, it is a python binding to a c++ library

soapy1 = SoapySDR.Device_enumerate()[-1]

#class Simulated:
#    def __init__(self, name='receiver'):
#        self.name = name
#        self.sim = simulation(name)
#        self.tuner = Simulated.Tuner(self, name + '_tuner')
#    def enable(self):
#    def disable(self):
#        self.sim = None
#    def read(self):
#        return self.sim.read()
#    class Tuner:
#        def __init__(self, fn):
#            self.sim = simulation(name)
#            self.last = self.sim.read()



class Receiver:
    def __init__(self, soapy_device_data = soapy1, channels = [0], rate = 1024 * 1024 * 2):
        self.soapy = SoapySDR.Device(soapy_device_data)
        self.channels = channels
        for channel in channels:
            self.soapy.setSampleRate(SOAPY_SDR_RX, channel, rate)
            # self.soapy.setGain(SOAPY_SDR_RX, channel, gain)
            # self.soapy.setBandwidth(SOAPY_SDR_RX, channel, bw)
            ## self.soapy.setFrequency(SOAPY_SDR_RX, channel, frequency)
        self.rate = self.soapy.getSampleRate(SOAPY_SDR_RX, 0)
        self.stream = self.soapy.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CF32, channels) 
        mtu = self.soapy.getStreamMTU(self.stream)
        self.buf = np.array([[0]*mtu]*len(channels), np.complex64) # np.complex64 is 2 floats
        self.tuner = Tuner(self)
    def __del__(self):
        self.soapy.closeStream(self.stream)
    def enable(self):
        self.soapy.activateStream(self.stream)
    def disable(self):
        self.soapy.deactivateStream(self.stream)
    def read(self):
        status = self.soapy.readStream(self.stream, self.buf, len(self.buf[0]))
        while status.ret == SoapySDR.SOAPY_SDR_OVERFLOW:
            print('warning: data dropped from radio')
            status = self.soapy.readStream(self.stream, self.buf, len(self.buf[0]))
            
        if status.ret <= 0:
            raise Exception(SoapySDR.errToStr(status.ret))
        if status.flags & SoapySDR.SOAPY_SDR_END_ABRUPT != 0:
            raiseException('Stream terminated prematurely (overflow)')
        # status.timeNs is accurate if SOAPY_SDR_HAS_TIME
        return self.buf
        
class Tuner:
    def __init__(self, receiver):
        self.receiver = receiver
        self.min = []
        self.max = []
        self.last = []
        for channel in self.receiver.channels:
            ranges = self.receiver.soapy.getFrequencyRange(SOAPY_SDR_RX, channel)
            self.min.append(min(range.minimum() for range in ranges))
            self.max.append(max(range.maximum() for range in ranges))
            self.last.append(self.receiver.soapy.getFrequency(SOAPY_SDR_RX, channel))
        self.min = np.array(self.min)
        self.max = np.array(self.max)
        self.last = np.array(self.last)
    def write_event(self, freqvec):
        for channel, frequency in zip(self.receiver.channels, freqvec):
            self.receiver.soapy.setFrequency(SOAPY_SDR_RX, channel, frequency)
