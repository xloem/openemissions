import os
import numpy as np
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'antenny','nyansat','station'))

from config.config import ConfigRepository
from antenny import esp32_antenna_api_factory as api_factory

raspi1 = 'antenny_raspi1.json'
raspi1_noimu = 'antenny_raspi1_noimu.json'

class Station:
    def __init__(self, antenny_configfile = raspi1):
        config = ConfigRepository(antenny_configfile)
        while True:
            try:
                self.antenny = api_factory(config)
                break
            except ValueError as e:
            #    if 'I2C' in e.args[0]:
            #        if '0x28' in e.args[0] or '0x29' in e.args[0]:
            #            print('Warning: bno055 imu not found.', e.args[0])
            #            newconfig = antenny_configfile + '_nobno055.json'
            #            print('Reconfiguring without imu as', newconfig)
            #            config.new(newconfig)
            #            config.set('use_imu', False)
            #            continue
                raise e
        if self.antenny.is_safemode():
            raise Exception('Failed to initialise motor driver')

        self.controllers = [self.antenny.antenna.elevation, self.antenny.antenna.azimuth]

        self.min = []
        self.max = []
        self.last = []
        for controller in self.controllers:
            self.min.append(controller.motor.min_duty)
            self.max.append(controller.motor.max_duty)
            self.last.append(controller.get_duty())
        self.min = np.array(self.min)
        self.max = np.array(self.max)
        self.last = np.array(self.last)
    def write_event(self, dutyvec):
        self.last = dutyvec
        self.antenny.antenna.elevation.set_duty(dutyvec[0])
        self.antenny.antenna.azimuth.set_duty(dutyvec[1])











    # it would be good to abstract the concept of 'something needing calibration' out
    # this could be 'motor angles', and a standard interface could enumerate it
    # this calibration needs a set of recordings that cover the full range of motion of the motors
    # it has a region it could narrow down for more information

    # both motors have a 180 degree maximum range
    # so calibration can be found by looking for pairs of values where elevation2 = 90-elevation1
    # and azimuth2 = azimuth1 + 180
    # such pairs will aim the antenna in the same direction.  if the antenna is symmetrical the
    # signal would be approximately the same.

    # consider data with pairs of duty information
    # [elevation_duty, azimuth_duty, spectrum]
    # we can cast a formula
    # duty = angle * constant + constant
    # both this and the pair differences are linear transforms, so we can cast the same for duty
    # elevation_duty2 = elevation_duty_vertical - elevation_duty1
    # azimuth_duty2 = 180 * azimuth_duty_per_degree + azimuth_duty1
    # we then basically autocorrelate the spectrums to solve for the two unknowns

    
    

#raspi1.pwm_calibration()
#  it looks like pwm_calibration assumes azimuth is euler angle 0, and elevation euler angle 2
#   the imu can be separately calibrated to itself first, it has internal memory for this
#  pwm_calibration is supposed to return duty-per-degree.  it doesn't look correct to me.
#  it looks to me like duty-per-degree should be 3-4 for both axes (even though they are differnt
#   servo models).
#  for quick use, it makes sense to assume 180-degree extent for either servo
