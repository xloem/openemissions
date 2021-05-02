import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'antenny','nyansat','station'))

from config.config import ConfigRepository
from antenny import esp32_antenna_api_factory as api_factory

raspi1 = api_factory(ConfigRepository('raspi1.json'))

raspi1.pwm_calibration()
