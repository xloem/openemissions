import receiver
import station
import numpy as np

s = station.Station(station.raspi1_noimu)

s.write_event(s.min + s.max * np.random.random(s.max.size))

r = receiver.Receiver(receiver.soapy1)
r.enable()
print(r.read())

