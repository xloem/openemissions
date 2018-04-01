#!/usr/bin/env python3

import argparse, os, sys, time

import numpy

from soapypower.writer import SoapyPowerBinFormat

parser = argparse.ArgumentParser(description='''
Watches specialized soapy power bin format files in the presence of square wave noise to monitor the noise level.

The format used is that of alexnask's recording settings branch with xloem's magnitude patch. TODO: provide links

To monitor the shielding effectiveness of an enclosure, place the noise generator some distance from the recorder
without the enclosure in place.  Record the noise dB level.  Then place the generator and recorder in the same
positions, but with the enclosure surrounding one of them.  Record the new noise dB level, and compare to the old.
''')

parser.add_argument('-i', '--file', default=[], type=argparse.FileType('rb'), nargs='*', help='soapy_power_bin files to monitor')
parser.add_argument('-d', '--dir', default=[], nargs='*', help='directories to watch for files')
parser.add_argument('-f', '--freq', default=[40], nargs='*', help='toggle frequencies of noise sources to watch')

args = parser.parse_args()

class WatchedDir:
    def __init__(self, path):
        self.path = path
        self.files = {}

    def more(self):
        ret = []
        for entry in os.scandir(self.path):
            if entry.is_file and entry.path not in self.files:
                try:
                    f = open(entry.path, 'rb')
                    try:
                        wf = WatchedFile(f)
                        self.files[entry.path] = wf
                        ret.append(wf)
                    except ValueError:
                        f.close()
                except FileNotFoundError:
                    pass
        return ret


class WatchedFile:
    spbfmt = SoapyPowerBinFormat()
    minsecs = 60 * 60

    def __init__(self, f):
        self.file = f
        self.lastSize = os.stat(self.file.fileno()).st_size
        self.header = WatchedFile.spbfmt.read_header(self.file)
        if self.header is None:
            raise ValueError('invalid file format')

    def grown(self):
        size = os.stat(self.file.fileno()).st_size
        if size > self.lastSize:
            self.lastSize = size
            return True
        else:
            return False

    def spectra(self):
        results = []
        result = WatchedFile.spbfmt.read(self.file)
        while result is not None:
            secs = result[0].time_stop - result[0].time_start
            if secs < WatchedFile.minsecs:
                WatchedFile.minsecs = secs
            results.append(result)
            result = WatchedFile.spbfmt.read(self.file)
        return results


class NoiseSource:
    def __init__(self, freq):
        self.freq = float(freq)

    def process(self, fil, header, array):
        # DC is at center array, equal to half length
        dc = int(len(array) / 2)

        # we care about magnitude every freq bins above that (the harmonics)
        num_harmonics = int((dc - 1) * header.step / self.freq)

        harmonics = [round(dc + i * self.freq / header.step) for i in range(1, num_harmonics + 1)]
        harmonics = [array[h] for h in harmonics]

        if not fil.header[0]['sweep'].log_scale:
            harmonics = 10 * numpy.log10(harmonics)

        # TODO: do something smarter than averaging the peak height
        avgDb = numpy.sum(harmonics) / num_harmonics

        tuned_freq = (header.stop - header.start) / 2 + header.start

        print('{} {} Hz {} dB {} MHz {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(header.time_start)), self.freq, avgDb, tuned_freq / 1000000, fil.header[1]))

dirs = [WatchedDir(d) for d in args.dir]

files = [WatchedFile(f) for f in args.file]
for d in dirs:
    print('Scanning {} ...'.format(d.path), end='\r')
    files.extend(d.more())

sources = [NoiseSource(freq) for freq in args.freq]

for f in files:
    for spectrum in f.spectra():
        for source in sources:
            source.process(f, *spectrum)

while True:
    growth = False
    deadline = time.time() + WatchedFile.minsecs / 2
    for d in dirs:
        print('Scanning {} ...'.format(d.path), end='\r')
        files.extend(d.more())
    for f in files:
        if f.grown():
            growth = True
            for spectrum in f.spectra():
                for source in sources:
                    source.process(f, *spectrum)
    now = time.time()
    if not growth and now < deadline:
        print('sleeping until {} ... ...'.format(time.strftime('%H:%M:%S', time.gmtime(deadline))), end='\r')
        time.sleep(deadline - now)
