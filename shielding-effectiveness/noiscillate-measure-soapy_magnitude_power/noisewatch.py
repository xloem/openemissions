#!/usr/bin/env python

import argparse, os, sys, time

import numpy as np

from soapypower.writer import SoapyPowerBinFormat

parser = argparse.ArgumentParser(description='''
Watches specialized soapy power bin format files in the presence of square wave noise to monitor the noise level.

The format used is that of alexnask's recording settings branch with xloem's magnitude patch. TODO: provide links

To monitor the shielding effectiveness of an enclosure, place the noise generator some distance from the recorder
without the enclosure in place.  Record the noise dB level.  Then place the generator and recorder in the same
positions, but with the enclosure surrounding one of them.  Record the new noise dB level, and compare to the old.
''')

parser.add_argument('-i', '--file', default=[], nargs='*', help='soapy_power_bin files to monitor')
parser.add_argument('-d', '--dir', default=[], nargs='*', help='directories to watch for files')
parser.add_argument('-f', '--freq', default=[40], nargs='*', help='toggle frequencies of noise sources to watch')
parser.add_argument('-e', '--expire', default=60, help='Number of minutes without data after which to close a file (default 60)')

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
                    wf = WatchedFile(entry.path)
                    self.files[entry.path] = wf
                    ret.append(wf)
                except ValueError:
                    pass
                except FileNotFoundError:
                    pass
        return ret


class WatchedFile:
    spbfmt = SoapyPowerBinFormat()
    minsecs = 60 * 60

    def __init__(self, path):
        self.file = open(path, 'rb')
        self.lastSize = os.stat(self.file.name).st_size
        self.lastTime = os.stat(self.file.name).st_mtime
        self.header = WatchedFile.spbfmt.read_header(self.file)
        if self.header is None:
            self.file.close()
            raise ValueError('invalid file format')

    def __del__(self):
        if hasattr(self, 'file'):
            self.file.close()

    def grown(self):
        size = os.stat(self.file.name).st_size
        if size > self.lastSize:
            self.lastSize = size
            return True
        else:
            return False

    def expired(self):
        global args
        return (time.time() - self.lastTime) > args.expire * 60

    def spectra(self):
        results = []
        #pos = self.file.tell()
        #sys.stdout.write('@{}\n'.format(pos))
        result = WatchedFile.spbfmt.read(self.file)
        while result is not None:
            self.lastTime = result[0].time_stop
            secs = result[0].time_stop - result[0].time_start
            if secs < WatchedFile.minsecs:
                WatchedFile.minsecs = secs
            results.append(result)
            #pos = self.file.tell()
            #sys.stdout.write('@{}\n'.format(pos))
            result = WatchedFile.spbfmt.read(self.file)

        
        return np.array(results)

def ExactPeakAverage(freq, array, step, log):
    # DC is at center array, equal to half length
    dc = int(len(array) / 2)

    # we care about magnitude every freq bins above that (the harmonics)
    num_harmonics = int((dc - 1) * step / freq)

    harmonics = [round(dc + i * freq / step) for i in range(1, num_harmonics + 1)]
    harmonics = [array[h] for h in harmonics]

    if not log:
        harmonics = 10 * np.log10(harmonics)

    return np.sum(harmonics) / num_harmonics

def AccumulateCurve(freq, data, step, log):
	periodSamp = freq / step
	periodBoxExtent = int(periodSamp / 2)

	periodArray = np.array([])

	# if true, will continue finding peaks through entire spectrum
	KEEPGOING = False
	
	CUTOFF = 100

	# remove data below DC
	data = data[int(len(data) / 2):]

	# convert to arithmetic domain
	if log:
		data = 10 ** (data / 10)
	#data = np.sqrt(data)

	numpeaks = int(len(data) / periodSamp - 0.5)

	# to accumulate blocks, we'll want to identify the mean peak height
	# and also the standard deviation of the noise contributing to it
	# we assume blocks are valid when the peak height is at least some multiple
	# of this standard distribution

	peakPeriods = []
	periodSum = np.zeros(periodBoxExtent * 2)
	periodCount = 0

	for peakid in range(numpeaks):
		peakSamp = int(round((peakid + 1) * periodSamp))
		peakPeriod = data[peakSamp - periodBoxExtent:peakSamp + periodBoxExtent]

		periodN = periodArray.shape[0]
		periodArray.resize((periodN + 1, peakPeriod.size), refcheck = False)
		periodArray[periodN] = peakPeriod

		if periodN == 0:
			continue

		peaks = periodArray[:,np.array([-1,0,1])+periodBoxExtent].mean(axis=1)
		antipeaks = periodArray[:,[-1,0,1]].mean(axis=1)
		
		noisefloor = antipeaks.mean()
		stddev = np.sqrt(np.sum((antipeaks - noisefloor) ** 2) / (antipeaks.size - 1))
		
		peakheight = peaks.mean() - noisefloor

		# TODO: we may be concerned with regard to a minimum or maximum peak height.
		#  in this case we may want to use 3*stddev to know whether peak is above or below
		#  whatever the cutoff is, given the recorded noise

		if periodN >= CUTOFF:
			if not KEEPGOING:
				break

		if peakheight > stddev * 2 or periodN >= CUTOFF:
			# 95% chance peak is not noise

			hz2 = (peakid + 1) * periodSamp * step
			hz1 = (peakid - periodN + 1) * periodSamp * step

			#print('Peaks {} - {}: nf={} sd={} height={} ratio={}'.format(peakid - periodN, peakid, noisefloor, stddev, peakheight, peakheight / stddev))
			#print('{} Hz: {}'.format((hz1+hz2)/2, peakheight))
			peakPeriods.append((hz1, hz2, periodArray))
			periodSum += np.sum(periodArray, axis=0)
			periodCount += periodN + 1
			periodArray = np.array([])

	if periodCount == 0:
		return float('-inf')

	periodSum /= periodCount
	noisefloor = periodSum[[-2,-1,0,1,2]].mean()
	strength = (periodSum - noisefloor).mean()

	power = strength #** 2
	dB = 10 * np.log10(power)

	return dB
 

class NoiseSource:
    def __init__(self, freq):
        self.freq = float(freq)

    def process(self, fil, header, array):

        tuned_freq = (header.stop - header.start) / 2 + header.start
        #print('tuned freq: {}'.format(tuned_freq))

        #dB = ExactPeakAverage(self.freq, array, header.step, fil.header[0]['sweep'].log_scale)
        dB = AccumulateCurve(self.freq, array, header.step, fil.header[0]['sweep'].log_scale)

        sys.stdout.write('{} {} Hz {} dB {} MHz {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(header.time_start)), self.freq, dB, tuned_freq / 1000000, fil.header[1]))

dirs = [WatchedDir(d) for d in args.dir]

files = set([WatchedFile(f) for f in args.file])
for d in dirs:
    sys.stdout.write('Scanning {} ...\r'.format(d.path))
    for f in d.more():
        files.add(f)

sources = [NoiseSource(freq) for freq in args.freq]

for f in files:
    for spectrum in f.spectra():
        for source in sources:
            source.process(f, *spectrum)

expiredfiles = []

while True:
    growth = False
    deadline = time.time() + WatchedFile.minsecs / 2
    for d in dirs:
        sys.stdout.write('Scanning {} ...\r'.format(d.path))
        for f in d.more():
            files.add(f)
    for f in files:
        if f.grown():
            growth = True
            for spectrum in f.spectra():
                for source in sources:
                    source.process(f, *spectrum)
        if f.expired():
            expiredfiles.append(f)
    for f in expiredfiles:
        files.remove(f)
        del f
    expiredfiles = []
    now = time.time()
    if not growth and now < deadline:
        sys.stdout.write('sleeping until {} ... ...\r'.format(time.strftime('%H:%M:%S', time.gmtime(deadline))))
        time.sleep(deadline - now)
