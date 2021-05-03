import numpy as np

# wraps an ndarray providing for in-place addition of rows
class npflexarray:
    def __init__(self, array = None):
        self._array = None
        self._view = None
        self.len = 0
        if array is not None:
            self.extend(array)
    def extend(self, more):
        if self._array is None:
            self._array = np.array(more)
            self._view = self._array
            self.len = len(more)
            return
        newlen = self.len + len(more)
        if newlen >= len(self._array):
            # double size of array
            self._array = np.append(self._array, self._array, axis=0)
        self._array[self.len:newlen] = more
        self._view = self._array[:newlen]
        self.len = newlen
    def array(self):
        return self._view
    def __getattr__(self, name):
        return getattr(self._view, name)
    def __getitem__(self, idx):
        return self._view[idx]
    def __len__(self):
        return self.len

def bin(y, num_bins, xlog = False):
    if not xlog:
        pts_per_bin = y.size // num_bins
        y = y[:num_bins * pts_per_bin]
    
        return np.arange(0, len(y), pts_per_bin), y.reshape((num_bins, pts_per_bin))
    else:
        # log bins
        # kind of a quick hack
        # this turns len=16,bins=3 into 1,2.57,6.66 instead of 1,3,7, and rounds
        binstarts = np.round(np.geomspace(1, len(y) + 1, num_bins + 1) - 1).astype(np.int)
        return binstarts[:-1], [
                y[idx1:idx2]
                for idx1, idx2
                in zip(binstarts[:-1],binstarts[1:])
        ]
    

def bin_stats(y, num_bins, **kwparams):

    x, y = bin(y, num_bins, **kwparams)

    if isinstance(y, np.ndarray):
        return np.array([
            y.mean(axis=1),
            y.std(axis=1),
            y.min(axis=1),
            y.max(axis=1),
            x
        ])
    else:
        print(x,y)
        result = np.array([
            [
                bin.mean(axis=0),
                bin.std(axis=0),
                bin.min(axis=0),
                bin.max(axis=0),
                idx
            ]
            for idx, bin in zip(x, y)
        ])
        return result.transpose((1,0))

import matplotlib.pyplot as plt
raspifb_canvas = plt.figure(figsize=(8,4.8),dpi=100).canvas
def plot_bin_stats(y, num_bins, **kwparams):
    import matplotlib.pyplot as plt
    if len(y) >= num_bins * num_bins:
        means, stds, mins, maxs, xs = bin_stats(y, num_bins, **kwparams)
        plt.errorbar(xs, means, stds, fmt='ok', lw=3)
        plt.errorbar(xs, means, [means - mins, maxs - means], fmt='.k', ecolor='gray', lw=1)
    else:
        plt.plot(y)#, fmt='ok')

def raspifb_spectrum(data):
    img = np.array([[[0]*2]*800]*480, dtype=np.uint8)
    data = data#np.log(data)
    datamax = max(data)
    datamin = min(data)
    #print(datamin,datamax)
    data = ((data / datamax) * 479).astype(np.uint16)
    #print('then', min(data), max(data))
    for idx, item in enumerate(data):
        #print(item)
        #item = 240
        item = 479-item
        img[:item,idx] = [0x00,0x00]
        img[item:-1,idx] = [0xff,0xff]
    with open('/dev/fb0', 'wb') as fb:
        fb.write(img)

def raspifb_plt_show():
    raspifb_canvas.draw()
    img = raspifb_canvas.buffer_rgba()

    # conversion to 16 bits.  most fbs are 32 and don't need this.
    img = np.array(img)
    img[:,:,1] = (img[:,:,1] << 4) | img[:,:,2]
    img = img[:,:,:2].tobytes()

    with open('/dev/fb0', 'wb') as fb:
        fb.write(img)

class simulation:
    # untested
    def __init__(self, name):
        fn = name + '.sim'
        try:
            self.file = open(fn, 'rb')
        except FileNotFoundError:
            with open(fn, 'wb'):
                pass
            self.file = open(fn, 'rb')
    def read(self):
        return np.load(self.file)
    def write(self, data):
        np.save(self.file, data)
