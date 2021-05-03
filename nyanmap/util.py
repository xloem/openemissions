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

def bin(y, num_bins):
    pts_per_bin = y.size // num_bins
    y = y[:num_bins * pts_per_bin]

    return np.arange(0, len(y), pts_per_bin), y.reshape((num_bins, pts_per_bin))

def bin_stats(y, num_bins):

    x, y = bin(y, num_bins)

    return np.array([
        y.mean(axis=1),
        y.std(axis=1),
        y.min(axis=1),
        y.max(axis=1),
        x
    ])

import matplotlib.pyplot as plt
raspifb_canvas = plt.figure(figsize=(8,4.8),dpi=100).canvas
def plot_bin_stats(y, num_bins):
    import matplotlib.pyplot as plt
    if len(y) >= num_bins * 2:
        means, stds, mins, maxs, xs = bin_stats(y, num_bins)
        plt.errorbar(xs, means, stds, fmt='ok', lw=3)
        plt.errorbar(xs, means, [means - mins, maxs - means], fmt='.k', ecolor='gray', lw=1)
    else:
        plt.plot(y)#, fmt='ok')

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
