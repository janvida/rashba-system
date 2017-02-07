import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm



def calculate_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    return dx / dy

# Read the array from disk. Note that this returns a 2D array.
data = np.loadtxt('400x10_sz_T0.05.txt')

data = data.reshape((10, 400, 81))

# Store the data into individual matrices
g_ul_p = data[0, :, :]
g_ul_m = data[1, :, :]
g_ur_p = data[2, :, :]
g_ur_m = data[3, :, :]
g_rl_p = data[4, :, :]
g_rl_m = data[5, :, :]
v_p = data[6, :, :]
v_m = data[7, :, :]
dv = data[8, :, :]
dos = data[9, :, :]

print type(dv)
print dv.transpose().shape
extent = (0, 400, 0, 81)
shape = dv.transpose().shape

plt.matshow(dv.transpose(),
           interpolation='nearest',extent=extent, aspect=calculate_aspect(shape, extent), origin = 'lower', cmap=cm.gist_rainbow)
plt.colorbar()
plt.show()
