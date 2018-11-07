import numpy as np

dx, dy = 0.05, 0.05

y, x = np.mgrid[slice(1, 5 + dy, dy), slice(1, 5 + dx, dx)]

z = np.sin(x)**10 + np.cos(10 + y*x) * np.cos(x)

z = z[:-1,:-1]

np.savetxt("trig.txt", z, header="{} {}".format(z.shape[0], z.shape[1]), comments='')
