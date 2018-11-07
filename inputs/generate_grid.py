import numpy as np

dx, dy = 0.05, 0.05

y, x = np.mgrid[slice(1, 5 + dy, dy), slice(1, 5 + dx, dx)]

z = np.sin(x)**10 + np.cos(10 + y*x) * np.cos(x)

z = z[:-1,:-1]

header = "ncols {}\n".format(z.shape[1])
header += "nrows {}\n".format(z.shape[0])
header += "xllcorner {}\n".format(0)
header += "yllcorner {}\n".format(0)
header += "cellsize {}".format(0.0001)
np.savetxt("trig.txt", z, header=header, comments='')
