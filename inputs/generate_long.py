import numpy as np

# Generates a grid that will result in a single contour, if interval is set to 0.6.
n = 100
m = 100
z = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        if i == 0:
            z[i][j] = 0.1
        elif i == n - 1:
            z[i][j] = 0.1
        elif j == 0:
            z[i][j] = 0.1
        elif j == 1:
            z[i][j] = 1
        elif j == m - 1:
            z[i][j] = 0.1
        elif i % 2 == 0:
            z[i][j] = 0.1
        else:
            z[i][j] = 1

header = "ncols {}\n".format(z.shape[1])
header += "nrows {}\n".format(z.shape[0])
header += "xllcorner {}\n".format(0)
header += "yllcorner {}\n".format(0)
header += "cellsize {}".format(0.0001)
np.savetxt("long.txt", z, header=header, comments='')
