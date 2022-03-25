#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import numpy as np

file = open("../user/out.dat", "r")

x, y, u, id = [], [], [], []

for l in file:
    row = l.split()
    x.append(float(row[0]))
    y.append(float(row[1]))
    u.append(float(row[2]))
    id.append(int(row[3]))

x = np.array(x)
y = np.array(y)
u = np.array(u)
id = np.array(id)

f, ax = plt.subplots(1, 1, sharex=True, sharey=True)
tpf = ax.scatter(x, y, s=20, c=u, marker='o', cmap=cm.rainbow)

# texts = []
# for i, txt in enumerate(id):
#     texts.append(ax.text(x[i], y[i], txt, fontsize=5))
#
# adjust_text(texts)

# yoffset = 0.1
# switch = -1
# for i, txt in enumerate(id):
#     ax.annotate(txt, xy=(x[i], y[i]), xytext=(x[i], y[i]+switch*yoffset))
#     switch*=-1

# tpf = ax.tripcolor(x, y, u, cmap="RdBu_r")
# ax[1].tricontourf(x, y, u, 10, cmap='rainbow')  # choose 20 contour levels, just to show how good its interpolation is
# ax[1].plot(x, y, 'ko ')
# ax[0].plot(x, y, 'ko ')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_trisurf(x, y, z, cmap='rainbow', linewidth=0)
# fig.colorbar(surf)

plt.colorbar(tpf)
plt.show()

file.close()
