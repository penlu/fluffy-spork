import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import subprocess
import sys

# plotting cosmetics
def set_settings():
  settings = {
    'font.size' : 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
  }
  plt.rcParams.update(settings)

  plt.figure(figsize=(8,4), dpi=240)

  plt.ylabel("steps")
  plt.xlabel("alpha")

N_values = [100, 300, 1000, 3000]
colors = [(0., 0.25, 0.5), (0., 0.5, 0.25), (0.5, 0.5, 0.), (0.5, 0., 0.25)]

# cloud and line per N, separately
for i in range(len(N_values)):
  cur_color = colors[i]

  out_scatter_a = np.load('plot/walk_steps_{}_scatter_a.npy'.format(i))
  out_scatter_s = np.load('plot/walk_steps_{}_scatter_s.npy'.format(i))
  out_a = np.load('plot/walk_steps_{}_avg_a.npy'.format(i))
  out_avg = np.load('plot/walk_steps_{}_avg.npy'.format(i))

  set_settings()
  plt.scatter(out_scatter_a, out_scatter_s, c=cur_color, marker='.', s=0.1)
  plt.plot(out_a, out_avg, color=cur_color, linewidth=1)
  plt.hlines(1, 0, 5)
  plt.savefig('plot/walk_steps_scatter_{}.png'.format(N_values[i]))
  plt.show()

# lines per N, all together
set_settings()
for i in range(len(N_values)):
  cur_color = colors[i]

  out_a = np.load('plot/walk_steps_{}_avg_a.npy'.format(i))
  out_avg = np.load('plot/walk_steps_{}_avg.npy'.format(i))

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  plt.plot(out_a, out_avg, color=cur_color, linewidth=1)

plt.hlines(1, 0, 5)
plt.savefig('plot/walk_steps_avgs.png')
plt.show()
