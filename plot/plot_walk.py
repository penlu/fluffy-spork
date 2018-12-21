import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import sys

# plotting cosmetics
settings = {
  'font.size' : 12,
  'xtick.labelsize': 8,
  'ytick.labelsize': 8
}
plt.rcParams.update(settings)

plt.figure(figsize=(8,4), dpi=240)

def folder_from_N(N):
  return 'data/walksat_n100_N{}_M{}_S{}_p30'.format(N, N*5, N*100)

N_values = [100, 300, 1000, 3000]
colors = [(0., 0.25, 0.5), (0., 0.5, 0.25), (0.5, 0.25, 0.), (0.5, 0., 0.25)]

# read input from each file, at each N
for i in range(len(N_values)):
  N = N_values[i]
  cur_color = colors[i]

  d = folder_from_N(N)

  # TODO sorting
  out_a = []    # alpha (x axis)
  out_sat = []  # proportion successfully sat
  with open(d + '/outcomes.txt') as out_f:
    out_reader = csv.DictReader(out_f, delimiter=' ')
    for row in out_reader:
      out_a += [float(row['M']) / N]
      out_sat += [float(row['sat']) / 100]

  out_a = np.array(out_a)
  out_sat = np.array(out_sat)

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  plt.plot(out_a, out_sat, color=cur_color, linewidth=1)

  plt.ylabel("P(sat)")
  plt.xlabel("alpha")

plt.savefig('plot/walk_sat.png')
plt.show()
