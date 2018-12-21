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

plt.ylabel("P(sat)")
plt.xlabel("alpha")

def folder_from_N(N):
  if N != 3000:
    return 'data/survey_n100_N{}_M{}_SS1000_SW{}_p30'.format(N, N*5, N*100)
  else:
    return 'data/survey_n50_N{}_M{}_SS1000_SW{}_p30'.format(N, N*5, N*100)

N_values = [100, 300, 1000, 3000]
colors = [(0., 0.25, 0.5), (0., 0.5, 0.25), (0.5, 0.25, 0.), (0.5, 0., 0.25)]

# read input from each file, at each N
for i in range(len(N_values)):
  N = N_values[i]
  cur_color_unk = colors[i]
  cur_color_sat = (cur_color_unk[0] * 2, cur_color_unk[1] * 2, cur_color_unk[2] * 2)

  d = folder_from_N(N)

  # TODO sorting
  out_a = []    # alpha (x axis)
  out_sat = []  # proportion successfully sat
  out_unc = []  # proportion exceeding 1000 iteration limit
  out_unk = []  # proportion that failed after walksat
  with open(d + '/outcomes.txt') as out_f:
    out_reader = csv.DictReader(out_f, delimiter=' ')
    for row in out_reader:
      out_a += [float(row['M']) / N]
      out_sat += [float(row['sat']) / 100]
      out_unk += [float(row['unknown']) / 100]

  out_a = np.array(out_a)
  out_sat = np.array(out_sat)
  out_unc = np.array(out_sat)
  out_unk = np.array(out_sat)

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  plt.plot(out_a, out_sat, color=cur_color_sat, linewidth=1)
  plt.plot(out_a, out_unk, color=cur_color_unk, linewidth=1)

plt.savefig('plot/survey_sat_unk.png')
plt.show()
