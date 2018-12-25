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

N_values = [100, 300, 1000, 3000]
colors = [(0., 0.25, 0.5), (0., 0.5, 0.25), (0.5, 0.25, 0.), (0.5, 0., 0.25)]

# read input from each file, at each N
N = 3000
cur_color_sat = (0., 0., 0.)

d = 'data/survey_phase_n50_N3000_SS1000_SW300000_p30'

# TODO sorting
out_a = []    # alpha (x axis)
out_sat = []  # proportion successfully sat
out_unc = []  # proportion exceeding 1000 iteration limit
out_unk = []  # proportion that failed after walksat
with open(d + '/outcomes.txt') as out_f:
  out_reader = csv.DictReader(out_f, delimiter=' ')
  for row in out_reader:
    out_a += [float(row['M']) / N]
    out_sat += [float(row['sat']) / 50]
    out_unk += [float(row['unknown']) / 50]

out_a = np.array(out_a)
out_sat = np.array(out_sat)
out_unc = np.array(out_sat)
out_unk = np.array(out_sat)

# plot with some color
# plt.plot arguments: x-values, y-values, color, width
plt.plot(out_a, out_sat, color=cur_color_sat, linewidth=1)
#plt.plot(out_a, out_unk, color=cur_color_unk, linewidth=1)

plt.savefig('plot/phase_sat.png')
plt.show()
