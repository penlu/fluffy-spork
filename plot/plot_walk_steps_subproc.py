import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import subprocess
import sys

# plotting cosmetics
settings = {
  'font.size' : 12,
  'xtick.labelsize': 8,
  'ytick.labelsize': 8
}
plt.rcParams.update(settings)

plt.figure(figsize=(8,4), dpi=240)

plt.ylabel("steps")
plt.xlabel("alpha")

def folder_from_N(N):
  return 'data/walksat_n100_N{}_M{}_S{}_p30'.format(N, N*5, N*100)

N_values = [100, 300, 1000, 3000]
colors = [(0., 0.25, 0.5), (0., 0.5, 0.25), (0.5, 0.25, 0.), (0.5, 0., 0.25)]

# step-finder regex
p_st = re.compile('walk: ([0-9]*) steps')
p_sat = re.compile('walk: sat')
devnull = open('/dev/null', 'w')

# go over directories
for i in range(len(N_values)):
  N = N_values[i]
  cur_color = colors[i]

  d = folder_from_N(N)

  # TODO sorting
  out_a = []  # alpha (x axis)
  out_s = []  # number of steps taken: collection for each alpha
  for M in range(0, N*5 + N/10, N/10):
    print N, M
    ss = []

    for run in range(1, 101):
      fn = d + '/out_' + str(M) + '/' + str(run) + '.txt'

      try:
        grep_out = subprocess.check_output(['grep', 'walk: [0-9]* steps', fn])
      except subprocess.CalledProcessError, e:
        grep_out = e.output

      if subprocess.call(['grep', 'walk: sat', fn], stdout=devnull, stderr=devnull) == 0:
        m = p_st.findall(str(grep_out))
        if m:
          steps = int(m[0])

          ss += [float(steps) / (N * 100)]
        else:
          print "??", run, N, M
          print str(grep_out)

    if ss:
      out_a += [float(M) / N]
      out_s += [ss]

  out_avg = list(map(lambda ss: sum(ss)/len(ss), out_s))

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  plt.plot(out_a, out_avg, color=cur_color, linewidth=1)

plt.hlines(1, 0, 5)
plt.savefig('plot/walk_steps.png')
plt.show()
