import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import subprocess
import sys

N_values = [100, 300, 1000, 3000]

def folder_from_N(N):
  return 'data/walksat_n100_N{}_M{}_S{}_p30'.format(N, N*5, N*100)

# step-finder regex
p_st = re.compile('walk: ([0-9]*) steps')
p_sat = re.compile('walk: sat')
devnull = open('/dev/null', 'w')

# go over directories
for i in range(len(N_values)):
  N = N_values[i]

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

  out_scatter_a = []
  out_scatter_s = []
  for ai in range(len(out_a)):
    for si in range(len(out_s[ai])):
      out_scatter_a += [out_a[ai]]
      out_scatter_s += [out_s[ai][si]]

  out_scatter_a = np.array(out_scatter_a)
  out_scatter_s = np.array(out_scatter_s)

  out_avg = list(map(lambda ss: sum(ss)/len(ss), out_s))

  out_a = np.array(out_a)
  out_avg = np.array(out_avg)

  np.save('plot/walk_pre/walk_steps_{}_scatter_a'.format(i), out_scatter_a)
  np.save('plot/walk_pre/walk_steps_{}_scatter_s'.format(i), out_scatter_s)
  np.save('plot/walk_pre/walk_steps_{}_avg_a'.format(i), out_a)
  np.save('plot/walk_pre/walk_steps_{}_avg'.format(i), out_avg)
