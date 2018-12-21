import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import subprocess
import sys

# step-finder regex
p_sigma = re.compile('survey: [0-9]* steps sigma/N ([+-.0-9e]*)')
p_polar = re.compile('survey: selected .*, pol ([+-.0-9e]*),.*')
devnull = open('/dev/null', 'w')

# only plot for 1k for now
N = 3000
d = 'data/survey_phase_n50_N3000_SS1000_SW300000_p30'

for M in range(12240, 12900, 60):
  # two lists of lists: first layer is by run, second is by step
  out_sigma = []  # instance complexity at each step
  out_polar = []  # polarization of last fixed var at each step
  out_lens = []  # lengths of these lists

  for run in range(1, 51):
    fn = d + '/out_' + str(M) + '/' + str(run) + '.txt'

    this_sigma = []
    this_polar = []
    if True or subprocess.call(['grep', 'survey: sat', fn], stdout=devnull, stderr=devnull) == 0:
      for l in open(fn):
        sig = p_sigma.findall(l)
        if sig:
          #print "found sigma {} {}".format(M, run)
          this_sigma += [float(sig[0])]

        pol = p_polar.findall(l)
        if pol:
          #print "found polar {} {}".format(M, run)
          this_polar += [float(pol[0])]

    #print "!!!", M, run
    assert(len(this_sigma) == len(this_polar) or len(this_sigma) == len(this_polar) + 1)

    # ending in trivial surveys, then no polarization prints
    # don't want to grep for more
    if len(this_sigma) == len(this_polar) + 1:
      this_sigma = this_sigma[:-1]

    if this_sigma:
      out_sigma += [this_sigma]
      out_polar += [this_polar]
      out_lens += [len(this_sigma)]

  print "progress:", M, len(out_lens)

  if len(out_lens) != 0:
    # transpose the lists
    step_sigma = [[] for i in range(max(out_lens))]
    step_polar = [[] for i in range(max(out_lens))]
    for step in range(max(out_lens)):
      for run in range(len(out_sigma)):
        if out_lens[run] > step:
          step_sigma[step] += [out_sigma[run][step]]
          step_polar[step] += [out_polar[run][step]]
    
    avg_sigma = list(map(lambda ss: sum(ss)/len(ss), step_sigma))
    avg_polar = list(map(lambda ss: sum(ss)/len(ss), step_polar))

    step_sigma = np.array(step_sigma)
    step_polar = np.array(step_polar)
    avg_sigma = np.array(avg_sigma)
    avg_polar = np.array(avg_polar)

    np.save('plot/phase_pre/phase_{}_step_sigma'.format(M), step_sigma)
    np.save('plot/phase_pre/phase_{}_step_polar'.format(M), step_polar)
    np.save('plot/phase_pre/phase_{}_avg_sigma'.format(M), avg_sigma)
    np.save('plot/phase_pre/phase_{}_avg_polar'.format(M), avg_polar)

