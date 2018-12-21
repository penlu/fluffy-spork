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

N = N_values[3]

M_min = 0
M_max = N*5

series_count = 0
for M in range(M_min, M_max + N/10, N/10):
  try:
    avg_sigma = np.load('plot/survey_pre/survey_{}_avg_sigma.npy'.format(M))
    series_count += 1
  except Exception, e:
    continue

series_count -= 1

# lines per M, all together
set_settings()
color_counter = 0
for M in range(M_min, M_max + N/10, N/10):
  if series_count != 0:
    color_ratio = float(color_counter) / series_count
  else:
    color_ratio = 0
  cur_color_avg = (0., math.cos(color_ratio * math.pi/2)/3, math.sin(color_ratio * math.pi/2)/3)
  cur_color_sample = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))
  #cur_color_avg = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

  try:
    step_sigma = np.load('plot/survey_pre/survey_{}_step_sigma.npy'.format(M))
    avg_sigma = np.load('plot/survey_pre/survey_{}_avg_sigma.npy'.format(M))
    color_counter += 1
  except Exception, e:
    continue

  points_x = []
  points_y = []
  # destack using numpy?
  for i in range(len(step_sigma)):
    for j in range(len(step_sigma[i])):
      points_x += [float(i) / N]
      points_y += [step_sigma[i][j]]

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  size = len(avg_sigma)
  plt.scatter(np.linspace(0, float(size - 1) / N, size), avg_sigma, c=cur_color_avg, marker='.', s=0.1)
  plt.scatter(np.array(points_x), np.array(points_y), c=cur_color_sample, marker='.', s=0.01)

plt.savefig('plot/survey_avg_sigma.png')
plt.show()

set_settings()
color_counter = 0
for M in range(M_min, M_max*5 + N/10, N/10):
  if series_count != 0:
    color_ratio = float(color_counter) / series_count
  else:
    color_ratio = 0
  cur_color_avg = (0., math.cos(color_ratio * math.pi/2)/3, math.sin(color_ratio * math.pi/2)/3)
  cur_color_sample = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))
  #cur_color_avg = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

  try:
    step_polar= np.load('plot/survey_pre/survey_{}_step_polar.npy'.format(M))
    avg_polar = np.load('plot/survey_pre/survey_{}_avg_polar.npy'.format(M))
    color_counter += 1
  except Exception, e:
    continue

  points_x = []
  points_y = []
  # destack using numpy?
  for i in range(len(step_polar)):
    for j in range(len(step_polar[i])):
      points_x += [float(i) / N]
      points_y += [step_polar[i][j]]

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  size = len(avg_polar)
  plt.scatter(np.linspace(0, float(size - 1) / N, size), avg_polar, c=cur_color_avg, marker='.', s=0.1)
  plt.scatter(np.array(points_x), np.array(points_y), c=cur_color_sample, marker='.', s=0.01)

plt.savefig('plot/survey_avg_polar.png')
plt.show()

