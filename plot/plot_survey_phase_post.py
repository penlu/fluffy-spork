import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import subprocess
import sys

# plotting cosmetics
def set_settings(ylabel, xlabel):
  settings = {
    'font.size' : 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8
  }
  plt.rcParams.update(settings)

  plt.figure(figsize=(8,4), dpi=240)

  plt.ylabel(ylabel)
  plt.xlabel(xlabel)

N = 3000

series_count = 0
for M in range(12240, 12900, 60):
  try:
    avg_sigma = np.load('plot/phase_pre/phase_{}_avg_sigma.npy'.format(M))
    series_count += 1
  except Exception, e:
    continue

series_count -= 1

RAINBOW=False

# lines per M, all together
set_settings('sigma (nats)', 'fixed variables')
if RAINBOW:
  color_counter = 0
  for M in range(12240, 12900, 60):
    if series_count != 0:
      color_ratio = float(color_counter) / series_count
    else:
      color_ratio = 0
    cur_color_sample = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

    try:
      step_sigma = np.load('plot/phase_pre/phase_{}_step_sigma.npy'.format(M))
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
    plt.scatter(np.array(points_x), np.array(points_y), c=cur_color_sample, marker='.', s=0.01)

color_counter = 0
for M in range(12240, 12900, 60):
  if series_count != 0:
    color_ratio = float(color_counter) / series_count
  else:
    color_ratio = 0

  if RAINBOW:
    cur_color_avg = (0., math.cos(color_ratio * math.pi/2)/3, math.sin(color_ratio * math.pi/2)/3)
  else:
    cur_color_avg = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

  try:
    avg_sigma = np.load('plot/phase_pre/phase_{}_avg_sigma.npy'.format(M))
    color_counter += 1
  except Exception, e:
    continue

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  size = len(avg_sigma)
  plt.scatter(np.linspace(0, float(size - 1) / N, size), avg_sigma, c=cur_color_avg, marker='.', s=0.1)

if RAINBOW:
  plt.savefig('plot/phase_avg_sigma_rainbow.png')
else:
  plt.savefig('plot/phase_avg_sigma.png')
plt.show()

set_settings('polarization', 'fixed variables')
plt.ylim(bottom=0.75)
if RAINBOW:
  color_counter = 0
  for M in range(12240, 12900, 60):
    if series_count != 0:
      color_ratio = float(color_counter) / series_count
    else:
      color_ratio = 0
    cur_color_sample = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

    try:
      step_polar= np.load('plot/phase_pre/phase_{}_step_polar.npy'.format(M))
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
    plt.scatter(np.array(points_x), np.array(points_y), c=cur_color_sample, marker='.', s=0.01)

color_counter = 0
for M in range(12240, 12900, 60):
  if series_count != 0:
    color_ratio = float(color_counter) / series_count
  else:
    color_ratio = 0

  if RAINBOW:
    cur_color_avg = (0., math.cos(color_ratio * math.pi/2)/3, math.sin(color_ratio * math.pi/2)/3)
  else:
    cur_color_avg = (0., math.cos(color_ratio * math.pi/2), math.sin(color_ratio * math.pi/2))

  try:
    avg_polar = np.load('plot/phase_pre/phase_{}_avg_polar.npy'.format(M))
    color_counter += 1
  except Exception, e:
    continue

  # plot with some color
  # plt.plot arguments: x-values, y-values, color, width
  size = len(avg_polar)
  plt.scatter(np.linspace(0, float(size - 1) / N, size), avg_polar, c=cur_color_avg, marker='.', s=0.1)

if RAINBOW:
  plt.savefig('plot/phase_avg_polar_rainbow.png')
else:
  plt.savefig('plot/phase_avg_polar.png')
plt.show()

