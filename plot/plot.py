### TEMPLATE FOR REFERENCE

import numpy as np
import matplotlib.pyplot as plt
import sys

# read input
in_data = ""
for line in sys.stdin:
  in_data += line

# convert to list
alldata = eval(in_data)

# go which proportion of the way to white
def lerp(color, amt):
  return (color[0]*(1 - amt) + amt, color[1]*(1 - amt) + amt, color[2]*(1 - amt) + amt)
allcolors = [(0., 0., 0.5), (0.5, 0.25, 0.), (0., 0.5, 0.)]

# plot with some color
### ALLDATA IS AN ARRAY OF SERIES
### CUR_SERIES IS ONE LONG ARRAY OF Y-VALUES TO BE PLOTTED
for condition in range(len(alldata)):
  cur_color = allcolors[condition]
  cur_series = alldata[condition]
  for i in range(len(cur_series)):
    s = np.array(cur_series[i])
    # plt.plot arguments: x-values, y-values, color, width
    plt.plot(np.linspace(1, s.shape[0], s.shape[0]), np.sqrt(s), color=lerp(cur_color, 0), linewidth=3)
plt.ylabel("eigenvalues")
plt.xlabel("h")
plt.show()
