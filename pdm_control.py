import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# This is just automates all the commands needed to run
# the pdm program, just to make things a little cleaner.
# When compiled, ./pdm takes four arguments:
# data_file num_bins seg_width max_freq
# This is just an easy way to run it (albeit very naively).

data_file = 'OGLE-BLG-CEP-027.dat'
num_bins = 10
seg_width = 10
max_freq = 20

#  Runs PDM
os.system("make")
os.system("./pdm {0} {1} {2} {3}".format(data_file, 
                                         num_bins, 
                                         seg_width,
                                         max_freq))

# Loads files and plots the relevant plots
data = np.loadtxt(data_file)
theta = np.loadtxt('output_theta.dat')

# Computes best period
sns.set_style("whitegrid")
min_ind = np.argmin(theta[:,1])
period = 1/(theta[min_ind,0])
plt.tick_params(axis='both', labelsize=14)

# A couple of plots to make sure everything worked
plt.scatter(theta[:,0], theta[:,1])
plt.xlabel("Frequency", fontsize=16)
plt.ylabel("$\Theta$", fontsize=16)
plt.show()
plt.close()

plt.gca().invert_yaxis()
plt.scatter(data[:,0], data[:,1])
plt.xlabel("Time (days)", fontsize=16)
plt.ylabel("Apparent Magnitude", fontsize=16)
plt.show()
plt.close()

plt.gca().invert_yaxis()
plt.scatter(np.mod(data[:,0],period)/period, data[:,1])
plt.xlabel("Phase", fontsize=16)
plt.ylabel("Apparent Magnitude", fontsize=16)
plt.show()
plt.close()
