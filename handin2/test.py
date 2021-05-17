from math import exp
from numpy.core.fromnumeric import shape
import cashMoneyFunctions as CMF
import MonteCarlo as MC
import math
import random
import numpy as np
import scipy.optimize as root
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import poisson
sns.set()


# Parameters from question text
r0 = 0.025
Ts = [1, 5, 10, 20, 30, 50] 
# Ts = [1, 10, 50]

# Parameters from Andreasen & Brøgger p. 9-10
theta = 0.0291
kappa = 0.0325
sigma = 0.0069
b_star = 0.0259

maxSigma = 0.0077
minSigma = 0.0001
numOfPoints = 50
from out8 import bs
bs2 = np.array(bs)
ts = (bs2[:,0])
bs = bs2[:,1]

fig, ax = plt.subplots()
ax.plot(ts, bs)
plt.xlabel('t')
plt.ylabel('Interest rate')
plt.legend()
plt.show()