from math import exp
import cashMoneyFunctions as CMF
import numpy as np
import scipy.optimize as root
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import MonteCarlo as MC
import time
sns.set()

# Part 1 - Cash money
print("Question 1")

# Parameters from question text
r0 = -0.005 
Ts = [1, 5, 10, 20, 30, 50] 
# Ts = [1, 10, 50]

# Parameters from Andreasen & Brøgger p. 9-10
theta = 0.0291
kappa = 0.0325
sigma = 0.0069
b_star = 0.0259

# Output related to question 1
for T in Ts:
    print("P(0, %s) = %s" %(T, CMF.P(r0, theta, sigma, kappa, T)))


print("Question 2")
paths = 100000

# The distribution of the integral is given by 
def simulateZCB(r0, theta, kappa, sigma, mat, paths):
    ekt = np.exp(-kappa*mat)
    W = np.random.normal(0, 1, paths)

    mean = mat*theta+(r0-theta)/kappa*(1-ekt)
    var = sigma*sigma/(kappa*kappa*kappa)*(kappa*mat+2*ekt-0.5*np.exp(-2*kappa*mat)-1.5)

    rates = np.exp(- (mean+np.sqrt(var)*W))
    return (np.mean(rates), np.std(rates), np.mean(rates)+1.96*np.std(rates)/np.sqrt(paths), np.mean(rates)-1.96*np.std(rates)/np.sqrt(paths))


maturities = range(1, 51)
t0 = time.time() 
data = [simulateZCB(r0, theta, kappa, sigma, x, paths) for x in maturities]
t1 = time.time()

total1 = t1-t0

for T in Ts:
    print("P(0, %s): mean: %s, 95%%-CI: [%s, %s] " %(T, data[T-1][0], data[T-1][3], data[T-1][2]))


# A lot of plotting
plotdata = pd.DataFrame([{'price': p, 'std': std,'+confidence': c1, '-confidence': c2} for p, std, c1, c2 in data])
plotdata['Maturities'] = maturities
plotdata['Analytical'] = [CMF.P(r0, theta, sigma, kappa, x) for x in maturities]

fig, ax = plt.subplots()
ax.plot(plotdata['Maturities'], plotdata['price'], label = "Simulated")
ax.plot(plotdata['Maturities'], plotdata['Analytical'], label = "Analytical")
ax.fill_between(plotdata['Maturities'], plotdata['-confidence'], plotdata['+confidence'], alpha=0.8)
plt.xlabel('Maturities')
plt.ylabel('Price of Zero Coupon Bond')
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(plotdata['Maturities'], plotdata['price']-plotdata['Analytical'])
ax.fill_between(plotdata['Maturities'], plotdata['-confidence']-plotdata['Analytical'], plotdata['+confidence']-plotdata['Analytical'], alpha=0.2)
plt.xlabel('Maturities')
plt.ylabel('Deviation from analytical solution')
plt.show()

fig, ax = plt.subplots()
ax.plot(plotdata['Maturities'], 100*(plotdata['Analytical']-plotdata['price'])/plotdata['Analytical'])
plt.xlabel('Maturities')
plt.ylabel('Relative Deviation(%)')
plt.show()


# The same thing using antithetic sampling
def simulateZCB_AT(r0, theta, kappa, sigma, mat, paths):
    ekt = np.exp(-kappa*mat)
    W = np.random.normal(0, 1, paths)

    mean = mat*theta+(r0-theta)/kappa*(1-ekt)
    var = sigma*sigma/(kappa*kappa*kappa)*(kappa*mat+2*ekt-0.5*np.exp(-2*kappa*mat)-1.5)
    r1 = np.exp(-(mean+np.sqrt(var)*W))
    r2 = np.exp(-(mean-np.sqrt(var)*W))
    rates = (r1 + r2) / 2
    return (np.mean(rates), np.std(rates), np.mean(rates)+1.96*np.std(rates)/np.sqrt(paths), np.mean(rates)-1.96*np.std(rates)/np.sqrt(paths))

maturities = range(1, 51)
t0 = time.time()
data2 = [simulateZCB_AT(r0, theta, kappa, sigma, x, paths) for x in maturities]
t1 = time.time()
total2 = t1-t0

print("Standard sampling time: %s, antithetic sampling time: %s, ratio: %s" %(total1, total2, total2 / total1))

for T in Ts:
    print("Standard sample variance: %s, antithetic sample variance: %s, efficiency ratio: %s" %(round(data[T-1][1]**2, 7), data2[T-1][1]**2, round(data[T-1][1]**2 / data2[T-1][1]**2, 5)))


# print("Question 3")

# plotdata1 = [100 * (CMF.CinfIntegral(r, theta, sigma, kappa, b_star, np.inf) - 1) for r in np.arange(-0.05, 0.05, 0.002)]
# sigma2 = 0.0075
# # b_hat = root.bisect(CMF.findB, -0.02, 0.05, args=(theta, sigma, kappa)) / Finds b_hat to be equal to value below
# b_hat = 0.03767795382955229
# plotdata2 = [100 * (CMF.CinfIntegral(r, theta, sigma2, kappa, b_hat, np.inf) - 1) for r in np.arange(-0.05, 0.05, 0.002)]

# # Create the figure
# fig, ax = plt.subplots()
# ax.plot(np.arange(-5, 5, 0.2), plotdata1, color='blue', label = "Cash premium, sigma = {}".format(sigma))
# ax.plot(np.arange(-5, 5, 0.2), plotdata2, color='orange', label = "Cash premium, sigma = {}".format(sigma2))
# plt.axvline(x = b_star * 100, linestyle='--', color='blue', label = "b* = {}, sigma = {}".format(b_star, sigma))
# plt.axvline(x = b_hat * 100, linestyle='--', color='orange', label = "b* = {}, sigma = {}".format(round(b_hat, 4), sigma2))
# plt.xlabel('Nominal short rate (%)')
# plt.ylabel('Premium (%)')
# plt.legend()
# plt.show()


# print("Question 4")

# maxSigma = 0.00779
# minSigma = 0.0001
# numOfPoints = 50
# sigmaRange = np.arange(minSigma, maxSigma, (maxSigma - minSigma)/numOfPoints)

# plotdata = [100 * (CMF.CinfIntegral(0, theta, sigma, kappa, root.bisect(CMF.findB, 0, 0.1, args=(theta, sigma, kappa)), np.inf) - 1) for sigma in sigmaRange]

# # Create the figure
# fig, ax = plt.subplots()
# ax.plot(sigmaRange*100, plotdata)
# plt.xlabel('Short rate volatility (%)')
# plt.ylabel('Premium (%)')
# plt.show()


# print("Question 5")

# rs = np.arange(-0.05, 0.03, 0.01)
# dt = 1/20
# paths = 100000 #should probably be lowered before running

# def MC_r_to_b(r):
#     discount, ts = CMF.sim_r_to_b(r, theta, kappa, sigma, b_star, paths, dt) # dt should be monthly - since parameters are yearly, dt should be 1/12
#     return (np.mean(discount), np.std(discount), np.mean(discount) + 1.96 * np.std(discount)/np.sqrt(paths), np.mean(discount) - 1.96 * np.std(discount)/np.sqrt(paths), ts)

# data = [MC_r_to_b(r) for r in rs]

# # Data processing
# plotdata = pd.DataFrame([{'price': p, 'std': std,'+confidence': c1, '-confidence': c2, 'ts': ts} for p, std, c1, c2, ts in data])
# plotdata['analytical'] = [CMF.CinfIntegral(r, theta, sigma, kappa, b_star, np.inf) for r in rs]
# fig, ax = plt.subplots()
# ax.plot(rs * 100, plotdata['price'])
# ax.plot(rs * 100, plotdata['analytical'])
# ax.fill_between(rs * 100, plotdata['-confidence'], plotdata['+confidence'], alpha=0.3)
# plt.xlabel('Nominal short rate (%)')
# plt.ylabel('Premium (%)')
# plt.show()

# i = 0
# for r in rs:
#     print("r(0) = %s, mean price: %s, 95%%-CI: [%s, %s], analytical: %s" \
#     %(round(r, 2), round(plotdata['price'][i], 5), round(plotdata['-confidence'][i], 5), round(plotdata['+confidence'][i], 5), round(plotdata['analytical'][i], 5)))    
#     i += 1


# print("Question 6")

# # Plots the individual distributions
# for ts in plotdata['ts']:
#     mean = np.mean(ts) * dt
#     std  = np.std(ts) * dt / np.sqrt(paths)
#     print("mean: %s, 95%%-CI: [%s, %s]" %(mean, mean - 1.96*std, mean + 1.96*std))
#     sns.distplot(np.array(ts) * dt)
#     plt.axvline(x = mean, linestyle='--', color='blue', label = "tau* = {}".format(mean))
#     plt.xlabel('Years')
#     plt.title('Distribution of when barrier is hit')
#     plt.legend()
#     plt.show()

# # Plots every other distribution in the same figure
# i = 0
# for ts in plotdata['ts']:
#     if (i % 2 == 0):
#         sns.distplot(np.array(ts) * dt)
#     i+=1

# plt.xlabel('Years')
# plt.title('Distribution of when barrier is hit')
# plt.show()


# print("Question 7")
# def MC_r_to_b_cv(r):
#     discount, control, ts = CMF.sim_r_to_b_cv(r, theta, kappa, sigma, b_star, paths, dt) # dt should be monthly - since parameters are yearly, dt should be 1/12
#     return (discount, control)

# data = [MC_r_to_b_cv(r) for r in rs]

# r_index = 0
# for discount, control in data:
#     mean1 = np.mean(discount)
#     mean2 = np.mean(control)

#     numerator = 0
#     denominator = 0
#     for i in range(len(discount)):
#         numerator += (discount[i] - mean1)*(control[i] - mean2)
#         denominator += (control[i] - mean2)**2

#     betaHat = numerator / denominator
#     analytical = CMF.P(rs[r_index], theta, sigma, kappa, 200)
#     pcv = list()
#     for i in range(len(discount)):
#         pcv.append(discount[i] + betaHat * (analytical - control[i]))

#     print("r(0) = %s, price = %s, 95%%-CI: [%s, %s]" %(round(rs[r_index], 2), np.mean(pcv), np.mean(pcv) - 1.96 * np.var(pcv)/paths, np.mean(pcv) + 1.96 * np.var(pcv)/paths))
#     print("Without control variable: %s, with control variable: %s, efficiency ratio: %s" %(np.var(discount), np.var(pcv), np.var(discount) / np.var(pcv)))
#     r_index += 1


# # # Part 2 - Finite difference goes American
# from CN import CrankNicolson as CN # Holds the Crank_Nicolson solver class
# # Constants given in question text
# r = 0.04  # risk free
# q = 0.02  # dividend yield
# S = 100  # underlying spot
# sigma = 2  # volatility
# T = 1.0  # maturity

# # print("Question 1")
# K = 100  # strike
# solver = CN(r, q, sigma, K, T) #Instance of the solver class
# solver.max_dt = 0.01 
# solver.N = 800 # Number of points
# solver.american = False # Solves the European


# price = solver.solve(S)
# print(price)

# # Plot the smile
# x = solver.S.flatten()
# y = solver.X.flatten()
# plt.xlabel('Strike')
# plt.ylabel('Price')
# plt.plot(x,y)
# plt.show()

# print("Question 2")
# solver.american = True # Solves the American
# price = solver.solve(S)
# print(price)

# # Plot the smile
# x = solver.S.flatten()
# y = solver.X.flatten()
# plt.xlabel('Strike')
# plt.ylabel('Price')
# plt.plot(x,y)
# plt.show()


# print("Question 3")
# S = 100
# K = 100  # strike
# prices = [45]
# def find_root(sigma, solver, S, price):
#     solver.sigma = sigma
#     return price - solver.solve(S)

# for price in prices:
#     impVol = root.newton(find_root, 1, args=(solver, S, price))
#     print("price = %s, impVol = %s" %(price, impVol))

