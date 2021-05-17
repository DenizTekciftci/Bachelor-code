import math
import numpy as np
from scipy.stats import norm  
from scipy.integrate import quad

def A(theta, sigma, kappa, T):
    temp_B = B(kappa, T)
    result = (theta - sigma * sigma / (2 * kappa * kappa)) * (temp_B - T) - sigma * sigma / (4 * kappa) * temp_B * temp_B
    return result

def B(kappa, T):
    result = (1 - math.exp(-kappa * T))/kappa
    return result

def P(r, theta, sigma, kappa, T):
    return math.exp(A(theta, sigma, kappa, T) - B(kappa, T) * r)

def v(u, sigma, kappa):
    return math.sqrt(sigma**2 / (2*kappa) * (1 - np.exp(-2 * kappa * u)))

def f(r, theta, sigma, kappa, u):
    Aprime = -(np.exp(-u * kappa) * (1 - np.exp(-u * kappa)) * sigma**2)/(2 * kappa**2) + (-1 + np.exp(-u * kappa)) * (theta - sigma**2/(2 * kappa**2))
    Bprime = np.exp(-u * kappa)
    return  Bprime * r - Aprime

def CinfIntegrand(u, r, theta, sigma, kappa, b):
    temp_v = v(u, sigma, kappa)
    temp_f = f(r, theta, sigma, kappa, u)
    quantile = (temp_f - b) / temp_v
    return P(r, theta, sigma, kappa, u) * (temp_f * norm.cdf(quantile) + temp_v * norm.pdf(quantile))

def CinfIntegral(r, theta, sigma, kappa, b, tau):
    return quad(CinfIntegrand, 0, tau, args=(r, theta, sigma, kappa, b))[0]

def findB(b, theta, sigma, kappa):
    return quad(CinfIntegrand, 0, np.inf, args=(b, theta, sigma, kappa, b))[0] - 1

def sim_r_to_b(r0, theta, kappa, sigma, b, paths, dt):
    rs = list()
    ts = list()
    for _ in range(paths):
        r = r0
        t = 0
        integral = 0
        while(True):
            dr = kappa * (theta - r) * dt + sigma * np.random.normal(0, np.sqrt(dt))
            r += dr
            t += dt
            integral += r * dt
            if (r >= b):
                break
        ts.append(t)
        rs.append(np.exp(-integral))
    return rs, ts


def sim_r_to_b_cv(r0, theta, kappa, sigma, b, paths, dt):
    rs = list()
    rcvs = list()
    ts = list()
    for i in range(paths):
        r = r0
        t = 0
        integral = 0
        barrier = False
        time = False
        while(not (barrier and time)):
            dr = kappa * (theta - r) * dt + sigma * np.random.normal(0, np.sqrt(dt))
            r += dr
            t += dt
            integral += r * dt
            if (r >= b and not barrier):
                barrier = True
                ts.append(t)
                rs.append(np.exp(-integral))
            if (t >= 200 and not time):
                time = True
                rcvs.append(np.exp(-integral))
    return rs, rcvs, ts