import numpy as np
import math 



def simulateR2(r0: float, theta: float, kappa: float, sigma: float, T: float, paths: float, steps: float):
    dt = T / steps
    normals = np.random.normal(0, dt, paths*steps).reshape(paths, steps)
    rs = list()
    for i in range(paths - 1):
        integral = 0
        for j in range(steps - 1):
            integral += math.exp(kappa * (j * dt)) * normals[i][j]
        rs.append(integral)
    return rs

def simulateIntegralR(r0: float, theta: float, kappa: float, sigma: float, T: float, paths: float, steps: float):
    dt = T / steps
    w = np.random.normal(0, np.sqrt(dt), paths*steps).reshape(paths, steps)
    irs = list()
    avgr = 0
    for i in range(paths - 1):
        r = r0
        integral_r = 0
        for j in range(steps - 1):
            r_1 = r
            r = r + kappa * (theta - r) * dt + sigma * w[i][j]
            integral_r += (r + r_1)/2 * dt
        avgr += r    
        irs.append(math.exp(-integral_r))
    # print("r = %s" %(avgr / paths))
    return irs

