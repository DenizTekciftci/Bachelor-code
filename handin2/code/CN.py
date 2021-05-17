import numpy as np
import matplotlib.pyplot as plt


class CrankNicolson:
    def __init__(self, r, q, sigma, K, T):
        self.r = r
        self.q = q
        self.sigma = sigma
        self.K = K
        self.T = T
        self.S_upper = 3 * K
        self.S = []
        self.X = []
        self.A = []
        self.b = []
        self.N = 200
        self.max_dt = 1/12.0
        self.american = False
        self.eps = 1e-5
        self.max_iter = 200
        self.omega = 1.2
        self.err = 0
        self.iter = 0

    def solve(self, S0):
        self.setInitialCondition()

        self.solvePDE()

        x = self.S.flatten()
        y = self.X.flatten()
        return np.interp(S0, x, y)

    def solvePDE(self):
        t = self.T
        while t > 0:
            dt = min(t, self.max_dt)
            self.setCoeff(dt)
            if self.american:
                self.solvePSOR()
            else:
                self.solveLinearSystem()
            t -= dt

    def setInitialCondition(self):
        self.S = np.linspace(0, self.S_upper, self.N)
        self.A = np.zeros((self.N, self.N))
        self.b = np.zeros((self.N, 1))
        self.X = np.maximum(self.K - self.S, 0)


    def setCoeff(self, dt):
        N = self.N
        r = self.r
        q = self.q
        S = self.S
        X = self.X
        sigma = self.sigma
        dS = S[1] - S[0]
        for i in range(0, N-1):
            alpha = 0.25 * dt * (np.square(sigma*S[i]/dS) - (r - q) * S[i]/dS)
            beta = 0.5 * dt * (r + np.square(sigma * S[i]/dS))
            gamma = 0.25 * dt * (np.square(sigma*S[i]/dS) + (r - q) * S[i]/dS)
            if i == 0:
                self.b[i] = X[i] * (1 - beta)
                self.A[i][i] = 1 + beta
            else:
                self.b[i] = alpha * X[i-1] + (1 - beta) * X[i] + gamma * X[i+1]
                self.A[i][i-1] = -alpha
                self.A[i][i] = 1 + beta
                self.A[i][i+1] = -gamma
        self.A[-1][N-4] = -1
        self.A[-1][N-3] = 4
        self.A[-1][N-2] = -5
        self.A[-1][N-1] = 2
        self.b[-1] = 0

    def solveLinearSystem(self):
        self.X = np.linalg.solve(self.A, self.b)

    def solvePSOR(self):
        N = self.N
        iter = 0
        omega = self.omega
        self.err = 1000
        while self.err > self.eps and iter < self.max_iter:
            iter += 1
            x_old = self.X.copy()
            for i in range(N-1):
                self.X[i] = (1 - omega) * self.X[i] + omega * self.b[i] / self.A[i][i]
                self.X[i] -= self.A[i][i+1] * self.X[i+1] * omega / self.A[i][i]
                self.X[i] -= self.A[i][i-1] * self.X[i-1] * omega / self.A[i][i]

            #for last row, use boundary condition
            self.X[N-1] = (1 - omega) * self.X[i] + omega * self.b[i] / self.A[i][i]
            for j in range(N-4, N):
                self.X[N-1] -= self.A[N-1][j] * self.X[j] * omega / self.A[N-1][N-1]

            self.applyConstraint()
            self.err = np.linalg.norm(x_old - self.X)
            self.iter = iter

    def applyConstraint(self):
        self.X = np.maximum(self.X, self.K - self.S)
