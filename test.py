import numpy as np
import matplotlib.pyplot as plt

# parameters
r = np.array([5, 1, 1, 1, 5], dtype=np.float64)

# initial condition
x = np.array([10, 10, 10, 10, 10], dtype=np.float64)

# function to compute residual
def Compute_R(x, r, n):
    xx = np.power(np.abs(x), n) * np.sign(x)
    R = np.array([
        r[1]*xx[1] - r[0]*xx[0] - r[2]*xx[2],
        r[4]*xx[4] - r[2]*xx[2] - r[3]*xx[3],
        -x[0] - x[1] + 10,
        x[0] - x[2] - x[3],
        x[1] + x[2] - x[4]
    ], dtype=np.float64)
    return R

# function to compute jacobian of residual
def Compute_J(x, r, n):
    xx = np.power(np.abs(x), n-1) * np.sign(x)
    J = np.array([
        [-n*r[0]*xx[0], n*r[1]*xx[1], -n*r[2]*xx[2], 0, 0],
        [0, 0, -n*r[2]*xx[2], -n*r[3]*xx[3], n*r[4]*xx[4]],
        [-1, -1, 0, 0, 0],
        [1, 0, -1, -1, 0],
        [0, 1, 1, 0, -1]
    ], dtype=np.float64)
    return J

# newton's method
N = 10000
R_r = np.empty((N, x.shape[0]), dtype=np.float64)
for i in range(N):
    R = Compute_R(x, r, 2)
    R_r[i] = R
    J = Compute_J(x, r, 2)
    dx = np.matmul(np.linalg.inv(J), R) # matmul will seem R as col vec
    x -= dx

print(x)
print(R_r[-1])

