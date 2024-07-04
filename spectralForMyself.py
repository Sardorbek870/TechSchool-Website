import numpy as np
from scipy.linalg import solve


def chebyshev_points(N):
    return np.cos(np.pi * np.arange(N + 1) / N)


def chebyshev_differentiation_matrix(N):
    x = chebyshev_points(N)
    c = np.ones((N + 1,))
    c[0] = 2
    c[-1] = 2
    c = c * (-1) ** np.arange(N + 1)

    X = np.tile(x, (N + 1, 1))
    dX = X - X.T

    D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
    D = D - np.diag(np.sum(D, axis=1))

    return D


def solve_diff_eq(N):
    x = chebyshev_points(N)
    D = chebyshev_differentiation_matrix(N)
    D2 = D @ D

    # Right-hand side of the equation
    f = (1 / 8) * (x + 1)

    # Matrix A representing the equation
    A = 0.01 * D2 - 0.5 * np.eye(N + 1)

    # Apply boundary conditions: U(-1) = 0, U(1) = 0
    A[0, :] = 0
    A[-1, :] = 0
    A[0, 0] = 1
    A[-1, -1] = 1
    f[0] = 0
    f[-1] = 0

    # Solve the linear system Au = f
    U = solve(A, f)

    return x, U


# Parameters
N = 10  # Number of Chebyshev points

# Solve the differential equation
x, U = solve_diff_eq(N)

# Print the solution
print("Chebyshev points:", x)
print("Solution at Chebyshev points:", U)

import matplotlib.pyplot as plt

# Plot the solution
plt.plot(x, U, 'o-', label='Numerical solution')
plt.xlabel('y')
plt.ylabel('U(y)')
plt.title('Solution of the differential equation using Chebyshev spectral method')
plt.legend()
plt.grid(True)
plt.show()

