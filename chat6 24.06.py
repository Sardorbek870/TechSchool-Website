import numpy as np
import scipy.linalg as la

# Parameters
epsilon = 0.01
N = 20  # Number of Chebyshev polynomials

# Chebyshev-Gauss-Lobatto points
y = np.cos(np.pi * np.arange(N + 1) / N)

# Function f(y) = 1/8 * (y + 1)
f = 1 / 8 * (y + 1)

# Construct the Chebyshev differentiation matrix
def chebyshev_diff_matrix(N):
    if N == 0:
        return np.array([[0]])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2
    c[-1] = 2
    c = c * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
    D = D - np.diag(np.sum(D, axis=1))
    return D

D = chebyshev_diff_matrix(N)
D2 = np.dot(D, D)  # Second differentiation matrix

# Construct the matrix for the PDE
A = epsilon * D2 - 0.5 * np.eye(N + 1)
b = f

# Apply boundary conditions
A[0, :] = 0
A[0, 0] = 1
A[-1, :] = 0
A[-1, -1] = 1
b[0] = 0
b[-1] = 0

# Solve the linear system
a = np.linalg.solve(A, b)

# Construct the solution U(y)
U = np.polynomial.chebyshev.chebval(y, a)

# Output the result
print("Chebyshev coefficients:", a)
print("Solution U(y) at Chebyshev points:", U)

# Compare with the exact solution
import matplotlib.pyplot as plt

def exact_solution(y, epsilon):
    return ((epsilon - 0.5) / (1 - np.exp(-1 / epsilon))) * (1 - np.exp(-(y + 1) / (2 * epsilon))) - (epsilon * (y + 1)) / 2 + (y + 1)**2 / 8

y_fine = np.linspace(-1, 1, 400)
U_exact = exact_solution(y_fine, epsilon)

plt.plot(y_fine, U_exact, label='Exact Solution')
plt.plot(y, U, 'o', label='Chebyshev Solution')
plt.legend()
plt.xlabel('y')
plt.ylabel('U(y)')
plt.title('Solution of the PDE using Chebyshev Spectral Method')
plt.show()
