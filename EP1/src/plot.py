import numpy as np
from matplotlib import pyplot as plt

V0 = 9
def f(x):
    return np.tan(np.sqrt(x + V0)) - np.sqrt(-x/(x + V0))

if __name__ == "__main__":
    x = np.linspace(-9+1e-3,0,1000)
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.plot(x, f(x))
    plt.title(r"$f(E) \times E$")
    plt.xlim(-9,0)
    plt.ylim(-20,20)
    plt.ylabel(r"$f(E)$")
    plt.xlabel(r"$E$ $\left(\frac{\hbar^2}{2ma^2}\right)$")
    plt.grid()
    plt.savefig("plot.png")
    plt.show()
