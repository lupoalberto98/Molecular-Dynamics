
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Parameters reading
par = np.loadtxt("parameters.txt")
m = par[0]

# Sequence reading
x = np.loadtxt("even.txt")/m*1.0
y = np.loadtxt("odd.txt")/m*1.0



plt.title("Linear Congruence Generator")
plt.xlabel("Even sequence")
plt.ylabel("Odd sequence")
plt.scatter(x,y)
plt.show()
