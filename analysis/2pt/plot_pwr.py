import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python plot <file_prefix>")
    sys.exit(1)

file_prefix = sys.argv[1]

data = []

filename = file_prefix + ".pwr"
with open(filename, "r") as f:
    counter = 0

    for i in f:
        tokens = i.strip().split()
        if counter < 0:
            pass
        elif counter == 2:
            title = tokens
            n = len(title)
        elif counter > 2:
            if len(tokens) == n:
                data.append([float(j) for j in tokens])
        counter += 1

data = np.array(data)
data = data.transpose()

for i in range(1, len(data)):
    if i in [5,7,8]:
        plt.plot(data[0], data[i], label=title[i], lw=2)

x_range = [min(data[0]), 5000] 
y_range = [min(data[1]), 40]

plt.xlim(x_range)
plt.ylim(y_range)

plt.xlabel("Frequency (cm$^{-1}$)", size="x-large")
plt.ylabel("N", size="x-large")
plt.legend()
plt.tight_layout()
save_filename = file_prefix + ".svg"
plt.savefig(save_filename, transparent=True)
plt.show()

