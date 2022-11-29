import numpy as np
stats = []
with open("./det/all.txt", "r") as f:
    f.readline()
    for line in f:
        tokens = line.split(",")[1:]
        stats.append((float(tokens[0]) + float(tokens[2]), float(tokens[1]) + float(tokens[3])))

stats = np.array(stats) / 2
print(np.mean(stats, axis=0))
print(np.median(stats, axis=0))
