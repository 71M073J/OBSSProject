import numpy as np
import time
start = time.time()
for i in range(1000000):
    a = (np.ones(100) > 0).astype(np.int32)
end = time.time()

print(end - start)
start = time.time()
for i in range(1000000):
    a = (np.ones(100) > 0) * 1
end = time.time()

print(end - start)