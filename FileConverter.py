import numpy as np
import os
from pathlib import Path

s = r"C:\Users\Public\Documents\190626"
s_path = Path(s)
a = os.listdir(s)# Список содержимого директории
s2 = r"C:\Users\Public\Documents\190626\Converted"

n = len(a)
print(n)
for i in range(n):
    if ".csv" in a[i]:
        file_path = s_path / '{}'.format(a[i])
        print(file_path)
        ar1 = np.loadtxt(str(file_path),  skiprows=6)
        t, u = ar1[:, 0], ar1[:, 1]
        print(t)
        print(u)
        #np.savetxt(s2+r'\{}'.format(a[i]), ar2)