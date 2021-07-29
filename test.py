

import numpy as np


q = np.array((1,5,4,6,8,9))

r = np.asarray(q)

print(q)
print(q[1:])
print(q[1:])

diff = q[1:] - q[0:-1]

diff_l = diff[:-2]
diff_r = diff[2:]
diff_c = diff[1:-1]


print(diff_l)
print(diff_r)
print(diff_c)

r[2:-1] = diff_l / diff_c

print(r)

