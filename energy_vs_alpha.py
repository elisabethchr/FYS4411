
import numpy as np
import matplotlib.pyplot as plt

filename = "alphasImp.txt"
f= open(filename,"r")

def readData(f):
    lines = f.readlines()[1:]

    n = len(lines)
    lines2 = [None] * n
    lines3 = [None] * 2 * n
#   data1 = [0] * 2 * n  # mantissa
#  data2 = [0] * 2 * n  # exponent
    #  data3 = [0]*2*n #Product
    col1 = [0] * n
    col2 = [0] * n
    col3 = [0] * n
    for i in range(0, n):
        a = lines[i]
        a = a[:-1]
        lines[i] = a
        lines2[i] = lines[i].split()

    k = 0
    for j in range(0, n):

        col1[j] = float(lines2[j][0])
        col2[j] = float(lines2[j][1])
        col3[j] = float(lines2[j][2])

    return (col1, col2,col3)

    print(n)

    f.close
print("hello")
col1,col2,col3 = readData(f)
#print(col1)
print(col2)
print(col3)

fig0 = plt.figure(0)
plt.title('1 particle, 1 dim, 1e5 steps')
plt.xlabel('\u03B1', size=12)
plt.ylabel('$<E_{loc}>$', size=12)
plt.plot(col1, col2, 'b', label='Numerical')
plt.plot(col1, col3, 'r', label='Analytical')

#ax = fig0.gca()
#ax.set_xticks(np.arange(20, 210, 10))
# ax.set_yticks(np.arange(0, 1., 30))
# plt.scatter(x, y)
plt.grid()
plt.show()
plt.legend()
plt.savefig("alphaImp1P1D1e5.png")
