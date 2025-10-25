import matplotlib.pyplot as plt

x = [1,2,3,4,5,6,1,2,3,4,5,6]
y = [1,2,3,4,5,6,1,2,3,4,5,6]

c = [1,1,1,1,1,1,1,1,1,1,1,1]

plt.bar(x, y)
plt.errorbar(x, y, yerr=c, fmt="o", color="r")
plt.show()