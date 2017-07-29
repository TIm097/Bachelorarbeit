import numpy as np

J = np.genfromtxt('test.txt', unpack = 'True')
print(J)

a = np.array([1+1j,1+2j,1-2j])
b = np.array([2+1j,3,1+4j])
z = np.conj(b).dot(a)
print(z)
