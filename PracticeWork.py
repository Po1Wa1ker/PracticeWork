import numpy as np
import cmath
from scipy.fft import ifft
import pandas as pd

def gaus_charact(y, z, sigma, gamma):
    xi = complex(y, z)
    sig2 = 0.5 * sigma * sigma
    
    temp1 = sig2 * xi**2 * 1j - gamma * xi
    
    return temp1

def gaus_ch_exp(y, z, t, sigma, gamma):
    temp1 = gaus_charact(y, z, sigma, gamma)
    
    temp1 = -t * temp1
    return np.exp(temp1)

m = 10
d = 0.01
# Модель Леви
sigma = 0.2
r = 0.05
gamma = 0
temp = 0j


S = 80
# Параметры опциона
K = 80
T = 0.5

temp = gaus_charact(0, -1, sigma, gamma)
print(f"Характеристическая функция в момент xi=i равна: {temp.imag:.8f} {temp.real:.8f}")
gamma = r - 0.5 * sigma * sigma
print(f"{gamma:.8f}")

M = int(2 ** m)
xi = -np.pi / d
dxi = -2 * xi
z = -2

v_F = np.zeros(M, dtype=np.complex128)
v_x = np.zeros(M, dtype=np.complex128)

sign_in = 1
sign_out = 1
print(f"Количество промежуктов времени: {M}")

# Подчет обратного интеграла
for k in range(M):
    # Вычисление коплексных чисел
    phi = gaus_ch_exp(xi, z, T, sigma, gamma)
    F = sign_out * phi
    
    g = -K / ((complex(xi, z) * complex(xi, z+1)))
    v_F[k] = F * g
    
    xi += dxi
    sign_out = -sign_out

v_x = ifft(v_F)

J_list = []
Price_List = []

for l in range(100):
    j = M // 2 + l
    x = -M * d / 2 + j * d
    price = -v_x[j].real * K * np.exp(-x * z) * np.exp(-r * T) / d
    J_list.append(j)
    Price_List.append(price)
    print(f"Цена опциона в момент S=Kexp({l:.1f}): {price:.8f}")

list_of_tuples = list(zip(J_list, Price_List))    
df = pd.DataFrame(list_of_tuples, columns=['J', 'Option price'])

df.to_csv('out.csv')
