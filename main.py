"""
Module with the actual simulation

"""

import time
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import seaborn as sns
from uncertainties import ufloat
from uncertainties import umath
from uncertainties import unumpy
import uncertainties
from env import Environment
from test_env import read_kayelaby_speed#, read_kayelaby_attenuation

TEMP_0 = 273.15 #K
ATM_P = 101325 #Pa
CP0 = 1.006E3 #j Kg^-1 K^-1
R_GAS = 8.31446261815324 #m^3 Pa K^-1 mol^-1

t = time.time()
Data = read_kayelaby_speed()
Temperatures = unumpy.uarray(Data['Temperature (°C)'] + TEMP_0, 0.5)
Humidities = unumpy.uarray(Data['Relative Humidity (%)'], 5)
rooms = [Environment(Temperatures[i],
                      Humidities[i],ATM_P)
          for i in range(len(Data))]
S0 = [r.sound_speed() for r in rooms]#-Data['Speed (m/s)']
G0 = [r.gamma_adiabatic() for r in rooms]
B0 = [r.b_mix_function(r.t_input) for r in rooms]
S = unumpy.uarray([s.n for s in S0],[s.s for s in S0])
G = unumpy.uarray([g.n for g in G0],[g.s for g in G0])
B = unumpy.uarray([b.n for b in B0],[b.s for b in B0])
fig,ax = plt.subplots(ncols = 2)

fig.set_figheight(6)
fig.set_figwidth(12)
plt.suptitle("Comparison of expected speed with evaluated speed")
ax[0].tick_params(direction = 'in')

ax[1].set(xlabel = 'Temperature (°C)')
ax[1].tick_params(direction = 'in')
sns.lineplot(x = 'Relative Humidity (%)', y = unumpy.nominal_values(S),
                data= Data, hue = 'Temperature (°C)', ax = ax[0])
sns.scatterplot(x = 'Relative Humidity (%)', y = 'Speed (m/s)',
                data= Data, hue = 'Temperature (°C)', ax = ax[0])
sns.lineplot(hue = 'Relative Humidity (%)', y = unumpy.nominal_values(S),
                data= Data, x = 'Temperature (°C)', ax = ax[1])
sns.scatterplot(hue = 'Relative Humidity (%)', y = 'Speed (m/s)',
                data= Data, x = 'Temperature (°C)', ax = ax[1])
elapsed = time.time() - t
print(elapsed)
