"""This module provides the functions for the analysis of wave speed, its
dependance on frequency and the regression to the environment conditions
"""
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.neighbors import KNeighborsRegressor
from env import Environment
import measure




for test in np.arange(3):
    start = time.time()
    measure.DRAW=True
    distance = 4000#1000*(test+1)
    period = (1+test)*4
    temperature = 300
    humidity = 100
    speed = measure.measure(distance=distance, period=period, method='pyroom',
                    humidity=humidity, temperature=temperature)
    def moving_average(x, w):
        return np.convolve(x, np.ones(w), 'valid') / w
    plt.figure()
    plt.semilogx(speed[0], speed[1], label='raw')
    for i in range(15,16):
        velocities = moving_average(speed[1], i)
        frequencies = moving_average(speed[0], i)
        plt.semilogx(frequencies, velocities, label=i)
    plt.title('Distance = {0}m, time = {1}s'.format(distance,period))
    plt.legend(title='mobile average')
    plt.show()
    end = time.time()
    print('Executed in {:.2f}s'.format(end-start))

# for h in humidities:
#     # fig, ax = plt.subplots()
#     # fig.set_figheight(6)
#     # ax.set_xscale('log')
#     # ax.set_title(h)
#     colors = sns.color_palette('deep',len(temperatures))
#     for t,c in zip(temperatures,colors):
#         room = Environment(t,h)
#         attenuations = room.attenuation_corrections(sweep)
#         speed = room.sound_speed_f(sweep)
#         delta_speed = speed-speed[0]
#         # plt.plot(sweep,delta_speed,color=c,label=t)
#         # fingerprint = [sweep[np.nonzero(delta_speed>dt)[0]][0]
#         #                 for dt in delta_thresholds
#         #                 if sweep[np.nonzero(delta_speed>dt)[0]].size>0]
#         # data.append({'Temperature':t, 'Humidity':h,'Fingerprint':fingerprint})
#         for dt in delta_thresholds:
#             frequencies_above_td = sweep[np.nonzero(delta_speed>dt)[0]]
#             # frequency_treshold = np.nan
#             if frequencies_above_td.size>0:
#                 frequency_treshold = frequencies_above_td[0]
#                 data.append({'Temperature':t, 'Humidity':h, 'Delta_c':dt,
#                              'Frequency':frequency_treshold})
#         # # delta_c = np.sum(speed[0]/(1/(attenuations*speed[0])-1),1)
#         # print(max(delta_speed))
#         # plt.plot(sweep,attenuations.T[0], color=c, linestyle='dotted')
#         # plt.plot(sweep,attenuations.T[1], color=c, linestyle='dashed')
#         # plt.plot(sweep,attenuations.T[0]+attenuations.T[1], color=c, label=t)
#         plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
#                    title='Temperature (K)')
# Database = pd.DataFrame(data)
# Database.to_csv('Fingerprint_db.txt')
# g = sns.relplot(data=Database, x='Delta_c', y = 'Frequency', hue='Humidity',
#             col='Temperature', col_wrap=4,)
# g.set(yscale="log")
