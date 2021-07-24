"""This module provides the functions for the analysis of wave speed, its
dependance on frequency and the regression to the environment conditions
"""
import numpy as np
import pandas as pd
# from matplotlib import pyplot as plt
# import seaborn as sns
from sklearn.neighbors import KNeighborsRegressor
from env import Environment

def generate_fingerprint(fingerprint_length: int = 10,#pylint: disable=R0914
                         humidity_n_samples : int = 21,
                         temperature_n_samples : int = 21):
    """
        This function generates a reference database, from a set of
    humidity_n_samples X temperature_n_samples environments.
    For each environment the function simulates a frequency sweep, and searches
    the frequencies f_i where the difference Δc = c(f_i)-c(10Hz) reaches a set
    of threshold values, called delta_thresholds.
    These frequencies, alongside with c(10Hz), constitute a fingerprint of the
    environment state.

    Parameters
    ----------
    humidity_n_samples : int, optional
        Dimension of the sampling of humidity values, in a range between
        0% and 100%. The default is 21.
    temperature_n_samples : int, optional
        Dimension of the sampling of temperature values, in a range between
        0°C and 40°C. The default is 21.
    fingerprint_length : int, optional
        Number of Δc values to take as reference, in a range between
        5mm/s and 75mm/s. The default is 10.

    Returns
    -------
    database : pd.Dataframe
        Table of fingerprints for the set of studied environments. It stores
        the temperature and humidity of the environment, the 0Hz sound speed
        and the frequencies f_i.

    """
    humidity_min = 0
    humidity_max = 100
    humidities = np.linspace(humidity_min, humidity_max,
                             humidity_n_samples)
    temperature_min = 273.15
    temperature_max = 313.15
    temperatures = np.linspace(temperature_min, temperature_max,
                             temperature_n_samples)
    frequency_min = 20
    frequency_max = 1E6
    frequency_n_samples = 1000
    sweep = np.geomspace(frequency_min, frequency_max, frequency_n_samples)
    delta_speed_min = 8E-3
    delta_speed_max = 75E-3
    delta_thresholds = np.linspace(delta_speed_min, delta_speed_max,
                                 fingerprint_length)
    data = []
    for h_i in humidities:
        for t_i in temperatures:
            room = Environment(t_i, h_i)
            speed_varying_f = room.sound_speed_f(sweep)
            speed_10_f = speed_varying_f[0]
            delta_speed = speed_varying_f-speed_10_f
            fingerprint = {dt:sweep[np.nonzero(delta_speed>dt)[0]][0]
                           for dt in delta_thresholds
                           if sweep[np.nonzero(delta_speed>dt)[0]].size>0}
            fingerprint['Temperature'] = t_i - 273.15
            fingerprint['Humidity'] = h_i
            fingerprint['Sound_speed_10'] = speed_10_f
            data.append(fingerprint)
    database = pd.DataFrame(data)
    database = database[['Temperature','Humidity', 'Sound_speed_10',
                        *delta_thresholds]]
    database = database.fillna(0)
    return database

def knn_regressor(database : pd.DataFrame, sample_to_fit : np.ndarray):
    """
    Employs a multidimensional KNN (Kernel Nearest-Neighbour) regressor to the
    simulated enviroments, to fit the model to a sample fingerprint, and
    infer its temperature and humidity.

    Parameters
    ----------
    database : pd.DataFrame
        Table of simulated environments and related fingerprints
    sample_to_fit : np.ndarray
        Fingerprint of an environment, obtained through generate_fingerprint()

    Returns
    -------
    environment_conditions : np.ndarray
        array of the [temperature, humidity] values yielded by the regression

    """
    sample_to_fit = sample_to_fit.fillna(0)
    x_vec = database.drop(['Temperature','Humidity'],axis=1)
    y_vec = database[['Temperature','Humidity']]
    def weight_gauss(dist, sig=2.0):
        return np.exp(-dist**2/(2*sig**2))
    neigh = KNeighborsRegressor(n_neighbors=6, weights=weight_gauss)
    neigh.fit(x_vec, y_vec)
    environment_conditions = neigh.predict(sample_to_fit)
    return environment_conditions

# sns.kdeplot(data=g,x='Sound_speed_10', hue='Temperature',
#             palette='icefire', shade=True, legend=False)
# sns.pairplot(data=g, hue='Temperature', palette='icefire',
#              y_vars=('Humidity','Sound_speed_10'))
# sns.scatterplot(data=g,x='Temperature',y=g.columns[4],hue='Humidity')
# data = []
# delta_thresholds = np.arange(5E-3,75E-3,4E-3)
# humidities = np.arange(0,100,5)
# sweep = np.geomspace(10,1E6,1000)
# temperatures = np.arange(273.15,300,4)
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
