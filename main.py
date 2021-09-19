"""This module provides the functions for the analysis of wave speed, its
dependance on frequency and the regression to the environment conditions
"""
import time
import re
import numpy as np
# import pandas as pd
# from matplotlib import pyplot as plt
# from matplotlib.lines import Line2D
# import seaborn as sns
# from env import Environment
import measure



def main():
    print("""This is temperat-here script, whose scope is to determine temperature and relative humidity of air in an environment, through a measure of sound speed.
All necessary data to perform a test are already into Databases/ and Data/ folders, but it is possible to generate new databases for other types of tests.""")
    go_on = True
    while go_on:
        commands = input("List of basic commands:\n-database: create a new database from user-defined inputs (NOTE: it can take up to 30 minutes for large arrays of data);\n-simulation: generate a virtual environment of constant temperature and humidity defined by user, and perform a simulated test;\n-experiment: perform a real-life experiment (requires professional microphone and speaker, still not tested);\n-quit: quit the script.\n")
        commands_splitted = commands.split(" ")
        if "database" in commands_splitted:
            database(commands_splitted[1:])
        elif "simulation" in commands_splitted:
            simulation(commands_splitted[1:])
        elif "experiment" in commands_splitted:
            experiment(commands_splitted[1:])
        elif "quit" in commands_splitted:
            go_on = False
        else:
            print("Command not recognized.\n")

def find_pattern_in_list(pattern, string_list):
    return list(filter(re.compile(pattern).match, string_list))

def database(args=None):
    if not args:
         print("The database will be stored into 'Databases/' folder, for later uses.")
         print("Provide database data as desired, in the following format:")
         print("D-l(integer)-M(float)-m(float) T(integer) H(integer)")
         print("l,M,m are the length, maximum value and minimum value of delta_thresholds,")
         print(" T and H the numbers of temperature values and humidity values;")
         print("If some arguments are not specified, defaults will replace them.")
         print("The arguments can be directly specified when calling database.")
         args = input(">>")
    deltas = find_pattern_in_list("D",args)
    if deltas:
        datas = deltas[0][6:].split('-') #cut 'delta' and split
        try:
            length = int(find_pattern_in_list("l",datas)[0][1:])
        except:
            length = 9
        try:
            max_delta = float(find_pattern_in_list("M",datas)[0][1:])
        except:
            max_delta=9
        try:
            min_delta = float(find_pattern_in_list("m",datas)[0][1:])
        except:
            min_delta = 0.05
        delta_thresholds = measure.generate_delta_thresholds(length, max_delta,
                                                             min_delta)
    else:
        delta_thresholds = measure.generate_delta_thresholds()
    temp_input = find_pattern_in_list("T",args)
    if temp_input:
        temperature_n_samples = int(temp_input[0][1:])
    else:
        temperature_n_samples = 21
    hum_input = find_pattern_in_list("H",args)
    if hum_input:
        humidity_n_samples = int(hum_input[0][1:])
    else:
        humidity_n_samples = 21
    print("Generating database")
    start = time.time()
    measure.generate_database(delta_thresholds, humidity_n_samples,
                             temperature_n_samples, load_and_save=False,
                             method='theory')
    end = time.time()
    print("Done in {}s".format(end-start))

def simulation(args=None):
    pass

def experiment(args=None):
    pass

if __name__ == "__main__":
    main()
# for h in humidities:
#     # ax.set_title(h)
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
