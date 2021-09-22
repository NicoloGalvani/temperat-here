"""This module provides the functions for the analysis of wave speed, its
dependance on frequency and the regression to the environment conditions
"""
import time
import re
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from env import Environment
import measure

def main():
    print("This is temperat-here script, whose scope is to determine temperature and \n",
          "relative humidity of air in an environment, through a measure of sound speed.\n",
          "All necessary data to perform a test are already into Databases/ and Data/ folders,\n",
          " but it is possible to generate new databases for other types of tests.")
    go_on = True
    while go_on:
        print("\nList of basic commands:\n",
              "-theory: hints to theory and examples of use of env package;\n",
              "-database: create a new database from user-defined inputs",
              "(NOTE: it can take up to 30 minutes for large arrays of data);\n",
              "-simulation: generate a virtual environment of constant temperature",
              "and humidity defined by user, and perform a simulated test;\n",
              "-experiment: perform a real-life experiment (requires large spaces, not tested deeply);\n",
              "-quit: quit the script.")
        commands = input(">>")
        commands_splitted = commands.split(" ")
        if "theory" in commands_splitted:
            theory()
        elif "simulation" in commands_splitted:
            simulation(commands_splitted[1:])
        elif "experiment" in commands_splitted:
            experiment(commands_splitted[1:])
        elif "database" in commands_splitted:
            database(commands_splitted[1:])
        elif "quit" in commands_splitted:
            go_on = False
        else:
            print("Command not recognized.\n")

def find_pattern_in_list(pattern, string_list):
    return list(filter(re.compile(pattern).match, string_list))

def database(args=None):
    if not args:
        print("The database will be stored into 'Databases/' folder, for later uses.\n",
        "Provide database data as desired, in the following format:\n",
        "D-l(integer)-M(float)-m(float) F T(integer) H(integer)\n",
        "If F is specified, the database will use a tones-fingerprint instead of delta-fingerprint,\n",
        "l,M,m are the length, maximum value and minimum value of the fingerprint,\n",
        " T and H the numbers of temperature values and humidity values;\n",
        "If some arguments are not specified, defaults will replace them.\n",
        "The arguments can be directly specified when calling database.\n")
        args = input(">>").split(' ')
    deltas = find_pattern_in_list("D", args)
    tone = bool(find_pattern_in_list("F", args))
    if deltas:
        datas = deltas[0][2:].split('-') #cut 'delta' and split
        try:
            length = int(find_pattern_in_list("l", datas)[0][1:])
        except:
            length = 9
        try:
            max_delta = float(find_pattern_in_list("M", datas)[0][1:])
        except:
            max_delta=9
        try:
            min_delta = float(find_pattern_in_list("m", datas)[0][1:])
        except:
            min_delta = 0.05
        if tone:
            input_array = measure.generate_tones(length, max_delta, min_delta)
        else:
            input_array = measure.generate_delta_thresholds(length, max_delta,
                                                             min_delta)
    else:
        if tone:
            input_array = measure.generate_tones()
        else:
            input_array = measure.generate_delta_thresholds()
    temp_input = find_pattern_in_list("T", args)
    if temp_input:
        temperature_n_samples = int(temp_input[0][1:])
    else:
        temperature_n_samples = 21
    hum_input = find_pattern_in_list("H", args)
    if hum_input:
        humidity_n_samples = int(hum_input[0][1:])
    else:
        humidity_n_samples = 21
    print("Generating database")
    start = time.time()
    measure.generate_database(input_array, humidity_n_samples,
                             temperature_n_samples, load_and_save=True,
                             method='simulation', tone=tone)
    end = time.time()
    print("Done in {}s".format(end-start))

def simulation(args=None):
    if not args:
        print("The simulation can be performed in two ways:\n",
        "A ray-tracing simulation performed through package pyroomacoustics,\n",
        "without obstacles between speaker and microphone or sorrounding walls.\n"
         # "which can be called with command 'pyroom';\n",
         # "-A built-in simulation performed decomposing the signal into its frequencies (slower),\n",
         # "which can be called with command 'decomposed'.\n",#NOT CORRECTLY WORKING
        "Provide the simulation conditions in the following format:\n",
        "pyroom/decomposed draw S-d(float)-p(float)-f(int) T(float)-l(int) H(float)-l(int) D-l(integer)-M(float)-m(float) F\n",
        "if draw is specified, the script will plot wave propagation characteristics,\n",
        "d,p,f are the travel distance, signal period and max studied frequency,\n",
        " T and H the temperature and humidity values in environment,\n",
        "If F is specified, the simulation will use a tones-fingerprint instead of delta-fingerprint,\n",
        "l,M,m are the length, maximum value and minimum value of the fingerprint;\n",
        "If some arguments are not specified, defaults will replace them.\n",
        "The arguments can be directly specified when calling simulation.\n")
        args = input(">>").split(' ')
    if 'decomposed' in args:
        method = 'decomposed'
    else:
        method = 'pyroom'
    tone = bool(find_pattern_in_list("F", args))
    humidity_n_samples = 21
    temperature_n_samples = 21
    print("Generating fingerprint")
    start = time.time()
    #Delta_threshold
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
            if tone:
                max_delta=1000
        try:
            min_delta = float(find_pattern_in_list("m",datas)[0][1:])
        except:
            min_delta = 0.05
            if tone:
                min_delta = 20
        delta_thresholds = measure.generate_delta_thresholds(length, max_delta,
                                                             min_delta)
    else:
        length=9
        max_delta=1000
        min_delta=20
        delta_thresholds = measure.generate_delta_thresholds()
    #Temperature
    temp_input = find_pattern_in_list("T",args)
    if temp_input:
        datas = temp_input[0].split('-')
        temperature = float(datas[0][1:])
        try:
            temperature_n_samples = float(find_pattern_in_list("l",datas)[0][1:])
        except:
            temperature_n_samples = 21
    else:
        temperature = 300
        temperature_n_samples = 21
    #Humidity
    hum_input = find_pattern_in_list("H",args)
    if hum_input:
        datas = hum_input[0].split('-')
        humidity = float(datas[0][1:])
        try:
            humidity_n_samples = float(find_pattern_in_list("l",datas)[0][1:])
        except:
            humidity_n_samples = 21
    else:
        humidity = 50
        humidity_n_samples = 21
    #Simulation
    finger_length = None
    if tone:
        try:
            finger_length = length
        except:
            finger_length = 9
    if 'draw' in args:
        measure.DRAW = True
    simulation = find_pattern_in_list("S", args)
    if simulation:
        datas = simulation[0][2:].split('-') #cut 'delta' and split
        try:
            distance = int(find_pattern_in_list("d",datas)[0][1:])
        except:
            distance = 4000
        try:
            period = float(find_pattern_in_list("p",datas)[0][1:])
        except:
            period = 5
        try:
            max_frequency = float(find_pattern_in_list("f",datas)[0][1:])
        except:
            max_frequency = 1000
        frequencies, velocities = measure.measure(distance, period,
                                                  method = method,
                                                  max_frequency = max_frequency,
                                                  temperature = temperature,
                                                  humidity = humidity,
                                                  finger_length = finger_length)
    else:
        frequencies, velocities = measure.measure(method = method,
                                          temperature = temperature,
                                          humidity = humidity,
                                          finger_length = finger_length)
    measure.DRAW = False
    #Fingerprint
    if tone:
        fingerprint = velocities
        input_array = measure.generate_tones(length, max_delta, min_delta)
    else:
        fingerprint = measure.generate_fingerprint(frequencies, velocities,
                                                   delta_thresholds)
        input_array = delta_thresholds
    print("Fingerprint done, searching for a compatible database (if absent, create new one)")
    #Database
    database = measure.generate_database(input_array, humidity_n_samples,
                             temperature_n_samples, load_and_save=True,
                             method='simulation', tone=tone)
    print("Database found, the results are the following:")
    #Classification
    result = measure.knn_regressor(database, fingerprint)[0]
    temp_accuracy = 1-abs(temperature-result[0]-273.15)/temperature
    hum_accuracy = 1-abs(humidity-result[1])/humidity
    print('Temperature: {0:.2f}°C | Accuracy: {1:.2f}%'.format(result[0], 100*temp_accuracy))
    print('Humidity: {0:.2f}% | Accuracy: {1:.2f}%'.format(result[1], 100*hum_accuracy))
    end = time.time()
    print("Done in {}s".format(end-start))

def experiment(args=None):
    if not args:
        print("The experiment requires to put the microphone and the speaker in an open space\n",
        "or in an anechoic room, placing them at the desired distance.\n",
        "At the present stage, the recommended value is on the order of km.\n",
        "Specify experimental conditions in the following format:\n",
        "draw E-d(float)-p(float)-f(int) T(float)-l(int) H(float)-l(int) D-l(integer)-M(float)-m(float)\n",
        "if draw is specified, the script will plot wave propagation characteristics,\n",
        "d,p,f are the travel distance, signal period and max studied frequency,\n",
        " T and H the temperature and humidity values in environment, to test accuracy,\n",
        "l,M,m are the length, maximum value and minimum value of delta_thresholds;\n",
        "If some arguments are not specified, defaults will replace them.\n",
        "The arguments can be directly specified when calling experiment.\n")
        args = input(">>").split(' ')
    method = 'experiment'
    if 'draw' in args:
        measure.DRAW = True
    humidity_n_samples = 21
    temperature_n_samples = 21
    temp_input = find_pattern_in_list("T",args)
    if temp_input:
        datas = temp_input[0].split('-')
        try:
            temperature = float(datas[0][1:])
        except:
            pass
        try:
            temperature_n_samples = float(find_pattern_in_list("l",datas)[0][1:])
        except:
            pass
    hum_input = find_pattern_in_list("H",args)
    if hum_input:
        datas = hum_input[0].split('-')
        try:
            humidity = float(datas[0][1:])
        except:
            pass
        try:
            humidity_n_samples = float(find_pattern_in_list("l",datas)[0][1:])
        except:
            pass
    else:
        humidity = 50
    experiment = find_pattern_in_list("E", args)
    if experiment:
        datas = experiment[0][2:].split('-') #cut 'delta' and split
        try:
            distance = float(find_pattern_in_list("d",datas)[0][1:])
        except:
            distance = 1000
        try:
            period = float(find_pattern_in_list("p",datas)[0][1:])
        except:
            period = 5
        try:
            max_frequency = float(find_pattern_in_list("f",datas)[0][1:])
        except:
            max_frequency = 1000
        try:
            frequencies, velocities = measure.measure(distance, period,
                                                  method = method,
                                                  max_frequency = max_frequency)
        except:
            raise ValueError("Analysis failed, try to increase the distance.")

    else:
        try:
            frequencies, velocities = measure.measure(method = method)
        except:
            raise ValueError("Analysis failed, try to increase the distance.")
    measure.DRAW = False
    print("Measurement done, generating fingerprint")
    start = time.time()
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
    fingerprint = measure.generate_fingerprint(frequencies, velocities,
                                               delta_thresholds)
    print("Fingerprint done, searching for a compatible database (if absent, create new one)")
    database = measure.generate_database(delta_thresholds, humidity_n_samples,
                             temperature_n_samples, load_and_save=True,
                             method='simulation')
    print("Database found, the results are the following:")
    result = measure.knn_regressor(database, fingerprint)[0]
    if temperature:
        temp_accuracy = 100-100*abs(temperature-result[0]-273.15)/temperature
    else:
        temp_accuracy = 0
    if humidity:
        hum_accuracy = max(100-100*abs(humidity-result[1])/humidity, 0)
    else:
        hum_accuracy = 0
    print('Temperature: {0:.2f}°C | Accuracy: {1:.2f}%'.format(result[0], temp_accuracy))
    print('Humidity: {0:.2f}% | Accuracy: {1:.2f}%'.format(result[1], hum_accuracy))
    end = time.time()
    print("Done in {}s".format(end-start))

def theory():
    print("Temperature, humidity and pressure concur to modify the virial coefficients\n",
          "in molecular many-body interaction, even though pressure effect is negligible.\n",
          "Adiabatic constant is proportional to the ratio between first and second derivative\n",
          "of second virial coefficient, and transfer its dependence to sound speed.\n",
          "On a second level, the sound absorption due to nitrogen and oxygen is frequency dependent,\n",
          "with an inpact on both sound attenuation and sound speed.\n",
          "The following plots, realized using methods of class env.Environment, show these effects.\n")
    humidities = np.arange(0, 101, 5)
    temperatures = np.arange(273.15, 314, 4)
    sweep = np.geomspace(20, 22_050, 1000)
    colors = sns.color_palette('deep', len(temperatures))
    legend_elements = [Line2D([0],[0], markersize=8, color=c, label=t-273.15)
                        for t, c in zip(temperatures, colors)]
    fig, axis = plt.subplots(2, 2, sharex = 'col')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    plt.suptitle("Temperature, humidity and frequency effects on sound")
    axis[0,0].set(ylabel='Adiabatic constant')
    axis[0,0].grid()
    axis[0,1].set(ylabel='Atmospheric attenuation at RH=50% (1/m)')
    axis[0,1].grid()
    axis[1,0].set(ylabel='0-frequency sound speed (m/s)',
                xlabel = 'Relative Humidity (%)')
    axis[1,0].grid()
    axis[1,1].set(ylabel='Sound speed frequency variation at RH=50% (m/s)',
                xlabel = 'Frequency (Hz)')
    axis[1,1].grid()
    axis[1,1].set_xscale('log')
    for temp, color in zip(temperatures,colors):
        rooms = [Environment(temp, hum) for hum in humidities]
        adiabatic_coefficients = [room.gamma_adiabatic() for room in rooms]
        speed_0 = [room.sound_speed_0() for room in rooms]
        attenuations = rooms[11].attenuation_factor(sweep)
        speed = rooms[11].sound_speed_f(sweep)
        delta_speed = speed-speed[0]
        axis[0,0].plot(humidities, adiabatic_coefficients, color=color)
        axis[1,0].plot(humidities, speed_0, color=color)
        axis[0,1].plot(sweep, attenuations, color=color)
        axis[1,1].plot(sweep, delta_speed, color=color)
    plt.legend(handles = legend_elements, bbox_to_anchor = (1.05, 1),
               loc = 'upper left', title = 'Temperature (°C)')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
