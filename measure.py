# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 14:03:02 2021

@author: nicolog
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.fftpack import rfft, fftfreq
import librosa
import sounddevice as sd
import pyroomacoustics as pra
from sklearn.neighbors import KNeighborsRegressor
from env import Environment

DRAW = False

def time_plot(speaker_signal, mic_signal):
    """
    Time domain plot of (top) the signal emitted from the speaker and (bottom)
    the signal acquired from the microphone.

    Parameters
    ----------
    speaker_signal : array-like
        It contains (time_list, intesity) arrays for the speaker.
    mic_signal : array-like
        It contains (time_list, intesity) arrays for the microphone.

    """
    fig, axis = plt.subplots(2, 1)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    plt.suptitle("Time domain visualization")
    axis[0].plot(speaker_signal[0], speaker_signal[1])
    axis[0].set(title='Emitted signal', xlabel="Time (s)", ylabel="Intensity")
    axis[1].plot(mic_signal[0], mic_signal[1])
    axis[1].set(title='Acquired signal', xlabel="Time (s)", ylabel="Intensity")
    plt.tight_layout()
    plt.show()

def spectra_plot(spectrum_emitted, spectrum_acquired):
    """
    Spectrogram of (top) the signal emitted from the speaker and (bottom)
    the signal acquired from the microphone.

    Parameters
    ----------
    spectrum_emitted : array-like
        It contains (time_list, frequency) arrays for the speaker spectrum.
    spectrum_acquired : array-like
        It contains (time_list, frequency) arrays for the microphone spectrum.

    """
    fig, axis = plt.subplots(2, 1)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    axis[0].semilogy(spectrum_emitted[0], spectrum_emitted[1],
                   linestyle='dotted', color='orange', label='maximum power')
    axis[0].set(title='Speaker', xlabel='Time (s)', ylabel='Frequency (Hz)')
    axis[0].legend()
    axis[0].grid()
    axis[1].semilogy(spectrum_acquired[0], spectrum_acquired[1],
                   linestyle='dotted', color='orange', label='maximum power')
    axis[1].set(title='Microphone', xlabel='Time (s)', ylabel='Frequency (Hz)')
    axis[1].legend()
    axis[1].grid()

def produce_signal(sampling_f:int = 48_000, period:float=1, gain:float=10):
    """
    Production of a up-chirp signal from 20Hz to 10.240kHz, modulated to
    compensate the asymmetric-in-frequency attenuation due to air.

    Parameters
    ----------
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    period : float, optional
        Period of the signal, in seconds. The default is 1.
    gain : float, optional
        Baseline gain of the signal. The default is 10.

    Returns
    -------
    time_series : np.ndarray
       Array of the time sequence of the signal.
    wave_series : np.ndarray
       Array of the intensities of the signal.

    """
    max_frequency = 10_240
    min_frequency = 20
    crescendo = gain*librosa.chirp(min_frequency, max_frequency,
                                   sampling_f, duration = period)
    modulation = np.geomspace(0.5, 50, len(crescendo))
    wave_series = modulation*crescendo
    time_series = np.arange(len(wave_series))/sampling_f
    return (time_series, wave_series)


def pyroom_simulation (signal:np.ndarray, sampling_f:int = 48_000,
                       distance:float = 10, temperature:float = -1,
                       humidity:float = -1):
    """
    Simulation of the propagation of a soundwave from a speaker to a
    microphone, using package pyroomacoustics. To simulate an environment
    free of obstacles, the generated room will be a corridor of
    4mx(4+distance)m, where the microphone and speaker are at the opposite
    sides of the room. The walls are chosen to be phonoabsorbent, to avoid
    effects of echo.

    Unfortunately this package uses a rough formula to detemrine the speed of
    sound, without a clear frequency dependance:
        c = 331.4 + 0.6 * t + 0.0124 * h.
    Combining it with simulations not perfectly reproducible, it follows that
    it will be used as a benchmark for the feature extraction test, without
    predictions on temperature and humidity.

    Parameters
    ----------
    signal : np.ndarray
        The signal which must be emitted from the speaker.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    distance : float
        Distance between speaker and microphone, in meters. The default is 10.
    temperature : float, optional
        Temperature of the room, in Kelvin. If negative (default -1) it is
        randomized in operative range [250, 330].
    humidity : float, optional
        Relative humidity of the room, in percentage. If negative (default -1),
        it is randomized in operative range [0 - 100].

    Returns
    -------
    mic_time : np.ndarray
       Array of the time sequence of the acquired signal.
    mic_signal : np.ndarray
       Array of the intensities of the acquired signal.

    """
    if temperature < 0:
        min_temp = 250
        max_temp = 330
        temperature = min_temp + np.random.rand()*(max_temp - min_temp)
    if humidity < 0:
        max_humidity = 100
        humidity = max_humidity * np.random.rand()
    room_dim = [1, 3 + distance]
    room = pra.ShoeBox(room_dim, fs=sampling_f, air_absorption=True,
                       materials=pra.Material(1., 0.15), max_order=0,
                       temperature = temperature-273.15, humidity = humidity)
    room.set_ray_tracing(receiver_radius=1, energy_thres=1e-5,
                         time_thres=13, hist_bin_size=0.002)
    mic_position = np.array([0.5,1])
    speaker_pos = mic_position+[0, distance+1]
    room.add_microphone(mic_position)
    room.add_source(speaker_pos, signal=signal)
    room.simulate()
    mic_signal = room.mic_array.signals[0]
    mic_time = np.arange(len(mic_signal))/ sampling_f
    return (mic_time, mic_signal)

def nkt_algorithm(signal:np.ndarray, sampling_f:int=48_000, beta:int = 1500,#pylint: disable=R0914
                  window_type = 'hamming', f_range:list = (24_000, 20)):
def corrected_simulation (signal:np.ndarray, sampling_f:int = 48_000,#pylint: disable=R0914
                       distance:float = 10, temperature:float = -1,
                       humidity:float = -1):
    """
    Simulation of the propagation of a soundwave from a speaker to a
    microphone. The original signal is transposed into a spectrogram
    representation, and for each frequency component the travel time is
    evaluated as Δt = distance / speed(φ). A new spectrogram is produced
    summing the contributes for each frequency at each arrival time, and it
    is inverted through griffinlim algorithm to recover a wave signal.
    It provides better results than 'pyroom' method, but griffinlim algorithm
    execution time scales worse:

    Parameters
    ----------
    signal : np.ndarray
        The signal which must be emitted from the speaker.
    frequencies : np.ndarray
        The signal which must be emitted from the speaker.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    distance : float
        Distance between speaker and microphone, in meters. The default is 10.
    temperature : float, optional
        Temperature of the room, in Kelvin. If negative (default -1) it is
        randomized in operative range [250, 330].
    humidity : float, optional
        Relative humidity of the room, in percentage. If negative (default -1),
        it is randomized in operative range [0 - 100].

    Returns
    -------
    mic_time : np.ndarray
       Array of the time sequence of the acquired signal.
    mic_signal : np.ndarray
       Array of the intensities of the acquired signal.

    """
    if temperature < 0:
        min_temp = 250
        max_temp = 330
        temperature = min_temp + np.random.rand()*(max_temp - min_temp)
    if humidity < 0:
        max_humidity = 100
        humidity = max_humidity * np.random.rand()
    room = Environment(temperature, humidity)
    max_frequency = 10_240
    min_frequency = 20
    freqs, times, spec_2d = cqt_algorithm(signal, sampling_f,
                                          (max_frequency, min_frequency))
    attenuation = 10 ** (room.attenuation_factor(freqs)*distance/10
                         + 2*np.log10(distance) + 1.1)
    speed_array = room.sound_speed_f(freqs)
    arrival_times=np.array([[time_0+distance/speed_f for time_0 in times]
                    for speed_f in speed_array], dtype = np.float32)
    t_start = arrival_times.min()
    t_end = arrival_times.max()
    t_step = times[1] - times[0]
    time_period_length = int((t_end - t_start) // t_step)
    final_spec_2d = np.zeros([len(freqs), time_period_length],
                             dtype = np.float32)
    arrival_index = ((arrival_times - t_start) // t_step).astype('int32') - 1
    for index, intensity in enumerate(spec_2d):
        final_spec_2d[index, arrival_index[index]
                          ] += intensity/attenuation[index]
    hop_length = int(0.5*max_frequency//min_frequency)
    mic_signal = librosa.griffinlim_cqt(final_spec_2d, sr = sampling_f,
                              fmin = min_frequency, hop_length = hop_length,
                              bins_per_octave = 40)
    mic_time = np.arange(len(mic_signal))/sampling_f+t_start
    return (mic_time, mic_signal)
    """
    Algorithm developed by Nisar, Khan, Tariq in "An Efficient Adaptive Window
    Size Selection Method for Improving Spectrogram Visualization". After a
    first analysis of the input signal, it produces a spectrogram using "Short
    Time Fourier Transform" or "Constant Q Transform" techniques, choosing the
    most efficient for the type of signal.
    The STFT is performed using scipy implementation, while CQT is performed
    using librosa, since scipy lacks it.

    Parameters
    ----------
    signal : np.ndarray
        The signal to analyze.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    beta : int, optional
        Parameter of the algorithm, used to choose the technique to employ.
        The default is 1500.
    window_type : {'rect', 'hamming', 'blackman'}, optional
        Parameter of the STFT. The default is 'hamming'.
    f_range : list, optional
        A priori maximum and minimum frequencies of the signal, used in CQT to
        determine the number of octaves. The default is [24_000, 20].

    Returns
    -------
    freqs : array-like
        List of frequencies for the spectrogram.
    times : array-like
        List of times for the spectrogram.
    spec_2d : 2D array-like
        Power spectral densities of the spectrogram pixels.

    """
    n_samples = len(signal)
    amplitudes = np.abs(rfft(signal))
    frequencies = fftfreq(n_samples, 1/sampling_f)[:n_samples//2+1]
    mean = (amplitudes*frequencies).sum()
    sigma = np.sqrt(((frequencies-mean)**2).sum()/(len(amplitudes-1)))
    if sigma<= beta:
        window_lobe = {'rect' : 2, 'hamming' : 4, 'blackman' : 6}
        width = int(3*sampling_f*window_lobe[window_type]/mean)
        freqs, times, spec_2d = spectrogram(signal, sampling_f, window_type,
                                            width, width//2)
    else:
        n_octaves = int(np.log2(f_range[1]/f_range[0]))
        bins_per_octave = 40
        n_bins = bins_per_octave*n_octaves
        hop_length = 2**(n_octaves-1)
        spec_2d = np.abs(librosa.cqt(signal, sampling_f, hop_length=hop_length,
                                 fmin = f_range[1], n_bins = n_bins,
                                 bins_per_octave = bins_per_octave,
                                 filter_scale = 0.8))
        freqs = librosa.cqt_frequencies(fmin = 20, n_bins = n_bins,
                                     bins_per_octave = bins_per_octave)
        times = hop_length/sampling_f*np.arange(len(spec_2d[0]))
    return freqs, times, spec_2d

def signal_processing(signal:np.ndarray, sampling_f:int=48_000):
    """
    Extract from the input signal the prevalent frequency in time, using the
    spectrogram generated through nkt_algorithm, and restitutes it along with
    a list of all the studied frequencies.

    Parameters
    ----------
    signal : np.ndarray
        The signal which is to process.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.

    Returns
    -------
    spectrum : array-like
        It contains the spectrum of the signal, as (times, frequencies)
    freqs : np.ndarray
        List of the frequencies studied with the spectrogram.

    """
    try:#trim the silent portion of the signal, if there is any
        tm0 = np.argwhere(abs(signal)>1E-3)[0][0] + 300
    except:#pylint: disable=W0702
        tm0 = 0
    time_list = times + time[tm0]
    main_component = np.array([freqs[np.argmax(spec_2d[:,i])]
                               for i in range(len(times))])
    spectrum = (time_list, main_component)
    return spectrum, freqs

def frequency_speed(spectrum_emitted, spectrum_acquired, freq_emitted,
                    distance:float=10):
    """
    Searches the first appearence in the signals of the frequencies listed in
    freqs, and compares their times to determine the speed of sound for that
    frequency. These data are merged to create the speed spectrum.

    Parameters
    ----------
    spectrum_emitted : array-like
        It contains (times, frequency) arrays for the speaker spectrum.
    spectrum_acquired : array-like
        It contains (times, frequency) arrays for the microphone spectrum.
    freq_emitted : np.ndarray
        List of the frequencies studied with the spectrogram.
    distance : float
        Distance between speaker and microphone, in meters. The default is 10.

    Returns
    -------
    speed_spectrum : 2D np.ndarray
        It contains (frequency, speed) arrays for the studied environment.

    """
    time_list = []
    frequencies = []
    for frequency in freq_emitted:
        try:
            t_emitted = spectrum_emitted[0][np.argwhere(
                                            spectrum_emitted[1]>=frequency)][0]
            t_acquired = spectrum_acquired[0][np.argwhere(
                                           spectrum_acquired[1]>=frequency)][0]
            delta_t = t_acquired[0]-t_emitted[0]
            time_list.append(delta_t)
            frequencies.append(frequency)
        except:# pylint: disable=bare-except
            pass
    speeds = distance/np.array(time_list)
    speed_spectrum = np.array([frequencies[30:], speeds[30:]])
    return speed_spectrum

def measure(distance:float=4000, period:float=20, sampling_f:int=153_600,#pylint: disable=R0913 disable=R0914
            method = 'pyroom', temperature:float = -1, humidity:float = -1):
    """
    Produce a varying-frequency signal and study the different speeds for the
    frequencies which compose it. The measurement can be performed in a
    simulation, with the library pyroomacoustics, or in a real-life experiment,
    which would require a high quality setup of speakers and microphones.
    If DRAW==True, it plots the produced and the acquired signals in time and
        frequency domains. The default is False.

    Parameters
    ----------
    distance : float, optional
        Distance between the speaker and the microphone, to know a priori,
        expressed in meters. The default is 4000.
    period : float, optional
        Period length of the signal, expressed in seconds. The default is 10.
    sampling_f : int, optional
        Sampling frequency of the signal, both in emission and acquisition,
        expressed in Hertz. The default is 153_600.
    method : { 'pyroom', 'corrected', 'experiment' }, optional
        Selects between a measurement conducted in a simulation, done through
        a simple formulation of sound propagation, through pyroom-acoustic
        package, or in a real-life experiment. The default is 'pyroom'.

    Returns
    -------
    speed_spectrum : 2D np.ndarray
        It contains (frequency, speed) arrays for the studied environment.

    """
    gain_modulation = 10+0.05*(distance-20)
    time, signal = produce_signal(sampling_f, period, gain_modulation)
    if method == 'corrected':
        mic_time, mic_signal = corrected_simulation(signal, sampling_f, distance,
                                                    temperature, humidity)
    elif method == 'pyroom':
        mic_time, mic_signal = pyroom_simulation(signal, sampling_f, distance,
                                                 temperature, humidity)
    elif method == 'experiment':
        print('Check that microphone and speaker are at {:.3f}m.'.format(
            distance))
        input("Press Enter to continue...")
        extended_signal = np.concatenate((signal[1], np.zeros(len(signal[1]))))
        mic_signal = sd.playrec(extended_signal, samplerate=sampling_f,
                                channels=1, out=np.zeros([3*len(signal[1]),1]))
        mic_time = np.arange(len(mic_signal))/sampling_f
        microphone = (mic_time, mic_signal)
    else:
        raise ValueError('Choose between corrected, pyroom or experiment')
    spectrum_emitted, freq_emitted = signal_processing(signal, time, sampling_f)
    spectrum_acquired, freq_received = signal_processing(mic_signal,# pylint: disable=unused-variable
                                                         mic_time, sampling_f)
    if DRAW:
        time_plot((time, signal), (mic_time, mic_signal))
        spectra_plot(spectrum_emitted, spectrum_acquired)
    speed_spectrum = frequency_speed(spectrum_emitted, spectrum_acquired,
                                     freq_emitted, distance)
    return speed_spectrum

######Stop measure, start analysis

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
    fingerprint_length : int, optional
        Number of Δc values to take as reference, in a range between
        5mm/s and 75mm/s. The default is 10.
    humidity_n_samples : int, optional
        Dimension of the sampling of humidity values, in a range between
        0% and 100%. The default is 21.
    temperature_n_samples : int, optional
        Dimension of the sampling of temperature values, in a range between
        0°C and 40°C. The default is 21.

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
    frequency_max = 22_500
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
            speed_20_f = speed_varying_f[0]
            delta_speed = speed_varying_f-speed_20_f
            fingerprint = {dt:sweep[np.nonzero(delta_speed>dt)[0]][0]
                            for dt in delta_thresholds
                            if sweep[np.nonzero(delta_speed>dt)[0]].size>0}
            fingerprint['Temperature'] = t_i - 273.15
            fingerprint['Humidity'] = h_i
            fingerprint['Sound_speed_20'] = speed_20_f
            data.append(fingerprint)
    database = pd.DataFrame(data)
    database = database[['Temperature','Humidity', 'Sound_speed_20',
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
