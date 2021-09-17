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

def simulation(signal:np.ndarray, sampling_f:int = 48_000, distance:float = 10,
               temperature:float = 25, humidity:float = 10):
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
        Temperature of the room, in Celsius degrees. The default is 25.
    humidity : float, optional
        Relative humidity of the room, in percentage. The default is 10.

    Returns
    -------
    mic_time : np.ndarray
       Array of the time sequence of the acquired signal.
    mic_signal : np.ndarray
       Array of the intensities of the acquired signal.

    """
    room_dim = [1, 4+distance]
    room = pra.ShoeBox(room_dim, fs=sampling_f, air_absorption=True,
                       materials=pra.Material(1., 0.15), max_order=0,
                       temperature = temperature, humidity = humidity)
    room.set_ray_tracing(receiver_radius=0.5, energy_thres=1e-5,
                         time_thres=13, hist_bin_size=0.002)
    mic_position = np.array([0.5,2])
    distance = [0, distance]
    room.add_microphone(mic_position)
    room.add_source(mic_position+distance, signal=signal)
    room.simulate()
    mic_signal = room.mic_array.signals[0]
    mic_time = np.arange(len(mic_signal)) / sampling_f
    return (mic_time, mic_signal)

def nkt_algorithm(signal:np.ndarray, sampling_f:int=48_000, beta:int = 1500,#pylint: disable=R0914
                  window_type = 'hamming', f_range:list = (24_000, 20)):
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
    tm0 = np.argwhere(abs(signal)>1E-1)[0][0] + 300
    freqs, times, spec_2d = nkt_algorithm(signal[tm0:], sampling_f, [10240,20])
    time_list = times + tm0/sampling_f
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

def measure(distance:float=4000, period:float=10, sampling_f:int=153_600,#pylint: disable=R0914
            method = 'simulation', draw:bool = False):
    """
    Produce a varying-frequency signal and study the different speeds for the
    frequencies which compose it. The measurement can be performed in a
    simulation, with the library pyroomacoustics, or in a real-life experiment,
    which would require a high quality setup of speakers and microphones.

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
    method : { 'simulation', 'experiment' }, optional
        Selects between a measurement conducted in a simulation, or a real-life
        experiment. The default is 'simulation'.
    draw : bool, optional
        If True, it plots the produced and the acquired signals in time and
        frequency domains. The default is False.

    Returns
    -------
    speed_spectrum : 2D np.ndarray
        It contains (frequency, speed) arrays for the studied environment.

    """
    gain_modulation = 10+0.05*(distance-20)
    signal = produce_signal(sampling_f, period, gain_modulation)
    if method == 'simulation':
        temperature = 25
        humidity = 10
        microphone = simulation(signal[1], sampling_f, distance,
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
        raise ValueError('Choose between simulation or experiment')
    spectrum_emitted, freq_emitted = signal_processing(signal[1], sampling_f)
    spectrum_acquired,  = signal_processing(microphone[1], sampling_f)# pylint: disable=unbalanced-tuple-unpacking
    if draw:
        time_plot(signal, microphone)
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
