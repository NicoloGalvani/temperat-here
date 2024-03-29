# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 14:03:02 2021

@author: nicolog
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, windows
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

def spectra_plot(spectrum_emitted, spectrum_acquired, tones=None):
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
    axis[0].set(title='Speaker', xlabel='Time (s)', ylabel='Frequency (Hz)')
    axis[0].legend()
    axis[0].grid()
    if isinstance(tones, np.ndarray):
        axis[0].bar(spectrum_emitted[0], spectrum_emitted[1], width=0.4,
                   linestyle='dotted', color='orange', label='maximum power')
        axis[1].bar(spectrum_acquired[0], spectrum_acquired[1], width=0.4,
                   linestyle='dotted', color='orange', label='maximum power')
    else:
        axis[0].semilogy(spectrum_emitted[0], spectrum_emitted[1],
                   linestyle='dotted', color='orange', label='maximum power')
        axis[1].semilogy(spectrum_acquired[0], spectrum_acquired[1],
                   linestyle='dotted', color='orange', label='maximum power')
    axis[1].set(title='Microphone', xlabel='Time (s)', ylabel='Frequency (Hz)')
    axis[1].legend()
    axis[1].grid()
    plt.tight_layout()
    plt.show()

def speed_plot(frequencies, velocities, tones=None):
    """
    Plot of speed of sound VS frequency in the studied environment.

    Parameters
    ----------
    frequencies : array-like
        Vector of frequencies, expressed in Hz.
    velocities : array-like
        Vector of velocities at different frequencies, expressed in m/s.

    """
    fig, axis = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(10)
    axis.set(title='Speed spectrum', ylabel='Speed (m/s)',
             xlabel='Frequency (Hz)')
    if isinstance(tones, np.ndarray):
        axis.bar(frequencies, velocities, label='filtered frequencies')
    else:
        axis.semilogx(frequencies, velocities,
                  label='mobile average over 15 points')
    axis.legend()
    axis.grid()
    plt.tight_layout()
    plt.show()

def produce_signal(sampling_f:int = 22_050, period:float=1,
                   max_frequency:int = 10_240, gain:float=10, tones = None):
    """
    Production of a up-chirp signal from 20Hz to 10.240kHz, modulated to
    compensate the asymmetric-in-frequency attenuation due to air.

    Parameters
    ----------
    sampling_f : int, optional
        Sampling frequency of the signal, in Hz. The default is 48_000.
    period : float, optional
        Period of the signal, in seconds. The default is 1.
    max_frequency: int, optional
        Maximum frequency of the signal produced, in Hz. The default is 10_240.
    gain : float, optional
        Baseline gain of the signal. The default is 10.
    tones : np.ndarray optional
        Pure tones corresponding to the frequencies inspected in following
        analysis. If different from None, a chord signal is produced.

    Returns
    -------
    time_series : np.ndarray
       Array of the time sequence of the signal.
    wave_series : np.ndarray
       Array of the intensities of the signal.

    """
    if isinstance(tones, np.ndarray):
        signal_length = period*sampling_f
        wave = np.zeros(signal_length)
        for index, tone in enumerate(tones):
            peak_width = int(signal_length/len(tones))
            time_shift = int(index*peak_width//5)
            pre = np.zeros(time_shift)
            post = np.zeros(signal_length-time_shift-peak_width)
            peak = librosa.tone(tone, sr = sampling_f, length = peak_width)
            window = windows.cosine(peak_width)
            wave += np.concatenate((pre, peak*window, post))
        signal = gain*wave
    else:
        min_frequency = 20
        signal = gain*librosa.chirp(min_frequency, max_frequency,
                                       sampling_f, duration = period)
    modulation = 1#np.geomspace(0.5, 50, len(signal))
    wave_series = modulation*signal
    time_series = np.arange(len(wave_series))/sampling_f
    return (time_series, wave_series)

def exp_record (signal : np.ndarray, sampling_f : int = 22_050,
                       distance:float = 1000):
    """
    Emit a sound through hardware speaker and record it through the microphone.
    This measurement condition is still under development, since it requires a
    professional hardware, with good performances in all the studied bandwidth,
    along with 'clean' experimental conditions: anechoic room, large distances,
    low environmental noise.

    Parameters
    ----------
    signal : np.ndarray
        The signal which must be emitted from the speaker.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    distance : float
        Distance between speaker and microphone, in meters. The default is 1000.

    Returns
    -------
    mic_time : np.ndarray
       Array of the time sequence of the acquired signal.
    mic_signal : np.ndarray
       Array of the intensities of the acquired signal.

    """

    print('Check that microphone and speaker are at {:.3f}m.'.format(distance))
    input("Press Enter to start acquisition...")
    period = len(signal)/sampling_f
    min_distance_to_avoid_superimposition = period*360
    extended_signal = np.concatenate((signal, np.zeros(len(signal))))
    if distance < min_distance_to_avoid_superimposition:
        mic_signal = sd.playrec(extended_signal, samplerate=sampling_f,
                            channels=1, out=np.zeros([2*len(signal),1]))
    else:
        extended_signal = np.concatenate((extended_signal, np.zeros(2*len(signal))))
        mic_signal = sd.playrec(extended_signal, samplerate=sampling_f,
                            channels=1, out=np.zeros([4*len(signal),1]))
    sd.wait()
    mic_time = np.arange(len(mic_signal))/sampling_f
    print('Signal acquired, start analysis (still under test).')
    return mic_time, mic_signal

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

def decomposed_simulation (signal:np.ndarray, sampling_f:int = 48_000,#pylint: disable=R0914 disable=R0913
                       max_frequency:int = 10_240, distance:float = 10,
                       temperature:float = -1, humidity:float = -1):
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
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    max_frequency: int, optional
        Maximum frequency of the signal produced, in Hz. The default is 10_240.
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
    min_frequency = 20
    freqs, times, spec_2d = nkt_algorithm(signal, sampling_f,
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
    final_spec_2d = np.zeros([len(freqs), time_period_length])
    arrival_index = ((arrival_times - t_start) // t_step).astype('int32') - 1
    for index, intensity in enumerate(spec_2d):
        final_spec_2d[index, arrival_index[index]
                          ] += intensity/attenuation[index]
    # hop_length = 2**(int(np.log2(max_frequency/min_frequency))-1)
    if (freqs[1]-freqs[0])==(freqs[101]-freqs[100]):
        mic_signal = librosa.griffinlim(final_spec_2d, #hop_length=hop_length,
                                        )
    else:
        mic_signal = librosa.griffinlim_cqt(final_spec_2d, sr = sampling_f,
                              fmin = min_frequency, #hop_length = hop_length,
                              bins_per_octave = 40)
    mic_time = np.arange(len(mic_signal))/sampling_f+t_start
    return (mic_time, mic_signal)

def cqt_algorithm(signal:np.ndarray, sampling_f:int=48_000,
                  f_range:list = (24_000, 20), tones=[]):
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
    f_range : list, optional
        A priori maximum and minimum frequencies of the signal, used in CQT to
        determine the number of octaves. The default is [24_000, 20].
    tones : np.ndarray optional
        Pure tones corresponding to the frequencies inspected in following
        analysis. If different from None, it superimposes f_range.

    Returns
    -------
    freqs : array-like
        List of frequencies for the spectrogram.
    times : array-like
        List of times for the spectrogram.
    spec_2d : 2D array-like
        Power spectral densities of the spectrogram pixels.

    """
    if isinstance(tones, np.ndarray):
        f_range = (tones[-1], tones[0])
        bins_per_octave = 30
    else:
        bins_per_octave = 45
    n_octaves = int(np.log2(f_range[0]/f_range[1]))
    n_bins = bins_per_octave*n_octaves
    hop_length = 2**(n_octaves-1)
    spec_2d = np.abs(librosa.cqt(signal, sampling_f, hop_length = hop_length,
                             fmin = f_range[1], n_bins = n_bins,
                             bins_per_octave = bins_per_octave,
                             filter_scale = 0.8))
    freqs = librosa.cqt_frequencies(fmin = 20, n_bins = n_bins,
                                 bins_per_octave = bins_per_octave)
    times = hop_length/sampling_f*np.arange(len(spec_2d[0]))
    return freqs, times, spec_2d

def nkt_algorithm(signal:np.ndarray, sampling_f:int=48_000,#pylint: disable=R0914
                  f_range:list = (24_000, 20), beta:int = 1500,
                  window_type = 'hamming', tones=[]):
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
    tones : np.ndarray optional
        Pure tones corresponding to the frequencies inspected in following
        analysis, gets passed to cqt_algorithm.

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
    amplitudes = np.abs(rfft(signal, n_samples//2+1))
    frequencies = fftfreq(n_samples, 1/sampling_f)[:n_samples//2+1]
    mean = (amplitudes*frequencies).sum()
    sigma = np.sqrt(((frequencies-mean)**2).sum()/(len(amplitudes)-1))
    if sigma<= beta:
        window_lobe = {'rect' : 2, 'hamming' : 4, 'blackman' : 6}
        width = int(3*sampling_f*window_lobe[window_type]/mean)
        freqs, times, spec_2d = spectrogram(signal, sampling_f, window_type,
                                            width, width//2)
    else:
        freqs, times, spec_2d = cqt_algorithm(signal, sampling_f, f_range, tones)
    return freqs, times, spec_2d

def signal_processing(signal:np.ndarray, time:np.ndarray,
                      sampling_f:int = 48_000, max_frequency:int = 10_240,
                      tones = None):
    """
    Extract from the input signal the prevalent frequency in time, using the
    spectrogram generated through nkt_algorithm, and restitutes it along with
    a list of all the studied frequencies.

    Parameters
    ----------
    signal : np.ndarray
        The signal which is to process.
    time : np.ndarray
        The lists of times for the signal.
    sampling_f : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    max_frequency : int, optional
        Maximum frequency of the signal produced, in Hz. The default is 10_240.    
    tones : np.ndarray optional
        Pure tones corresponding to the frequencies inspected in following
        analysis. If different from None, filters only the tones maxima.

    Returns
    -------
    spectrum : array-like
        It contains the spectrum of the signal, as (times, frequencies)
    freqs : np.ndarray
        List of the frequencies studied with the spectrogram.

    """
    try:#trim the silent portion of the signal, if there is any
        tm0 = np.argwhere(abs(signal)>1E-3)[0][0] + int(sampling_f//160)
    except:#pylint: disable=W0702
        tm0 = 0
    min_frequency=20
    freqs, times, spec_2d = nkt_algorithm(signal[tm0:], sampling_f,
                                        f_range=[max_frequency, min_frequency],
                                        tones = tones)
    time_list = times + time[tm0]
    if isinstance(tones, np.ndarray):
        p_time = np.array([])
        p_intensity = np.array([])
        for tone in tones[:-1]:#merge  
            start_index = np.argwhere(freqs>=tone)[0]
            min_index = int(max(start_index-1, 0))
            max_index = int(min(start_index+2, len(freqs)))
            indeces = range(min_index, max_index)
            intensity_merged = np.array([spec_2d[ind] 
                                         for ind in indeces]).sum(axis=0)
            peak_pos = [np.argmax(intensity_merged)]
            p_time = np.append(p_time, time_list[peak_pos])
            p_intensity = np.append(p_intensity, intensity_merged[peak_pos])
        spectrum = (p_time, tones[:-1])
        return spectrum, tones[:-1]
    else:
        main_component = np.array([freqs[np.argmax(spec_2d[:,i])]
                           for i in range(len(times))])
        spectrum = (time_list, main_component)
        return spectrum, freqs

def frequency_speed(spectrum_emitted:np.ndarray, spectrum_acquired:np.ndarray,
                    freq_emitted:np.ndarray, distance:float=10, tones=None):
    """
    Searches the first appearence in the signals of the frequencies listed in
    freqs, and compares their times to determine the speed of sound for that
    frequency. These data are merged to create the speed spectrum.

    Parameters
    ----------
    spectrum_emitted : np.ndarray
        It contains (times, frequency) arrays for the speaker spectrum.
    spectrum_acquired : np.ndarray
        It contains (times, frequency) arrays for the microphone spectrum.
    freq_emitted : np.ndarray
        List of the frequencies studied with the spectrogram.
    distance : float
        Distance between speaker and microphone, in meters. The default is 10.
    tones : np.ndarray optional
        Pure tones corresponding to the frequencies inspected in following
        analysis. If different from None, evaluates the time difference from
        the peaks positions, with a faster procedure.
        
    Returns
    -------
    speed_spectrum : 2D np.ndarray
        It contains (frequency, speed) arrays for the studied environment.

    """
    time_list = []
    frequencies = []
    if isinstance(tones, np.ndarray):
        frequencies = tones[:-1]
        time_list = spectrum_acquired[0] - spectrum_emitted[0]
        speeds = distance/time_list
    else:
        for frequency in freq_emitted:
            try:
                t_emitted = spectrum_emitted[0][np.argwhere(
                                            spectrum_emitted[1]>=frequency)][0]
                t_acquired = spectrum_acquired[0][np.argwhere(
                                            spectrum_acquired[1]>=frequency)][0]
                delta_t = t_acquired[0]-t_emitted[0]
                if delta_t > 0:
                    time_list.append(delta_t)
                    frequencies.append(frequency)
            except:# pylint: disable=bare-except
                pass

        speeds = distance/np.array(time_list)
        def smoothe(vector, n_points_averaging):
            return np.convolve(vector, np.ones(n_points_averaging), 'valid'
                               ) / n_points_averaging
        speeds = smoothe(speeds[30:], 15)
        frequencies = smoothe(frequencies[30:], 15)
    speed_spectrum = np.array([frequencies, speeds])
    return speed_spectrum


def measure(distance:float=1000, period:float=5, #pylint: disable=R0913 disable=R0914
            sampling_f:int=22_050, max_frequency:int = 1000, finger_length:int = None,
            method = 'decomposed', temperature:float = -1, humidity:float = -1):
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
        expressed in Hertz. The default is 22_050.
    finger_length : int, optional
        Length of the fingerprint which will be produced for the analysis.
        If specified, the signal will be limited to this number of frequencies,
        and the processing will filter them out of received signal.
    method : { 'pyroom', 'decomposed', 'experiment' }, optional
        Selects between a measurement conducted in a simulation, done through
        a simple formulation of sound propagation, through pyroom-acoustic
        package, or in a real-life experiment. The default is 'pyroom'.
    temperature : float
        If the method is not 'experiment', it will be passed to the simulation
        as the environment temperature: if negative, it will become a random
        value. The default is -1.
    humidity : float
        If the method is not 'experiment', it will be passed to the simulation
        as the environment humidity: if negative, it will become a random
        value. The default is -1.
    Returns
    -------
    speed_spectrum : 2D np.ndarray
        It contains (frequency, speed) arrays for the studied environment.
        If finger_length is specified, the frequencies and corresponding speed
        are accordingly discretized.

    """
    gain_modulation = 10+0.05*(distance-20)
    if finger_length:
        tones = generate_tones(finger_length, max_frequency)
    else:
        tones = None
    time, signal = produce_signal(sampling_f, period, max_frequency,
                                  gain_modulation, tones)
    if method == 'decomposed':
        mic_time, mic_signal = decomposed_simulation(signal, sampling_f,
                                                    max_frequency, distance,
                                                    temperature, humidity)
    elif method == 'pyroom':
        mic_time, mic_signal = pyroom_simulation(signal, sampling_f, distance,
                                                 temperature, humidity)
    elif method == 'experiment':
        mic_time, mic_signal = exp_record(signal, sampling_f, distance)
    else:
        raise ValueError('Choose between decomposed, pyroom or experiment')
    spectrum_emitted, freq_emitted = signal_processing(signal, time, sampling_f,
                                                       max_frequency, tones)
    spectrum_acquired, freq_received = signal_processing(mic_signal, mic_time,# pylint: disable=unused-variable
                                                    sampling_f, max_frequency, tones)
    
    speed_spectrum = frequency_speed(spectrum_emitted, spectrum_acquired,
                                     freq_emitted, distance, tones)
    if DRAW:
        time_plot((time, signal), (mic_time, mic_signal))
        spectra_plot(spectrum_emitted, spectrum_acquired, tones)
        speed_plot(speed_spectrum[0], speed_spectrum[1], tones)
    return speed_spectrum

######Stop measure, start analysis

def generate_delta_thresholds(length:int = 9, max_delta:float= 9,
                              min_delta:float = 0.05):
    """
    Baseline function to generate an appropriate delta_threshold vector

    Parameters
    ----------
    length : int, optional
        Number of points taken into account. The default is 9.
    max_delta : float, optional
        Maximum variation for sound speed. The default is 9.
    min_delta : float, optional
        Minimum variation for sound speed. The default is 0.05.

    Returns
    -------
    delta_threshold : np.ndarray
        Vector used for fingerprint preparation of a signal.

    """
    return np.linspace(min_delta, max_delta, length)

def generate_tones(length=9, max_frequency:int = 1000, min_frequency:int = 20):
    return np.geomspace(min_frequency, max_frequency, length)

def generate_fingerprint(frequency_data : np.ndarray, speed_data : np.ndarray,
                         delta_thresholds : np.ndarray):
    """
        This function generates a reference database, from a set of
    humidity_n_samples X temperature_n_samples environments.
    For each environment the function simulates a frequency sweep, and searches
    the frequencies f_i where the difference Δc = c(f_i)-c(20Hz) reaches a set
    of threshold values, called delta_thresholds.
    These frequencies, alongside with c(10Hz), constitute a fingerprint of the
    environment state.

    Parameters
    ----------

    frequency_data : np.ndarray
        Array of frequencies of the characteristic curve, expressed in Hz.
        Must have the same length of frequency_data.

    speed_data : np.ndarray
        Array of speeds of the characteristic curve, expressed in m/s.
        Must have the same length of frequency_data.

    delta_thresholds : np.ndarray
        Number of Δc values to take as reference, in a range between
        5mm/s and 75mm/s. The default is 10.

    Returns
    -------
    fingerprint : pd.Dataframe
        Fingerprints of the sample under study. It stores the temperature and
        humidity of the environment, the 0Hz sound speed and the
        frequencies f_i.

    """

    if len(frequency_data)!=len(speed_data):
        print(len(frequency_data), len(speed_data))
        raise ValueError("frequency_data and speed_data must have the same length")
    speed_20_f = speed_data[0]
    delta_speed = speed_data-speed_20_f
    fingerprint = np.array([frequency_data[np.argwhere(abs(delta_speed)>=dt)][0][0]
                    if frequency_data[np.argwhere(abs(delta_speed)>=dt)].size>0
                    else 0 for dt in delta_thresholds])
    fingerprint = np.insert(fingerprint, 0, speed_20_f)
    return fingerprint

def generate_database(input_array : np.ndarray,#pylint: disable=R0914
                          humidity_n_samples : int = 21,
                          temperature_n_samples : int = 21,
                          method = 'simulation', tone : bool = False,
                          load_and_save:bool = True):
    """
        This function generates a reference database, from a set of
    humidity_n_samples X temperature_n_samples environments. Parameter 'tone'
    determines one functioning way among the following:
    -tone == False:
        For each environment the function simulates a frequency sweep, and searches
        the frequencies f_i where the difference Δc = c(f_i)-c(20Hz) reaches a set
        of threshold values, called delta_thresholds, which enters as input.
        These frequencies, alongside with c(20Hz), constitute a fingerprint of the
        environment state.
    -tone == True:
        For each environment the function simulates a combination of pure tones,
        listed by input, and searches for each of them the speed c(f_i): these
        constitute a fingerprint of the environment state.
    

    Parameters
    ----------
    input_array : np.ndarray
        Array of Δc values to take as reference, expressed in m/s. OR
        Array of pure tones composing the signal, expressed in Hz.
    humidity_n_samples : int, optional
        Dimension of the sampling of humidity values, in a range between
        0% and 100%. The default is 21.
    temperature_n_samples : int, optional
        Dimension of the sampling of temperature values, in a range between
        0°C and 40°C. The default is 21.
    method : { 'simulation', 'theory'}, optional
        Selects between a generation of the database through a simulation, done
        with 'decomposed' measure, or through a direct calculation of c(φ) with
        environment class, faster but less accurate. The default is 'theory'.
    tone : bool optional
        Selects between a discrete frequencies produced signal and a continuos
        frequency sweep; determines also the type of fingerprint adopered.
    load_and_save : bool, optional
        If True, the function checks if a database file with the same inputs
        is already in the folder, and loads it; otherwise it evaluates it and
        save it in the folder.
        If False, it evaluates the database as a variable.

    Returns
    -------
    database : pd.Dataframe
        Table of fingerprints for the set of studied environments. It stores
        the temperature and humidity of the environment, the 0Hz sound speed
        and the frequencies f_i.

    """
    if load_and_save:
        prefix = 'sim'
        if method == 'theory':
            prefix = 'the'
        if tone:
            prefix = 'tone_'+prefix
        rows_meta = 'H{0}T{1}'.format(humidity_n_samples,temperature_n_samples)
        columns_meta = 'len{0}max{1}.csv'.format(len(input_array),
                                                 input_array[-1])
        name = prefix+rows_meta+columns_meta
        folder = 'Databases/'
        try:
            database = pd.read_csv(folder+name, sep=',')
            if 'Unnamed: 0' in database:
                database = database.drop('Unnamed: 0', axis=1)
            return database
        except:#pylint: disable=W0702
            print('Database {} not found, creating one.'.format(name))
    humidity_min = 0
    humidity_max = 100
    humidities = np.linspace(humidity_min, humidity_max, humidity_n_samples)
    temperature_min = 273.15
    temperature_max = 313.15
    temperatures = np.linspace(temperature_min, temperature_max,
                               temperature_n_samples)
    if tone:
        finger_length = len(input_array)
    else:
        finger_length = None
    data = []
    for h_i in humidities:
        for t_i in temperatures:
            if method == 'theory':
                frequency_min = 20
                frequency_max = 10000
                frequency_n_samples = 1000
                frequencies = np.geomspace(frequency_min, frequency_max,
                                           frequency_n_samples)
                room = Environment(t_i, h_i)
                speed_varying_f = room.sound_speed_f(frequencies)
            elif method == 'simulation':
                speed_spectrum = measure(temperature=t_i, humidity=h_i,
                                         method='pyroom',
                                         finger_length=finger_length)
                speed_varying_f = speed_spectrum[1]
                frequencies = speed_spectrum[0]
            else:
                raise ValueError('Choose between theory or simulation')
            if tone:
                fingerprint_array = speed_varying_f
                columns = frequencies
            else:
                fingerprint_array = generate_fingerprint(frequencies,
                                                         speed_varying_f,
                                                         input_array)
                fingerprint_array = np.insert(fingerprint_array, 0,
                                              speed_varying_f[0])
                columns = input_array.astype('object')
                columns = np.insert(columns, 0, 'Sound_speed_20')
            fingerprint = {column : fingerprint_array[i]
                           for i, column in enumerate(columns)}
            fingerprint['Temperature'] = t_i - 273.15
            fingerprint['Humidity'] = h_i
            # fingerprint['Sound_speed_20'] = speed_varying_f[0]
            data.append(fingerprint)
    database = pd.DataFrame(data)
    database = database[['Temperature','Humidity', #'Sound_speed_20',
                        *columns]]
    database = database.fillna(0)
    if load_and_save:
        database.to_csv(folder+name)
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
    sample_to_fit = sample_to_fit.reshape(1, -1)
    x_vec = database.drop(['Temperature','Humidity'], axis=1)
    y_vec = database[['Temperature','Humidity']]
    def weight_gauss(dist, sig=2.0):
        return np.exp(-dist**2/(2*sig**2))
    neigh = KNeighborsRegressor(n_neighbors=6, weights=weight_gauss)
    neigh.fit(x_vec, y_vec)
    environment_conditions = neigh.predict(sample_to_fit)
    return environment_conditions
