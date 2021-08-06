# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 14:03:02 2021

@author: nicolog
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.fft import rfft, fftfreq
import librosa
from librosa import display
import sounddevice as sd
import pyroomacoustics as pra

def time_plot(speaker_signal, mic_signal):
    """
    Time domain plot of (top) the signal emitted from the speaker and (bottom)
    the signal acquired from the microphone. 

    Parameters
    ----------
    speaker_signal : array-like
        It contains (times, intesity) arrays for the speaker.
    mic_signal : array-like
        It contains (times, intesity) arrays for the microphone.

    """
    fig, ax = plt.subplots(2, 1)
    plt.suptitle("Time domain visualization")
    ax[0].plot(speaker_signal[0], speaker_signal[1])
    ax[0].set(title='Emitted signal', xlabel="Time (s)", ylabel="Intensity")
    ax[1].plot(mic_signal[0], mic_signal[1])
    ax[1].set(title='Acquired signal', xlabel="Time (s)", ylabel="Intensity")
    plt.tight_layout()
    plt.show()

def spectra_plot(spectrum_emitted, spectrum_acquired):
    """
    Spectrogram of (top) the signal emitted from the speaker and (bottom)
    the signal acquired from the microphone. 

    Parameters
    ----------
    spectrum_emitted : array-like
        It contains (times, frequency) arrays for the speaker spectrum.
    spectrum_acquired : array-like
        It contains (times, frequency) arrays for the microphone spectrum.

    """
    fig, ax = plt.subplots(2, 1)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    ax[0].semilogy(spectrum_emitted[0], spectrum_emitted[1],
                   linestyle='dotted', color='orange', label='maximum power')
    ax[0].set(title='Speaker', xlabel='Time (s)', ylabel='Frequency (Hz)')
    ax[0].legend()
    ax[0].grid()
    ax[1].semilogy(spectrum_acquired[0], spectrum_acquired[1],
                   linestyle='dotted', color='orange', label='maximum power')
    ax[1].set(title='Microphone', xlabel='Time (s)', ylabel='Frequency (Hz)')
    ax[1].legend()
    ax[1].grid()

def produce_signal(fs:int = 48_000, T:float=1, gain:float=10):
    """
    Production of a up-chirp signal from 20Hz to 10.240kHz, modulated to
    compensate the asymmetric-in-frequency attenuation due to air.

    Parameters
    ----------
    fs : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    T : float, optional
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
    crescendo = gain*librosa.chirp(min_frequency, max_frequency, fs, duration=T)
    modulation = np.geomspace(0.5, 50, len(crescendo))
    wave_series = modulation*crescendo
    time_series = np.arange(len(wave_series))/fs
    return (time_series, wave_series)

def simulation(signal:np.ndarray, fs:int = 48_000, distance:float = 10,
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
    fs : int, optional
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
    room = pra.ShoeBox(
        room_dim, fs=fs, materials=pra.Material(1., 0.15), max_order=0,
        air_absorption=True, temperature = temperature, humidity = humidity
        )
    room.set_ray_tracing(receiver_radius=0.5, energy_thres=1e-5,
                         time_thres=13, hist_bin_size=0.002)
    R = np.array([0.5,2])
    distance = [0, distance]
    room.add_microphone(R)
    room.add_source(R+distance, signal=signal)
    room.simulate()
    mic_signal = room.mic_array.signals[0]
    mic_time = np.arange(len(mic_signal)) / fs
    return (mic_time, mic_signal)

def nkt_algorithm(signal:np.ndarray, fs:int=48_000, beta:int = 1500,
                  window_type = 'hamming', max_frequency:int = 24_000,
                  min_frequency:int = 20):
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
    fs : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    beta : int, optional
        Parameter of the algorithm, used to choose the technique to employ.
        The default is 1500.
    window_type : TYPE, optional
        Parameter of the STFT. The default is 'hamming'.
    max_frequency : int, optional
        A priori maximum frequency of the signal, used in CQT to determine the
        number of octaves. The default is 24_000.
    min_frequency : int, optional
        A priori minimum frequency of the signal, used in CQT to determine the
        number of octaves. The default is 20.

    Returns
    -------
    ff : array-like
        List of frequencies for the spectrogram.
    tt : array-like
        List of times for the spectrogram.
    Sxx : 2D array-like
        Power spectral densities of the spectrogram pixels.

    """
    n_samples = len(signal)
    amplitudes = np.abs(rfft(signal, norm='ortho'))
    frequencies = fftfreq(n_samples, 1/fs)[:n_samples//2+1]
    mu = (amplitudes*frequencies).sum()
    sigma = np.sqrt(((frequencies-mu)**2).sum()/(len(amplitudes-1)))
    if sigma<= beta:
        window_lobe = {'rect' : 2, 'hamming' : 4, 'blackman' : 6}
        width = int(3*fs*window_lobe[window_type]/mu)
        ff, tt, Sxx = spectrogram(signal, fs, window_type, width, width//2)
    else:
        n_octaves = int(np.log2(max_frequency/min_frequency))
        bins_per_octave = 40
        n_bins = bins_per_octave*n_octaves
        hop_length = 2**(n_octaves-1)
        Sxx = np.abs(librosa.cqt(signal, fs, hop_length=hop_length,
                                 fmin=min_frequency, n_bins=n_bins,
                                 bins_per_octave=bins_per_octave,
                                 filter_scale=0.8))
        ff = librosa.cqt_frequencies(fmin=20, n_bins=n_bins,
                                     bins_per_octave=bins_per_octave)
        tt = hop_length/fs*np.arange(len(Sxx[0]))
    return ff, tt, Sxx

def signal_processing(signal:np.ndarray, fs:int=48_000, T:float=1):
    """
    Extract from the input signal the prevalent frequency in time, using the
    spectrogram generated through nkt_algorithm, and restitutes it along with
    a list of all the studied frequencies.

    Parameters
    ----------
    signal : np.ndarray
        The signal which is to process.
    fs : int, optional
        Sampling frequency of the signal, in Hertz. The default is 48_000.
    T : float, optional
        Period of the signal, in seconds. The default is 1.

    Returns
    -------
    spectrum : array-like
        It contains the spectrum of the signal, as (times, frequencies)
    ff : np.ndarray
        List of the frequencies studied with the spectrogram.

    """
    tm0 = np.argwhere(abs(signal)>1E-1)[0][0]+300
    ff, tt, Sxx = nkt_algorithm(signal[tm0:], fs, max_frequency=10240)
    times = tt+tm0/fs
    main_component = np.array([ff[np.argmax(Sxx[:,i])]for i in range(len(tt))])
    spectrum = (times, main_component)
    return spectrum, ff

def frequency_speed(spectrum_emitted, spectrum_acquired, freq_emitted,
                    distance:float=10):
    """
    Searches the first appearence in the signals of the frequencies listed in
    ff, and compares their times to determine the speed of sound for that
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
    times = []
    frequencies = []
    for frequency in freq_emitted:
        try:
            t_emitted = spectrum_emitted[0][np.argwhere(
                                            spectrum_emitted[1]>=frequency)][0]
            t_acquired = spectrum_acquired[0][np.argwhere(
                                           spectrum_acquired[1]>=frequency)][0]
            delta_t = t_acquired[0]-t_emitted[0]
            times.append(delta_t)
            frequencies.append(frequency)
        except:
            pass
    speeds = distance/np.array(times)
    speed_spectrum = np.array([frequencies[30:], speeds[30:]])
    return speed_spectrum

def measure(distance:float=4000, T:float=10, fs:int=153_600,
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
    T : float, optional
        Period length of the signal, expressed in seconds. The default is 10.
    fs : int, optional
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
    signal = produce_signal(fs, T, gain_modulation)
    if method == 'simulation':
        temperature = 25
        humidity = 10
        microphone = simulation(signal[1], fs, distance, temperature, humidity)
    elif method == 'experiment':
        print('Check that microphone and speaker are at {:.3f}m.'.format(
            distance))
        input("Press Enter to continue...")
        extended_signal = np.concatenate((signal[1], np.zeros(len(signal[1]))))
        mic_signal = sd.playrec(extended_signal, samplerate=fs, channels=1,
                                out = np.zeros([3*len(signal[1]), 1]) )
        mic_time = np.arange(len(mic_signal))/fs
        microphone = (mic_time, mic_signal)
    else:
        raise ValueError('Choose between simulation or experiment')
    spectrum_emitted, freq_emitted = signal_processing(signal[1], fs, T)
    spectrum_acquired, freq_acquired = signal_processing(microphone[1], fs, T)
    if draw:
        time_plot(signal, microphone)
        spectra_plot(spectrum_emitted, spectrum_acquired)
    speed_spectrum = frequency_speed(spectrum_emitted, spectrum_acquired,
                                     freq_emitted, distance)
    return speed_spectrum
 
for test in np.arange(10):
    start = time.time()
    speed = measure(draw=True)
    plt.figure()
    plt.semilogx(speed[0], speed[1])
    plt.title('Test = {}'.format(test))
    plt.show()
    end = time.time()


########################OLD FUNCTIONS NOT USEFUL ANYMORE#######################
from scipy.signal import get_window

def cqt_algorithm_old(signal, fs, b_zero:float = 10, min_frequency:float = 20,
                   max_frequency:float = 24_000, scaling:float = 1.4142):
     def k_component(width, signal, number_filters):
         arg = -2*np.pi*number_filters/width
         window = get_window('hamming', width)
         series = window*signal[:width]*complex(np.cos(arg),np.sin(arg))
         return series.sum()/width
     f_zero = b_zero/2 + min_frequency
     b_list = [b_zero,]
     f_list = [f_zero,]
     k = 1
     while b_list[-1]/2+f_list[-1] < max_frequency:
         b_list.append(scaling*b_list[k-1])
         f_list.append(f_zero + np.sum(b_list[:k]) + (b_list[k]-b_zero)/2)
         k += 1
     number_filters = len(b_list)
     transform = [k_component(fs/b_k, signal, number_filters) for b_k in b_list]
     return transform

def undersample(signal, fs, fs_target):
    rateo = fs//fs_target
    new_signal = np.array([signal[i] for i in range(len(signal))
                           if i % rateo == 0])
    return new_signal

def find_stop(signal, fs, t0, threshold: float = 2, delta_t: float = 0.04):
    delta_i = int(fs*delta_t)
    under_threshold = np.argwhere(abs(signal[t0:]) < threshold)
    t_no_signal = [under_threshold[i]
                  for i in range(len(under_threshold)-delta_i)
                  if under_threshold[i+delta_i] - under_threshold[i] == delta_i
                  ][0][0] + t0
    return t_no_signal

def make_spectrogram(signal, fs, t_inset : float = 0, is_low : bool = False,
                    fs_reduced : int = 400, nperseg : int = 512):
    ff, tt, Sxx = nkt_algorithm(signal, fs)
    main_component = np.array([ff[np.argmax(Sxx[:,i])]for i in range(len(tt))])
    time_series = tt + t_inset
    spectrum = [time_series, main_component]
    return spectrum

def signal_processing_old(source_signal, mic_signal, fs, fs_reduced, T, delay):
    ts1 = int(T*fs)
    ts2 = ts1 + int(T*delay*fs)
    tm0 = np.argwhere(abs(mic_signal)>0.01)[0][0]
    tm1_s = find_stop(mic_signal, fs, tm0+ts1, 0.5, T*0.04)
    try:
        tm1_e = find_stop(mic_signal, fs, tm1_s+fs//100, 3, T*0.04)
    except:
        tm1_e = tm1_s + fs//200
    tm2 = len(mic_signal) - fs//200
    times = np.array([tm0, tm1_s, tm1_e, tm2])
    spectrum_0 = make_spectrogram(source_signal[:ts1], fs)
    spectrum_1 = make_spectrogram(mic_signal[tm0:tm1_s], fs, tm0/fs)
    spectrum_2 = make_spectrogram(source_signal[ts2:], fs, ts2/fs, is_low=True,
                                        fs_reduced=fs_reduced,
                                        nperseg=64)#*fs_reduced//1000
    spectrum_3 = make_spectrogram(mic_signal[tm1_e:tm2], fs, tm1_s/fs,
                                  is_low=True, fs_reduced=fs_reduced,
                                  nperseg=64)#*fs_reduced//1000
    spectra = (spectrum_0, spectrum_1, spectrum_2, spectrum_3)
    return spectra, times

def spectra_comparison_old(spectra):
    fig, ax = plt.subplots(2, 2)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    plot_spectrogram('High wave speaker', ax[0,0], spectra[0])
    plot_spectrogram('High wave microphone', ax[1,0],spectra[1])
    plot_spectrogram('Low wave speaker', ax[0,1], spectra[2])
    plot_spectrogram('Low wave microphone', ax[1,1], spectra[3])
    plt.show()

def plot_spectrogram(title, axis, signal):
    axis.semilogy(signal[0], signal[1], linestyle='dotted', color='orange',
                  label='maximum power')
    axis.set_title(title)
    axis.set_xlabel('Time (s)')
    axis.set_ylabel('Frequency (Hz)')
    axis.legend()
    axis.grid()
