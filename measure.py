# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 14:03:02 2021

@author: marcosam
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp, spectrogram, get_window
from scipy.fft import rfft, fftfreq
import pyroomacoustics as pra
import sounddevice as sd
import librosa
# from scipy.io import wavfile
# from scipy.signal import fftconvolve
# import IPython


def time_plot(method, speaker_signal, mic_signal, mic_time, fs,
              times):
    plt.figure()
    plt.suptitle(method)
    plt.subplot(2, 1, 1)
    plt.plot(speaker_signal[0], speaker_signal[1])
    plt.title('Original signal')
    plt.xlabel("Time [s]")
    plt.subplot(2, 1, 2)
    time_start = 0#np.argwhere(abs(mic_signal)>1E-9)[0][0]
    time_end = -1#time_start+int(3*fs)
    plt.plot((np.arange(len(mic_signal))/fs),#[time_start:time_end],
             mic_signal#[time_start:time_end]
             )
    # plt.vlines(times/fs, min(mic_signal), max(mic_signal), colors='orange',
    #            linestyles = 'dotted')
    plt.title("Microphone signal")
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()


def plot_spectrogram(title, axis, signal):
    axis.semilogy(signal[0], signal[1], linestyle='dotted', color='orange',
                  label='maximum power')
    axis.set_title(title)
    axis.set_xlabel('time (s)')
    axis.set_ylabel('Frequency (Hz)')
    # if 'High' in title:
    # axis.set_yscale('log')
    axis.legend()
    axis.grid()

def find_stop(signal, fs, t0, threshold: float = 2, delta_t: float = 0.04):
    delta_i = int(fs*delta_t)
    under_threshold = np.argwhere(abs(signal[t0:]) < threshold)
    t_no_signal = [under_threshold[i]
                  for i in range(len(under_threshold)-delta_i)
                  if under_threshold[i+delta_i] - under_threshold[i] == delta_i
                  ][0][0] + t0
    return t_no_signal

def undersample(signal, fs, fs_target):
    rateo = fs//fs_target
    new_signal = np.array([signal[i] for i in range(len(signal))
                           if i % rateo == 0])
    return new_signal

def produce_signal(fs:int = 48_000, T:float=0.75,
                   delay:float = 0.2, gain:float=10):
    # base_gain = 10
    max_frequency = 10_240
    # mid_frequency = 200
    min_frequency = 20
    # n_per_period = int(T*fs)
    crescendo = gain*librosa.chirp(min_frequency, max_frequency, fs,
                                duration=T)
    modulation = np.geomspace(80,1,len(crescendo))
    wave_series = modulation*crescendo#np.flip(crescendo)*modulation
    time_series = np.arange(len(wave_series))/fs
    # tH = np.arange(0, n_per_period) / fs
    # tD = tH.max() + np.arange(0, int(delay*n_per_period))/fs
    # tL = tD.max() + np.arange(0, 3*n_per_period) / fs

    # wave_high = base_gain*chirp(tH, f0=max_frequency, f1=mid_frequency, t1=T,
    #                    method='logarithmic')
    # wave_low= modulation*base_gain*chirp(np.arange(0, 3*n_per_period)/fs,
    #                    f0=mid_frequency, f1=min_frequency, t1=2.98*T,
    #                   method='linear')
    # time_series = np.concatenate((tH, tD, tL))
    # wave_series = np.concatenate((wave_high, np.zeros(tD.shape), wave_low))
    signal = (time_series, wave_series)
    return signal

def simulation(signal, fs, distance, temperature, humidity):
    room_dim = [3, 4+distance]  # meters
    room = pra.ShoeBox(
        room_dim, fs=fs, materials=pra.Material(1., 0.15), max_order=0,
        air_absorption=True, temperature = temperature, humidity = humidity,
        # ray_tracing = True
        )
    # Set the ray tracing parameters
    room.set_ray_tracing(receiver_radius=1.5, energy_thres=1e-5)
    # add source and set the signal to WAV file content
    R = np.array([1.5,2])
    room.add_microphone(R)
    distance = [0, distance]
    room.add_source(R+distance, signal=signal)
    room.image_source_model()
    room.simulate()
    sig = room.mic_array.signals[0]
    mic_time = np.arange(len(sig))/fs
    return sig, mic_time

# def cqt_algorithm(signal, fs, b_zero:float = 10, min_frequency:float = 20,
#                   max_frequency:float = 24_000, scaling:float = 1.4142):
#     def k_component(width, signal, number_filters):
#         arg = -2*np.pi*number_filters/width
#         window = get_window('hamming', width)
#         series = window*signal[:width]*complex(np.cos(arg),np.sin(arg))
#         return series.sum()/width
#     f_zero = b_zero/2 + min_frequency
#     b_list = [b_zero,]
#     f_list = [f_zero,]
#     k = 1
#     while b_list[-1]/2+f_list[-1] < max_frequency:
#         b_list.append(scaling*b_list[k-1])
#         f_list.append(f_zero + np.sum(b_list[:k]) + (b_list[k]-b_zero)/2)
#         k += 1
#     number_filters = len(b_list)
#     transform = [k_component(fs/b_k, signal, number_filters) for b_k in b_list]
#     return transform

def nkt_algorithm(signal, fs, beta = 1500, window_type = 'hamming',
                  max_frequency = 24_000, min_frequency = 20
                  ):
    n_samples = len(signal)
    # print(n_samples)
    amplitudes = np.abs(rfft(signal, norm='ortho'))
    frequencies = fftfreq(n_samples, 1/fs)[:n_samples//2+1]
    mu = (amplitudes*frequencies).sum()
    sigma = np.sqrt(((frequencies-mu)**2).sum()/(len(amplitudes-1)))
    window_lobe = {'rect' : 2, 'hamming' : 4, 'blackman' : 6}
    width = int(3*fs*window_lobe[window_type]/mu)
    if sigma<= beta:
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
        # b_zero = 10
        # cqt_algorithm(signal, fs, b_zero, min_frequency, max_frequency)
    return ff, tt, Sxx

def make_spectrogram(signal, fs, t_inset : float = 0, is_low : bool = False,
                         fs_reduced : int = 400, nperseg : int = 512
                         ):
    # if is_low:
    #     signal = undersample(signal, fs, fs_reduced)
    #     fs = fs_reduced
    ff, tt, Sxx = nkt_algorithm(signal, fs)
    #spectrogram(signal, fs=fs, nperseg=nperseg)
    main_component = np.array([ff[np.argmax(Sxx[:,i])]for i in range(len(tt))])
    time_series = tt + t_inset
    spectrum = [time_series, main_component]
    return spectrum

def signal_processing_divided(source_signal, mic_signal, fs, fs_reduced, T, delay):
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

def spectra_comparison(spectrum_emitted, spectrum_acquired):
    fig, ax = plt.subplots(2, 1)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    plot_spectrogram('Speaker', ax[0], spectrum_emitted)
    plot_spectrogram('Microphone', ax[1],spectrum_acquired)
    # plot_spectrogram('High wave speaker', ax[0,0], spectra[0])
    # plot_spectrogram('High wave microphone', ax[1,0],spectra[1])
    # plot_spectrogram('Low wave speaker', ax[0,1], spectra[2])
    # plot_spectrogram('Low wave microphone', ax[1,1], spectra[3])
    plt.show()

def signal_processing(signal, fs, T):
    tm0 = np.argwhere(abs(signal)>1E-1)[0][0]+370
    # tm1 = [i for i in range(tm0,len(signal)-1,1)
    #        if signal[i+1]-signal[i]>2][0] + tm0
    ff, tt, Sxx = nkt_algorithm(signal[tm0:], fs, max_frequency=10240)
    times = tt+tm0/fs
    main_component = np.array([ff[np.argmax(Sxx[:,i])]for i in range(len(tt))])
    spectrum = [times, main_component]
    return spectrum, ff

def measure(distance:float=10, T:float=0.75, delay:float=0.1, fs:int=48_000,
            fs_reduced:int=8_000, method:str = 'pyroom', draw:bool = False):
    gain_modulation = 10+0.05*(distance-20)
    signal = produce_signal(fs, T, delay, gain_modulation)
    if method == 'pyroom':
        temperature = 15#####tocheck
        humidity = 100#80#####tocheck
        mic_signal, mic_time = simulation(signal[1], fs,
                                          distance, temperature, humidity)
        print('simulation done')
    elif method == 'hardware':
        print('Check that microphone and speaker are at {:.3f}m.'.format(
            distance))
        input("Press Enter to continue...")
        extended_signal = np.concatenate((signal[1], np.zeros(len(signal[1]))))
        mic_signal = sd.playrec(extended_signal, out = np.zeros([2*len(signal[1]),1]),
                                samplerate=fs, channels=1)
        mic_time = np.arange(len(mic_signal))/fs
    else:
        raise ValueError('Choose between pyroom or hardware')
    # spectra, times = signal_processing_old(signal[1], mic_signal, fs,
    #                                    fs_reduced, T, delay)
    spectrum_emitted, freq_emitted = signal_processing(signal[1], fs, T)
    spectrum_acquired, freq_acquired = signal_processing(mic_signal, fs, T)
    times=[]
    if draw:
        time_plot('Pyroomacoustics simulation', signal, mic_signal, mic_time,
                  fs, times)
        spectra_comparison(spectrum_emitted, spectrum_acquired)
    speed_spectrum = frequency_speed(spectrum_emitted, spectrum_acquired,
                                     freq_emitted, distance)
    return spectrum_acquired, speed_spectrum

def frequency_speed(spectrum_emitted, spectrum_acquired, freq_emitted,
                    distance):
    # spectrum_emitted = np.concatenate((spectra[0],spectra[2]),1)
    # spectrum_acquired = np.concatenate((spectra[1],spectra[3]),1)
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
    # plt.figure()
    # plt.plot(times)
    speeds = distance/np.array(times)
    delta_speed = speeds - speeds[0]
    speed_spectrum = np.array([frequencies[30:], speeds[30:]])
    return speed_spectrum

# for multiplier in [2.5,3,3.2,3.4]:
for distance in [4000,]:#np.arange(10,38,1):
    start = time.time()
    T = 10#412/np.sqrt(multiplier)
    fs = int(48_000*3.2)
    fs_reduced = 1_000#*multiplier
    spec, speed = measure(distance, T, fs=fs, fs_reduced=fs_reduced, draw=True)
    plt.semilogx(speed[0], speed[1])
    # plt.ylim(339,345)
    plt.title('mult={0}, dist={1}'.format(fs, distance))
    plt.show()
    end = time.time()
    print(fs, distance, end-start)
# plt.figure()
# plt.plot(spectra[0],spectra[1])
# plt.ylim(0,5)
# plt.xscale('log')
# micro = undersample(mic, 1_000_000, 22_000)
# signal = undersample(sig[1], 1_000_000, 22_000)
# sd.play(signal, 22_000)
# sd.play(micro, 22_000)
###Signal processing
# import scipy.io.wavfile
# scipy.io.wavfile.write('microphone_sound.wav', 22_000, micro)
# scipy.io.wavfile.write('speaker_sound.wav', 22_000, signal)
#View




#######Scrivere test
# variare distance (alti valori), T, H, delay? tenerlo nei limiti
# freq low < 400
# freq high > 150
#%%

# import numpy as np
# from scipy import signal
# # from scipy.fft import fftshift
# import matplotlib.pyplot as plt


# freq_max = 10000#00
# fs = freq_max*2#44100
# N = 600
# n = 1000
# amp = 100#200 * np.sqrt(2)
# noise_power = 0.01 * fs / 2
# # time = np.arange(N) / float(fs)
# dt = 1/fs
# #freq = 440


# x=np.array([0])
# y=np.array([0])
# f_real = np.array([0])

# for freq in np.geomspace(freq_max,20,N):
# #    samples = dt*sampling_rate
#     dx = dt*(1+np.arange(n))+x[-1]#*(1+(freq_max-freq)//freq_max))
#     dy = amp*np.sin(2 * np.pi * freq * dx )#/ sampling_rate
#     x = np.concatenate((x,dx))
#     y = np.concatenate((y,dy))
#     f_real = np.concatenate((f_real,np.array([freq])))
# # plt.figure()
# # plt.plot(x,y)

# # rng = np.random.default_rng()
# # fs = 12e4
# # N = 1e6
# # amp = 200 * np.sqrt(2)
# # noise_power = 0.01 * fs / 2
# # time = np.arange(N) / float(fs)
# # # mod = 80*np.cos(2*np.pi*0.25*time**2)
# # carrier = amp * np.sin(2*np.pi*3e3*time*(1+time))#+mod
# # noise = rng.normal(scale=np.sqrt(noise_power), size=time.shape)
# # noise *= np.exp(-time/5)
# # x = carrier + noise


# plt.figure()
# f, t, Sxx = signal.spectrogram(y, fs, window='nuttall',
#                                nperseg=8000, noverlap=20,)
# peak_frequency = [np.argmax(Sxx[:,i]) for i in range(len(t))]
# plt.scatter(t,f[peak_frequency], label='fft')
# plt.plot(np.linspace(0,max(t),N+1),f_real, color='orange', label='generated')
# # plt.pcolormesh(t, f, Sxx, shading='gouraud')
# plt.ylabel('Frequency [Hz]')
# plt.yscale('log')
# # plt.ylim([10,1000000])
# plt.xlabel('Time [sec]')
# # plt.legend()
# plt.show()

# #%%
# import numpy as np
# from scipy.fftpack import fft
# import time
# import scipy.io.wavfile
# import matplotlib.pyplot as plt
# import sounddevice as sd
# import os

# SAMPLE_FREQ = 44100 # Sampling frequency of the recording
# SAMPLE_DUR = 2  # Duration of the recoding

# print("Grab your guitar!")
# time.sleep(1) # Gives you a second to grab your guitar ;)

# myRecording = sd.rec(SAMPLE_DUR * SAMPLE_FREQ, samplerate=SAMPLE_FREQ,
#                      channels=1,dtype='float64')
# print("Recording audio")
# sd.wait()

# sd.play(myRecording, SAMPLE_FREQ)
# print("Playing audio")
# sd.wait()

# scipy.io.wavfile.write('example1.wav', SAMPLE_FREQ, myRecording)


# sampleFreq, myRecording = scipy.io.wavfile.read("example1.wav")
# sampleDur = len(myRecording)/sampleFreq
# timeX = np.arange(0,sampleDur, 1/sampleFreq)

# plt.plot(timeX, myRecording)
# plt.ylabel('x(k)')
# plt.xlabel('time[s]')
# plt.title('time domain wave')
# plt.show()

# sampleFreq, myRecording = scipy.io.wavfile.read("example1.wav")
# sampleDur = len(myRecording)/sampleFreq
# timeX = np.arange(0, sampleFreq/2, sampleFreq/len(myRecording))
# absFreqSpectrum = abs(fft(myRecording))
# print(absFreqSpectrum)

# plt.plot(timeX, absFreqSpectrum[:len(myRecording)//2])
# plt.ylabel('|X(n)|')
# plt.xlabel('frequency[Hz]')
# plt.title('frequency domain wave')
# plt.show()
