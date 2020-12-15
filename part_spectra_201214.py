import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
import numpy as np
import csv
from scipy.stats import linregress
import openpyxl as xl
import matplotlib.gridspec as gridspec
from scipy.fftpack import rfft, rfftfreq, irfft

folder_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\201124')

proc = ProcessSignal('201124')
csv_types = proc.read_type_file()
csv_signals = csv_types['signal_files']
csv_signal_nums = csv_types['signal_nums']
excel_dicts = proc.read_excel(csv_signal_nums)['numbers']
noise_nums = excel_dicts['noise']
magnetron_nums = excel_dicts['magnetron']
list_3_4 = ['str013.csv', 'str039.csv']
list_1_4 = ['str043.csv', 'str071.csv']

files = os.listdir(folder_path)
magnetron_files = [files[i] for i in range(len(files)) if files[i][3:6] in magnetron_nums and 'csv' in files[i]]
noise_files = [files[i] for i in range(len(files)) if files[i][3:6] in noise_nums and 'csv' in files[i]]


def open_file(file_name):
    if 'csv' in file_name and '_1314' not in file_name:
        file_path = folder_path / f'{file_name}'
        t = []
        u = []
        with open(file_path) as File:
            reader = csv.reader(File)
            for row in reader:
                t_val = float(row[3])
                t.append(t_val)
                u_val = float(row[4])
                u.append(u_val)
        t = np.asarray(t)
        u = np.asarray(u)
        file_dict = {'time': t,
                     'voltage': u}
        return file_dict
    else:
        pass


def fft_amplitude(t, u):
    dt = np.abs(t[1] - t[0])
    len_t = len(t)
    u_fft = rfft(u)
    freq_fft = rfftfreq(len_t, dt)
    n = len(u)
    p_freq = freq_fft[1:n:2]
    n_freq = len(p_freq)
    amp_fft = np.zeros(n_freq)
    amp_fft[0] = u_fft[0] / n
    if n % 2 == 0:
        amp_fft[n_freq - 1] = u_fft[-1] / n
        j = 1
        for i in range(1, n_freq - 2):
            amp_fft[i] = 2 * np.sqrt(u_fft[j] ** 2 + u_fft[j + 1] ** 2) / n
            j += 2
    if n % 2 != 0:
        j = 1
        for i in range(1, n_freq - 1):
            amp_fft[i] = 2 * np.sqrt(u_fft[j] ** 2 + u_fft[j + 1] ** 2) / n
            j += 2

    ind_fake_freqs_1 = np.logical_or(p_freq <= 1.249e9, 1.251e9 <= p_freq)
    p_freq_1 = p_freq[ind_fake_freqs_1]
    amp_fft_1 = amp_fft[ind_fake_freqs_1]

    ind_fake_freqs_2 = np.logical_or(p_freq_1 <= 2.499e9, 2.501e9 <= p_freq_1)
    p_freq_2 = p_freq_1[ind_fake_freqs_2]
    amp_fft_2 = amp_fft_1[ind_fake_freqs_2]

    fft_amp_dict = {'frequency': p_freq_2,
                    'amplitude': amp_fft_2}
    return fft_amp_dict


def mean_frequency(freq, magnitude):
    m_sqre = magnitude * magnitude
    d_freq = freq[1] - freq[0]
    integ_arr = np.zeros(freq.size - 1)
    freq_arr = np.zeros(freq.size - 1)
    for i in range(1, freq.size):
        ind = freq <= freq[i]
        freqs = freq[ind]
        freq_el = freqs[-1]
        m_sqre_el = m_sqre[ind]
        integ_arr_el = np.trapz(m_sqre_el, x=freqs, dx=d_freq)
        integ_arr[i-1] = integ_arr_el
        freq_arr[i-1] = freq_el
    half_integral = 0.5 * np.max(integ_arr)
    try:
        ind_mean_left = integ_arr <= half_integral
        mean_freq_left = freq_arr[ind_mean_left][-1]
        int_left = integ_arr[ind_mean_left][-1]

        ind_mean_right = integ_arr >= half_integral
        mean_freq_right = freq_arr[ind_mean_right][0]
        int_right = integ_arr[ind_mean_right][0]

        if mean_freq_left != mean_freq_right:
            x = np.array([mean_freq_left, mean_freq_right])
            y = np.array([int_left, int_right])
            k, b, r_value, p_value, std_err = linregress(x, y)
            freq_interp = np.linspace(mean_freq_right, mean_freq_left, num=200)
            square_interp = k * freq_interp + b
            linreg_inds = square_interp <= half_integral
            linreg_mean_freq = freq_interp[linreg_inds][0]
        #plt.plot(freq_arr, integ_arr)
        #plt.plot(linreg_mean_freq, square_interp[linreg_inds][0], marker='.', color='red', ms=0.5, linestyle='')
        #plt.hlines(y=half_integral, xmin=2.72e9, xmax=2.76e9, linewidth=0.7)
        #plt.vlines(x=linreg_mean_freq, ymin=0, ymax=square_interp[0], linewidth=0.7)
        #plt.show()
        mean_freq = np.round((linreg_mean_freq / 1e9), 3)
        return mean_freq
    except IndexError:
        pass


def part_mean_freq(file_name):
    file_data = open_file(file_name)
    t, u = file_data['time'], file_data['voltage']
    time_bonds = [120e-9, 280e-9, 440e-9, 600e-9]
    gs = gridspec.GridSpec(3, 1)
    fig = plt.figure(num=1, dpi=300, figsize=[11.69, 14.27])
    for i in range(len(time_bonds) - 1):
        print(time_bonds[i], time_bonds[i+1])
        segment_inds = np.logical_and(time_bonds[i] >= t, t <= time_bonds[i+1])
        t_segment = t[segment_inds]
        u_segment = u[segment_inds]

        fft_data = fft_amplitude(t_segment, u_segment)
        freqs, amps = fft_data['frequency'], fft_data['amplitude']

        if file_name in magnetron_files:
            freq_inds = np.logical_or(np.logical_and(freqs >= 1e9, 2.6e9 >= freqs), np.logical_and(freqs >= 2.8e9, freqs <= 4e9))
            filt_freqs, filt_amps = freqs[freq_inds], amps[freq_inds]
        else:
            freq_inds = np.logical_and(freqs >= 1e9, freqs <= 4e9)
            filt_freqs, filt_amps = freqs[freq_inds], amps[freq_inds]

        mean_freq = mean_frequency(filt_freqs, filt_amps)

        ax = fig.add_subplot(gs[i, 0])
        line, = ax.plot(freqs / 1e9, amps, color='dodgerblue')
        line.set_label(f'{mean_freq} ГГц')
        ax.legend()
        ax.set_title(f'{np.round((time_bonds[i] / 1e-9), 0)}-{np.round((time_bonds[i+1] / 1e-9), 0)} нс')
        ax.plot(filt_freqs / 1e9, filt_amps, color='orange')
        ax.set_xlim(left=1, right=4.5)
        ax.set_xlabel('f, ГГц')
    small_graphs_path = folder_path / 'Pictures' / f'{file_name}_1314.png'
    fig.savefig(small_graphs_path)
    plt.close(fig)


def all_spectra():
    for file in noise_files:
        part_mean_freq(file)

all_spectra()

#part_mean_freq('str039.csv')