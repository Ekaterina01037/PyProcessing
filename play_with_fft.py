import numpy as np
import openpyxl as excel
import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
from scipy.fftpack import rfft, irfft, rfftfreq


def fft_filter(t, u):
    len_t = len(t)
    dt = np.abs(t[1] - t[0])
    fft_u = rfft(u)
    freqs_fft = rfftfreq(len_t, dt)
    ind_mask = np.logical_and(2.695e9 < freqs_fft, freqs_fft < 2.725e9)
    u_empty = np.zeros(len(fft_u))
    bandpass_filter = np.zeros(len(fft_u))
    bandpass_filter[ind_mask] = 1
    b_filt_u = irfft(fft_u * bandpass_filter)
    u_empty[ind_mask] = fft_u[ind_mask]
    fft_filtered_u = irfft(u_empty)
    #plt.plot(freqs_fft, bandpass_filter)
    plt.plot(t, u)
    plt.plot(t, fft_filtered_u)
    plt.plot(t, b_filt_u)
    #plt.xlim(left=2.675e9, right=2.75e9)
    plt.show()
    return fft_filtered_u


test = ProcessSignal('201008')
csv_files_list = os.listdir(test.csv_files_path)
csv_types = test.read_type_file()
csv_signals = csv_types['signal_files']
csv_signal_nums = csv_types['signal_nums']
excel_dicts = test.read_excel(csv_signal_nums)['numbers']
magnetron_nums = excel_dicts['magnetron']
reb_nums = excel_dicts['noise']
print(reb_nums)
reb_dict = {}
ambigous_nums =[]
dt_list = []
for csv_signal in csv_signals:
    if csv_signal in csv_files_list:
        num = csv_signal[3:6]
        if num in magnetron_nums or num in reb_nums:
            full_file = test.open_file(csv_signal, reduced=False)
            t_full, u_full = full_file['time'], full_file['voltage']
            file = test.open_file(csv_signal, reduced=True)
            t, u = file['time'], file['voltage']
            dt = np.abs(t[1] - t[0])
            pl_density = test.read_excel(csv_signal)['dicts'][num]['Ток плазмы, А']
            full_integral = np.round(test.e_square(t, u) / 1e-8, 2)
            full_file_full_integral = np.round(test.e_square(t_full, u_full) / 1e-8, 2)

            filtered_u = fft_filter(t, u)
            filtered_u_full = fft_filter(t_full, u_full)
            filt_integral = np.round(test.e_square(t, filtered_u) / 1e-8, 3)
            filt_integral_full = np.round(test.e_square(t_full, filtered_u_full) / 1e-8, 2)

            b_filtered_u = test.bandpass_filter(t, u, 2.695e9, 2.725e9)
            b_filtered_u_full = test.bandpass_filter(t_full, u_full, 2.695e9, 2.725e9)
            b_filt_integral = np.round(test.e_square(t, b_filtered_u) / 1e-8, 2)
            b_filt_integral_full = np.round(test.e_square(t_full, b_filtered_u_full) / 1e-8, 2)
            '''
            fft = test.fft_amplitude(t, u, dt)
            freqs = fft['frequency']
            amps = fft['amplitude']
            n_dt_2 = (t[-1] - t[0]) ** 2
            amp_full_integral = np.round(2 * n_dt_2 * test.e_square(freqs, amps) / 1e-8, 2)
            #peak_fft = test.fft_amplitude(t, filtered_u, dt)
            peak_inds = np.logical_and(freqs > 2.695e9, freqs < 2.725e9)
            freq_peak = freqs[peak_inds]
            amp_peak = amps[peak_inds]
            #p_freqs = peak_fft['frequency']
            #p_amps = peak_fft['amplitude']
            #n_dt_2_peak = (t[-1] - t[0]) ** 2
            peak_amp_integral = np.round(2 * n_dt_2 * test.e_square(freq_peak, amp_peak) / 1e-8, 3)

            integral_ratio = np.round(full_integral / amp_full_integral, 2)
            #peak_integral_ratio = np.round(filt_integral / peak_amp_integral, 2)
            print(f'num = {num}, pl_dens={pl_density}, integral_ratio = {peak_amp_integral}')
            if integral_ratio > 1.1:
                ambigous_nums.append(num)
                dt_list.append(n_dt_2)
            #plt.plot(num, integral_ratio, color='red', marker='o')
            plt.plot(pl_density, amp_full_integral, color='blue', marker='o')
            plt.plot(pl_density, peak_amp_integral, color='red', marker='o')
plt.show()
print(ambigous_nums)
print(dt_list)
'''
'''
        if num in magnetron_nums:
            len_t = len(t)
            dt = np.abs(t[1] - t[0])
            fft_u = rfft(u)
            freqs_fft = rfftfreq(len_t, dt)
            print(len(freqs_fft))

            ind_mask = np.logical_and(2.695e9 < freqs_fft, freqs_fft < 2.725e9)
            u_empty = np.zeros(len(ind_mask))
            u_empty[ind_mask] = fft_u[ind_mask]
            print(len(u_empty))
            fft_filtered_u = irfft(u_empty)
            print(max(u), max(fft_filtered_u))
            plt.plot(t, u, color='k')
            plt.plot(t, fft_filtered_u)
            plt.show()
            filtered_integral = np.round(test.e_square(t, fft_filtered_u) / 1e-7, 2)
            print('e_sqare_integral = ', full_integral, '\n',
                  'filtered_integral =', filtered_integral)
        '''