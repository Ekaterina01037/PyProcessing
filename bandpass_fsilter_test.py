import numpy as np
import csv
import openpyxl as excel
import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
from scipy.fftpack import rfft, irfft, rfftfreq

path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\fft_test')
origin_filt_file = path / 'str015_filt.csv'
file = path / 'str015.csv'

def open_file(file_path, filtered=False):
    t = []
    u = []
    t_filt = []
    u_filt = []
    with open(str(file_path)) as File:
        reader = csv.reader(File)
        for row in reader:
            t.append(float(row[3]))
    with open(str(file_path)) as File:
        reader = csv.reader(File)
        for row in reader:
            u.append(float(row[4]))
    #file_dict = {'t': t, 'u': u}
    if filtered:
        with open(str(file_path)) as File:
            reader = csv.reader(File)
            for row in reader:
                t_filt.append(float(row[5]))
        with open(str(file_path)) as File:
            reader = csv.reader(File)
            for row in reader:
                u_filt.append(float(row[6]))
    file_dict = {'t': t,
                 'u': u,
                 't_filt': t_filt,
                 'u_filt': u_filt}
    return file_dict


def fft_filter(t, u, low_freq=2.695e9, high_freq=2.725e9):
    len_t = len(t)
    dt = np.abs(t[1] - t[0])
    fft_u = rfft(u)
    freqs_fft = rfftfreq(len_t, dt)
    ind_mask = np.logical_and(low_freq < freqs_fft, freqs_fft < high_freq)
    bandpass_filter = np.zeros(len(fft_u))
    bandpass_filter[ind_mask] = 1
    b_filt_u = irfft(fft_u * bandpass_filter)
    return b_filt_u
'''
file_data = open_file(file)
t, u = file_data['t'], np.asarray(file_data['u'])
fft_u = fft_filter(t, u)
origin_data = open_file(origin_filt_file, filtered=True)
t_filt_or, u_filt_or = origin_data['t_filt'], origin_data['u_filt']
'''
#plt.plot(t, u)
#plt.plot(t, fft_u)
#plt.plot(t_filt_or, u_filt_or)
#plt.show()

proc = ProcessSignal('201111')
file = 'str019.csv'
file_data = proc.open_file(file)
t, u = file_data['time'], file_data['voltage']
fft_filt_u = proc.fft_filter(t, u)
bandpass_u = proc.bandpass_filter(t, u, 2.695e9, 2.725e9)
new_file_path = path / 'new_file_19_1.csv'
new_file = open(str(new_file_path), 'w', newline='')
with new_file:
    writer = csv.writer(new_file)
    for i in range(0, len(t)):
        writer.writerow([0, 0, 0, t[i], u[i], fft_filt_u[i], bandpass_u[i]])