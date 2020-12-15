import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from ProcessClass_10 import ProcessSignal
import matplotlib.pyplot as plt

def fft_amplitude_new(t, u):
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

    # ind_cutoff = p_freq_2 <= cutoff_frequency
    # cut_freq = p_freq_2[ind_cutoff]
    # cut_fft_amp = amp_fft_2[ind_cutoff]
    fft_amp_new_dict = {'frequency': p_freq,
                        'amplitude': amp_fft}
    return fft_amp_new_dict

def fft_amplitude_old(t, u):
    dt = np.abs(t[1] - t[0])
    len_t = len(t)
    u_fft = rfft(u)
    freq_fft = rfftfreq(len_t, dt)
    n = len(u)
    p_freq = freq_fft[1:n:2]
    n_freq = len(p_freq)
    amp_fft = np.zeros(n_freq)
    if n % 2 == 0:
        amp_fft[n_freq - 1] = u_fft[-1] / n
        j = 1
        for i in range(0, n_freq - 1):
            amp_fft[i] = np.sqrt(u_fft[j] ** 2 + u_fft[j + 1] ** 2) / n
            j += 2
    if n % 2 != 0:
        j = 1
        for i in range(0, n_freq - 1):
            amp_fft[i] = np.sqrt(u_fft[j] ** 2 + u_fft[j + 1] ** 2) / n
            j += 2

    ind_fake_freqs_1 = np.logical_or(p_freq <= 1.249e9, 1.251e9 <= p_freq)
    p_freq_1 = p_freq[ind_fake_freqs_1]
    amp_fft_1 = amp_fft[ind_fake_freqs_1]

    ind_fake_freqs_2 = np.logical_or(p_freq_1 <= 2.499e9, 2.501e9 <= p_freq_1)
    p_freq_2 = p_freq_1[ind_fake_freqs_2]
    amp_fft_2 = amp_fft_1[ind_fake_freqs_2]

    # ind_cutoff = p_freq_2 <= cutoff_frequency
    # cut_freq = p_freq_2[ind_cutoff]
    # cut_fft_amp = amp_fft_2[ind_cutoff]
    fft_amp_old_dict = {'frequency': p_freq,
                        'amplitude': amp_fft}
    return fft_amp_old_dict


test = ProcessSignal('201026')
csv_types = test.read_type_file()
csv_signals = csv_types['signal_files'][1::]
csv_signal_nums = csv_types['signal_nums']
excel_dicts = test.read_excel(csv_signal_nums)['numbers']
magnetron_nums = excel_dicts['magnetron']
reb_nums = excel_dicts['noise']
for csv_signal in csv_signals:
    num = csv_signal[3:6]
    file = test.open_file(csv_signal, reduced=False)
    t, u = file['time'], file['voltage']
    cut_file_inds = np.logical_and(t > 70e-9, t < 332e-9)
    t_cut, u_cut = t[cut_file_inds], u[cut_file_inds]
    simple_cut_integral = np.round(test.e_square(t_cut, u_cut) / 1e-8, 2)
    pl_density = test.read_excel(csv_signal)['dicts'][num]['Ток плазмы, А']
    #plt.plot(t, u)
    #plt.title(f'{u[0], u[-1], t[-1]-t[0]}')
    #plt.show()
    dt_2 = (t[-1] - t[0])**2
    fft_old = fft_amplitude_old(t, u)
    amps_old, freqs_old = fft_old['amplitude'], fft_old['frequency']
    fft_new = fft_amplitude_new(t, u)
    amps_new, freqs_new = fft_new['amplitude'], fft_new['frequency']
    integral_u = np.round(test.e_square(t, u) / 1e-8, 2)
    fft_signal = test.fft_signal(t, u, np.abs(t[1]-t[0]))
    simple_amps_fft, simple_freqs_fft = fft_signal['amplitude'], fft_signal['frequency']
    simple_fft_integral = np.round(dt_2 * test.e_square(simple_freqs_fft, simple_amps_fft) / 2e-8, 2)
    integral_amp_old = np.round(dt_2 * test.e_square(freqs_old, amps_old) / 2e-8, 2)
    integral_amp_new = np.round(dt_2 * test.e_square(freqs_new, amps_new) / 2e-8, 2)
    integral_amp = test.e_square(freqs_new, amps_new)
    plt.plot(freqs_old, amps_old, color='k')
    plt.plot(freqs_new, amps_new, color='red')
    plt.xlim(left=2.7e9, right=2.78e9)
    plt.title(f'num = {num}, e_2 = {integral_u}, amp_2 = {simple_cut_integral}, peak={np.max(amps_new)}')
    plt.show()
    print(f'pl_d = {pl_density}, simple_cut_int = {simple_cut_integral}, e_2 = {integral_u},amp_2 = {integral_amp_new}, simple_fft_int = {simple_fft_integral}')