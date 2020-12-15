from ProcessClass import ProcessSignal
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def max_freq(exp):
    test = ProcessSignal(str(exp))
    csv_files = test.csv_files()
    signal_dict = {}
    for csv_file in csv_files:
        csv_num = int(csv_file[3:6])
        if csv_num % 2 == 0 and csv_num > 2:
            file = test.open_file(csv_file)
            t = file['time']
            v = file['voltage']
            dt = file['time_resolution']
            t_0 = t[0]
            if csv_num in [4, 6]:
                t_1 = t_0 + 150e-9
                t_step = 29e-9
            else:
                t_1 = t_0 + 80e-9
                t_step = 40e-9
            part_dict = {}
            for i in range(12):
                t_2 = t_1 + 262e-9
                part_inds = np.logical_and(t > t_1, t < t_2)
                t_part = t[part_inds]
                v_part = v[part_inds]

                fft_dict = test.fft_amplitude(t_part, v_part, dt)
                ampls = fft_dict['amplitude']
                freqs = fft_dict['frequency']
                #peak
                peak = max(ampls)
                max_ind = np.argmax(ampls)
                freq_peak = freqs[max_ind]
                #half_peak_lvl
                half_peak = peak / 2

                f_1_ind_bottom = np.logical_and(ampls <= half_peak, freqs < freq_peak)
                f_1_bottom = freqs[f_1_ind_bottom][-1]
                amp_1_bottom = ampls[f_1_ind_bottom][-1]

                f_1_ind_top = ampls >= half_peak
                f_1_top = freqs[f_1_ind_top][0]
                amp_1_top = ampls[f_1_ind_top][0]

                x_1 = np.array([f_1_bottom, f_1_top])
                y_1 = np.array([amp_1_bottom, amp_1_top])
                k_1, b_1, r_value_1, p_value_1, std_err_1 = linregress(x_1, y_1)
                freq_interp_1 = np.linspace(f_1_bottom, f_1_top, num=200)
                amp_interp_1 = k_1 * freq_interp_1 + b_1
                linreg_inds_1 = amp_interp_1 <= half_peak
                linreg_freq_1 = freq_interp_1[linreg_inds_1][-1]

                f_2_top = freqs[f_1_ind_top][-1]
                amp_2_top = ampls[f_1_ind_top][-1]

                f_2_ind_bottom = np.logical_and(ampls <= half_peak, freqs > freq_peak)
                f_2_bottom = freqs[f_2_ind_bottom][0]
                amp_2_bottom = ampls[f_2_ind_bottom][0]

                x_2 = np.array([f_2_bottom, f_2_top])
                y_2 = np.array([amp_2_bottom, amp_2_top])
                k_2, b_2, r_value_2, p_value_2, std_err_2 = linregress(x_2, y_2)
                freq_interp_2 = np.linspace(f_2_bottom, f_2_top, num=200)
                amp_interp_2 = k_2 * freq_interp_2 + b_2
                linreg_inds_2 = amp_interp_2 <= half_peak
                linreg_freq_2 = freq_interp_2[linreg_inds_2][-1]

                delta_f = linreg_freq_2 - linreg_freq_1
                '''
                plt.plot(freqs, ampls)
                plt.xlim(left=2.73e9, right=2.75e9)
                plt.hlines(half_peak, xmin=0, xmax=4e9)
                plt.hlines(peak, xmin=0, xmax=4e9)
                plt.plot(freq_interp_1, amp_interp_1, marker='.', color='red', ms=0.5)
                plt.vlines(linreg_freq_1, ymin=0, ymax=peak)
                plt.plot(freq_interp_2, amp_interp_2, marker='.', color='green', ms=0.5)
                plt.vlines(linreg_freq_2, ymin=0, ymax=peak)
                plt.title('num = {}, peak = {}, delta = {}, freq = {}'.format(csv_num, np.round(peak, 3), np.round(delta_f / 1e6, 3), np.round(freq_peak / 1e9, 4)))
                plt.show()
                '''
                part_dict[i+1] = {'t_1': np.round((t_1 / 1e-9), 5), 't_2': t_2,
                                  'd_t': t_2 - t_1, 't_step': t_step,
                                  'peak': peak, 'freq_peak': np.round((freq_peak / 1e9), 3),
                                  'delta_f': np.round((delta_f / 1e6), 3)}
                print('t_1 =', np.round((t_1 / 1e-9), 5))

                t_1 = t_1 + t_step
            signal_dict[csv_num] = part_dict
    return signal_dict
#max_freq(191126)

def peak_table(exp, plot_type='peak'):
    test = ProcessSignal(str(exp))
    peak_dict = max_freq(exp)
    keys = peak_dict.keys()
    fig = plt.figure(num=1, dpi=300)
    ax = fig.add_subplot(111)
    for key in keys:
        part_dicts = peak_dict[key]
        part_keys = part_dicts.keys()
        t = []
        peaks = []
        freq_peaks = []
        delta_fs = []
        signal_nums = []
        for part_key in part_keys:
            signal_nums.append(key)
            file = part_dicts[part_key]
            time = file['t_1']
            if time not in t:
                t.append(time)
                peaks.append(file['peak'])
                freq_peaks.append(file['freq_peak'])
                delta_fs.append(file['delta_f'])
        sorted_inds = np.argsort(np.array(t))
        t = np.array(t)
        peaks = np.array(peaks)
        delta_fs = np.array(delta_fs)
        freq_peaks = np.array(freq_peaks)
        signal_nums = np.array(signal_nums)
        t_sorted = t[sorted_inds]
        peak_sorted = peaks[sorted_inds]
        delta_f_sorted = delta_fs[sorted_inds]
        freq_peaks = freq_peaks[sorted_inds]
        signal_nums_sorted = signal_nums[sorted_inds]

        if plot_type is 'peak':
            line, = ax.plot(t_sorted, peak_sorted, marker='^', linewidth=0.9,  ms=3)
            line.set_label('{}'.format(key))
            ax.set_ylabel('Peak amplitude', fontsize=12)
            png_name = test.pics_path / 'Peak'
        if plot_type is 'peak_freq':
            line, = ax.plot(t_sorted, freq_peaks, marker='^', linewidth=0.9,  ms=3)
            line.set_label('{}'.format(key))
            ax.set_ylabel('Peak frequency', fontsize=12)
            ax.set_ylim(bottom=2.735, top=2.745)
            png_name = test.pics_path / 'Peak_frequency'
        if plot_type is 'delta':
            line, = ax.plot(t_sorted, delta_f_sorted, marker='^', linewidth=0.9, ms=3)
            line.set_label('{}'.format(key))
            ax.set_ylabel('Peak width, MHz', fontsize=12)
            png_name = test.pics_path / 'delta_f'
    ax.set_xlabel('Time, ns', fontsize=12)
    ax.set_title('{}. Magnetron'.format(exp))
    ax.legend(loc='best')
    ax.grid(which='both', axis='both')
    fig.savefig(png_name)
peak_table(191126, plot_type='delta')