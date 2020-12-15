from ProcessClass import ProcessSignal
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def max_freq_plasma(exp):
    test = ProcessSignal(str(exp))
    csv_types = test.read_type_file()
    csv_signals = csv_types['signal_files']
    csv_signal_nums = csv_types['signal_nums']
    print(csv_signals)
    excel_dicts = test.read_excel(csv_signal_nums)['numbers']
    magnetron_nums = excel_dicts['magnetron']
    pl_densities = np.zeros(len(magnetron_nums))
    magnetron_dicts = {}
    for i, num in enumerate(magnetron_nums):
        for csv_signal in csv_signals:
            csv_num = csv_signal[3:6]
            if csv_num == num and int(num) >= 126:
                magnetron_nums[i] = int(num)
                print('I am working on {} signal'.format(csv_num))
                pl_density = test.read_excel(csv_signal_nums)['dicts'][num]['Ток плазмы, А']
                pl_densities[i] = pl_density
                file = test.open_file(csv_signal, reduced=True)
                time = file['time']
                voltage = file['voltage']
                dt = file['time_resolution']

                fft_dict = test.fft_amplitude(time, voltage, dt)
                ampls = fft_dict['amplitude']
                freqs = fft_dict['frequency']
                # peak
                peak = max(ampls)
                max_ind = np.argmax(ampls)
                freq_peak = freqs[max_ind]
                # half_peak_lvl
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

                magnetron_dicts[csv_num] = {'peak': peak, 'freq_peak': np.round((freq_peak / 1e9), 3),
                                            'delta_f': np.round((delta_f / 1e6), 3), 'pl_d': pl_density}
    return(magnetron_dicts)
#max_freq_plasma(191120)

def peak_graph(exp, plot_type='peak'):
    test = ProcessSignal(str(exp))
    peak_dicts = max_freq_plasma(exp)
    keys = peak_dicts.keys()
    peaks = []
    freq_peaks = []
    delta_fs = []
    signal_nums = []
    pl_ds = []
    fig = plt.figure(num=1, dpi=300)
    ax = fig.add_subplot(111)
    for key in keys:
        peak_dict = peak_dicts[key]
        signal_nums.append(key)
        peaks.append(peak_dict['peak'])
        freq_peaks.append(peak_dict['freq_peak'])
        delta_fs.append(peak_dict['delta_f'])
        pl_ds.append(peak_dict['pl_d'])
    peaks = np.array(peaks)
    freq_peaks = np.array(freq_peaks)
    delta_fs = np.array(delta_fs)
    signal_nums = np.array(signal_nums)
    pl_ds = np.array(pl_ds)
    inds = np.argsort(pl_ds)

    peaks = peaks[inds]
    freq_peaks = freq_peaks[inds]
    delta_fs = delta_fs[inds]
    signal_nums = signal_nums[inds]
    pl_ds = pl_ds[inds]

    if plot_type is 'peak':
        ax.plot(pl_ds, peaks, marker='^', linewidth=0.9, ms=3)
        ax.set_ylabel('Амплитуда пика', fontsize=12)
        png_name = test.pics_path / 'Peak'
    if plot_type is 'peak_freq':
        ax.plot(pl_ds, freq_peaks, marker='^', linewidth=0.9, ms=3)
        ax.set_ylabel('Частота пика, ГГц', fontsize=12)
        png_name = test.pics_path / 'Peak_frequency'
    if plot_type is 'delta':
        ax.plot(pl_ds, delta_fs, marker='^', linewidth=0.9, ms=3)
        ax.set_ylabel('Ширина пика, МГц', fontsize=12)
        png_name = test.pics_path / 'delta_f'

    ax.set_xlabel('Плотность плазмы, отн.ед.', fontsize=12)
    #ax.set_xlim(right=11)
    ax.set_title('{}, поглотители №2, 3, А'.format(exp))
    ax.grid(which='both', axis='both')
    fig.savefig(png_name)
peak_graph(191106,  plot_type='delta')