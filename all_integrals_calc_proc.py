import numpy as np
import openpyxl as excel
import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal


def all_integrals_proc(exp_num, central_freq=2.714, band_half_width=0.05):
    proc = ProcessSignal(str(exp_num))
    csv_types = proc.read_type_file()
    csv_signals = csv_types['signal_files']
    csv_signal_nums = csv_types['signal_nums']
    excel_dicts = proc.read_excel(csv_signal_nums)['numbers']
    noise_nums = excel_dicts['noise']
    magnetron_nums = excel_dicts['magnetron']

    print('Noise nums:', noise_nums)
    print('Magsnetron nums:', magnetron_nums)
    nums, n_s, full_ints, peak_ints = [], [], [], []
    noise_ints, noise_ns, n_nums = [], [], []
    for signal in csv_signals:
        num = signal[3:6]
        if num in magnetron_nums:
            file = proc.open_file(signal, reduced=True, red_type=130)
            u = file['voltage']
            t = file['time']
            f_low = (central_freq - band_half_width) * 1e9
            f_high = (central_freq + band_half_width) * 1e9
            u_filt = proc.fft_filter(t, u, f_low, f_high)
            dt_2 = (t[-1] - t[0]) ** 2

            pl_density = proc.read_excel(signal)['dicts'][num]['Ток плазмы, А']

            integral = np.round(proc.e_square(t, u) / 1e-8, 3)
            filt_int = np.round(proc.e_square(t, u_filt) / 1e-8, 3)
            freqs, amps = proc.fft_amplitude(t, u)['frequency'], proc.fft_amplitude(t, u)['amplitude']
            peak_inds = np.logical_and(freqs >= f_low, freqs <= f_high)
            p_freqs, p_amps = freqs[peak_inds], amps[peak_inds]

            full_int = np.round(dt_2 * proc.e_square(freqs, amps) / 2e-8, 3)
            peak_int = np.round(dt_2 * proc.e_square(p_freqs, p_amps) / 2e-8, 3)

            n_s.append(pl_density)
            nums.append(num)
            full_ints.append(full_int)
            peak_ints.append(peak_int)

            print('peak_int = ', np.round(integral, 2), 'noise_int =', np.round(integral - filt_int, 2),
                  'noise_fft =', np.round(full_int - peak_int, 2))
        if num in noise_nums:
            file = proc.open_file(signal, reduced=True, red_type=130)
            u = file['voltage']
            t = file['time']
            dt = file['time_resolution']
            dt_2 = (t[-1] - t[0]) ** 2

            pl_density = proc.read_excel(signal)['dicts'][num]['Ток плазмы, А']

            integral = np.round(proc.e_square(t, u) / 1e-8, 3)
            freqs, amps = proc.fft_amplitude(t, u)['frequency'], proc.fft_amplitude(t, u)['amplitude']
            noise_int = np.round(dt_2 * proc.e_square(freqs, amps) / 2e-8, 3)
            if num in noise_nums:
                noise_ints.append(noise_int)
                noise_ns.append(pl_density)
                n_nums.append(num)

    ind_sort = np.argsort(np.asarray(n_s))
    ns = np.asarray(n_s)[ind_sort]
    nums = np.asarray(nums)[ind_sort]
    full_ints, peak_ints = np.asarray(full_ints)[ind_sort], np.asarray(peak_ints)[ind_sort]

    n_ind_sort = np.argsort(np.asarray(noise_ns))
    noise_ns, n_nums, noise_ints = np.asarray(noise_ns)[n_ind_sort], np.asarray(n_nums)[n_ind_sort], np.asarray(noise_ints)[n_ind_sort]
    ex_table = excel.Workbook()
    ex_table.create_sheet(title='Integral', index=0)
    sheet = ex_table['Integral']
    sheet['A1'] = 'Номер(шум)'
    sheet['B1'] = 'Плотность плазмы, отн.ед.'
    sheet['C1'] = 'W_2, *10-8'

    sheet['E1'] = 'Номер'
    sheet['F1'] = 'Плотность плазмы, отн.ед.'
    sheet['G1'] = 'Полный интеграл, *10-8'
    sheet['H1'] = 'W_f0, *10-8'
    sheet['I1'] = 'W_1, *10-8'

    for z in range(noise_ints.size):
        cell = sheet.cell(row=z + 2, column=1)
        cell.value = int(n_nums[z])
        cell = sheet.cell(row=z + 2, column=2)
        cell.value = noise_ns[z]
        cell = sheet.cell(row=z + 2, column=3)
        cell.value = noise_ints[z]

    for k in range(full_ints.size):
        cell = sheet.cell(row=k + 2, column=5)
        cell.value = int(nums[k])
        cell = sheet.cell(row=k + 2, column=6)
        cell.value = ns[k]
        cell = sheet.cell(row=k + 2, column=7)
        cell.value = full_ints[k]
        cell = sheet.cell(row=k + 2, column=8)
        cell.value = peak_ints[k]
        cell = sheet.cell(row=k + 2, column=9)
        cell.value = full_ints[k] - peak_ints[k]

    path = proc.excel_folder_path / f'Integrals_{exp_num}_130.xlsx'
    ex_table.save(path)


all_integrals_proc(210211)