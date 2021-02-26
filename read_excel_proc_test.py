import openpyxl as xl
from openpyxl.utils import get_column_letter
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

main_excel_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\210220\210220.xlsx')
type_file_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\210220\Excel\types.xlsx')

def read_excel():
    wb = xl.load_workbook(main_excel_path)
    sheet = wb['Лист1']
    max_row = sheet.max_row
    max_col = sheet.max_column
    min_row = 1

    for i in range(1, max_col):
        for j in range(1, max_row):
            cell = sheet.cell(row=j, column=i)
            cell_val = str(cell.value)
            if 'файл' in cell_val and min_row == 1:
                min_row, min_col = j, i
    min_col_letter, max_col_letter = get_column_letter(min_col), get_column_letter(max_col)
    #cells_data = [sheet[f'{min_col_letter}{min_row}':f'{max_col_letter}{max_row}'][0][i][j].value for i, j in range(max_col, max_row)]
    data_keys = [sheet[f'{min_col_letter}{min_row}':f'{max_col_letter}{min_row}'][0][i].value for i in range(max_col)]
    rows_data = []
    for row in range(min_row, max_row):
        row_data = [sheet[f'{min_col_letter}{row}':f'{max_col_letter}{row}'][0][i].value for i in range(max_col)]
        print(row_data)
    print(cells_data)
    keys = ['file_nums', 'plasma_source_heat', 'plasma_current', 'pressure', 'magnetic_field', 'magnetron_delay',
            'magnetron_start', 'GIN_voltage', 'lamp_detector', 'magnetron_charge_voltage', 'magnetron_absorbers',
            'magnetron_voltage', 'U_in', 'comment']
    for k in range(1, max_col + 1):
        cell = sheet.cell(row=min_row, column=k)
        cell_name = cell.value
        if 'файл' in cell_name:
            cell_key = 'num'
        keys.append(cell_key)
read_excel()
'''
    row_dicts = []
    for l in range(min_row + 1, max_row + 1):
        row_dict = {}
        vals = []
        for m in range(1, max_col + 1):
            cell = sheet.cell(row=l, column=m)
            cell_val = cell.value
            vals.append(cell_val)
        for i, key in enumerate(keys):
            row_dict[key] = vals[i]
        row_dicts.append(row_dict)

    plasma_dicts = []
    reb_dicts = []
    magnetron_dicts = []
    for row_dict in row_dicts:
        fnum = row_dict['Номер файла']
        gin_voltage = row_dict['Напряжение ГИНа']
        d_plasma = row_dict['Ток плазмы, А']
        heating = row_dict['Накал']
        magnetron_delay = row_dict['Задержка магнетрона, нс']
        # magnetron_in_voltage = row_dict['Входное напряжение магнетрона, В']
        comment = 0
        if isinstance(fnum, int):
            if isinstance(d_plasma, float) or isinstance(d_plasma, int):
                if comment != 'except':
                    if isinstance(magnetron_delay, float) or isinstance(magnetron_delay, int):
                        magnetron_dicts.append(row_dict)
                    else:
                        plasma_dicts.append(row_dict)
            elif heating is None and isinstance(gin_voltage, float):
                reb_dicts.append(row_dict)

    useful_dicts = [plasma_dicts, reb_dicts, magnetron_dicts]

    proc_signal_dicts = {}
    use_nums = {}
    list_signals = ['noise', 'reb', 'magnetron']
    for i, use_dict in enumerate(useful_dicts):
        # use_csv_files = []
        fnums = []
        for row_dict in use_dict:
            fnum = row_dict['Номер файла']
            if fnum < 100:
                num_1 = fnum // 10
                num_2 = fnum % 10
            elif fnum > 10000:
                num_1 = fnum // 1000
                num_2 = fnum % 1000
            else:
                num_1 = fnum // 100
                num_2 = fnum % 100
            fnum_1 = '{:03d}'.format(num_1)
            fnum_2 = '{:03d}'.format(num_2)

            if fnum_1 in csv_signal_nums:
                row_dict['Номер файла'] = [fnum_1]
                fnums.append(fnum_1)
                # use_csv_files.append(row_dict)
                proc_signal_dicts[fnum_1] = row_dict
            elif fnum_2 in csv_signal_nums:
                row_dict['Номер файла'] = [fnum_2]
                fnums.append(fnum_2)
                # use_csv_files.append(row_dict)
                proc_signal_dicts[fnum_2] = row_dict
        # proc_signal_dicts[list_signals[i]] = use_csv_files
        use_nums[list_signals[i]] = fnums

    excel_results = {'dicts': proc_signal_dicts,
                     'numbers': use_nums}
    return excel_results
'''


def part_fft(self, csv_signals, interest_nums,
             part_nums, fft_type='part',
             block_full=False, block_part=True,
             peak=False, noise=False):
    ex_table = xl.Workbook()
    ex_table.create_sheet(title='Full FFT', index=0)
    ex_table.create_sheet(title='Part FFT', index=1)
    sheet1 = ex_table['Full FFT']
    sheet2 = ex_table['Part FFT']
    sheet2['A1'] = 'File Number'
    sheet1['A1'] = 'File Number'
    signal_dict = {}
    row_f = 1
    row_p = 1
    for j, csv_signal in enumerate(csv_signals):
        signal_num = csv_signal[3:6]
        if signal_num in interest_nums or signal_num in part_nums:
            use_signal = self.open_file(csv_signal, reduced=True)
            use_t = use_signal['time']
            use_u = use_signal['voltage']
            dt = use_signal['time_resolution']
            fig = plt.figure(num=1, dpi=300)
            ax = fig.add_subplot(111)
            line = None
            ax.set_prop_cycle(color=['mediumseagreen', 'dodgerblue', 'indigo'])
            if fft_type == 'full' or fft_type == 'both':
                if signal_num in interest_nums:
                    print('{} am in interest_nums'.format(signal_num))
                    if block_full is True:
                        filtered_u = self.fft_filter(use_t, use_u, 2.69e9, 2.74e9, filt_type='bandstop')
                        fft_results = self.fft_amplitude(use_t, filtered_u)
                    else:
                        fft_results = self.fft_amplitude(use_t, use_u)
                    pl_density = self.read_excel(interest_nums)['dicts'][signal_num]['Ток плазмы, А']
                    freq = fft_results['frequency']
                    amp = fft_results['amplitude']

                    if peak:
                        peak_inds = np.logical_and(2.69e9 < freq, freq < 2.74e9)
                        peak_freqs = freq[peak_inds]
                        peak_amps = amp[peak_inds]
                        mean_freq = self.mean_frequency(peak_freqs, peak_amps)
                        # line, = ax.plot(peak_freqs, peak_amps, linewidth=1.2)
                        line, = ax.plot(freq, amp, linewidth=1.2)
                        peak_max = np.round(np.max(peak_amps), 2)
                    else:
                        mean_freq = self.mean_frequency(freq, amp)
                        line, = ax.plot(freq, amp, linewidth=0.7)
                        line, = ax.plot(freq, amp, linewidth=0.7)
                    try:
                        spectrum_mean_freq = mean_freq['mean_freq']
                        line.set_label(r'$f = {} GHz$'.format(spectrum_mean_freq))

                        value_f = str('f, GHz')
                        cell_f = sheet1.cell(row=1, column=2)
                        cell_f.value = value_f

                        cell_pl = sheet1.cell(row=1, column=3)
                        cell_pl.value = str('n, arb.units')

                        row_f = row_f + 1
                        cell_name = 'A{}'.format(row_f)
                        sheet1[str(cell_name)] = '{}'.format(signal_num)
                        value = spectrum_mean_freq
                        cell = sheet1.cell(row=row_f, column=2)
                        cell.value = value

                        value_pl = pl_density
                        cell = sheet1.cell(row=row_f, column=3)
                        cell.value = value_pl
                    except TypeError:
                        pass

            if fft_type == 'part' or fft_type == 'both':
                if signal_num in part_nums:
                    print('{} in part_nums'.format(signal_num))
                    if block_part is True:
                        filtered_u = self.bandstop_filter(use_t, use_u, 2.709e9, 2.769e9)
                        part_signal = self.signal_parts(use_t, filtered_u, dt)
                    else:
                        part_signal = self.signal_parts(use_t, use_u, dt)
                    part_keys = part_signal.keys()
                    part_freq_dict = {}

                    row_p = row_p + 1
                    cell_name = 'A{}'.format(row_p)
                    sheet2[str(cell_name)] = '{}'.format(signal_num)
                    for i, part_key in enumerate(part_keys):
                        k = i + 1
                        col = k + 1
                        part_time = part_signal[part_key]['time']
                        part_voltage = part_signal[part_key]['voltage']
                        fft_results = self.fft_signal(part_time, part_voltage, dt)
                        freq = fft_results['frequency']
                        amp = fft_results['amplitude']
                        line, = ax.plot(freq, amp, linewidth=0.7)

                        mean_freq = self.mean_frequency(freq, amp)
                        part_freq = mean_freq['mean_freq']

                        value_f = str('f{}'.format(k))
                        cell_f = sheet2.cell(row=1, column=col)
                        cell_f.value = value_f

                        value = part_freq
                        cell = sheet2.cell(row=row_p, column=col)
                        cell.value = value
                        part_freq_dict['f{}, GHz'.format(k)] = part_freq
                        line.set_label(r'$f_{}= {} GHz$'.format(k, part_freq))
                    pl_density = self.read_excel(part_nums)['dicts'][signal_num]['Ток плазмы, А']
                    signal_dict[signal_num] = part_freq_dict
            if fft_type == 'part' or fft_type == 'both':
                if signal_num in part_nums:
                    ax.set_title(r'$File\/Number = {} (Part FFT, n={})$'.format(signal_num, pl_density))
            if fft_type == 'full' or fft_type == 'both':
                if signal_num in interest_nums:
                    if noise:
                        ax.set_title(r'$File\/Number = {},\/ (n={}) $'.format(signal_num, pl_density))
                    else:
                        # absorbers = self.read_excel(part_nums)['dicts'][signal_num]['Поглотители в тракте магнетрона']
                        ax.set_title(r'$File\/Number = {},\/ \/ n={} $'.format(signal_num, pl_density))

            if line is not None:
                ax.set_ylim(bottom=0)
                if peak:
                    ax.set_xlim(left=2.68e9, right=2.74e9)
                else:
                    ax.set_xlim(left=0, right=4e9)
                ax.grid(which='both', axis='both')
                ax.set_xlabel(r'$Frequency, GHz$')
                ax.set_ylabel(r'$Amplitude$')
                ax.legend()
                if fft_type == 'part' or fft_type == 'both':
                    if signal_num in part_nums:
                        png_name = self.fft_part_path / signal_num
                if fft_type == 'full' or fft_type == 'both':
                    if signal_num in interest_nums:
                        if noise:
                            png_name = self.fft_reb / 'reb_{}'.format(signal_num)
                            table_path = self.fft_reb / 'fft_reb_{}.xlsx'.format(signal_num)
                        else:
                            if block_full:
                                png_name = self.fft_noise_base / 'noise_base_{}'.format(signal_num)
                                table_path = self.fft_noise_base / 'fft_noise_base_{}.xlsx'.format(signal_num)
                            else:
                                png_name = self.fft_magnetron / 'magnetron_{}'.format(signal_num)
                                table_path = self.fft_magnetron / 'fft_magnetron_{}.xlsx'.format(signal_num)
                        if peak:
                            png_name = self.fft_peak / 'peak_{}_1'.format(signal_num)
                            table_path = self.fft_peak / 'peak_{}.xlsx'.format(signal_num)
                fig.savefig(png_name)
                plt.close(fig)
        print('Circle {} complete'.format(j))
    ex_table.save(table_path)


def fft_full(self, magnetron_full=False, magnetron_noise_base=True, peak=False, reb_noise=False, peak_freq=2.71e9,
             peak_gate=50e6):
    types = self.read_type_file()
    csv_signals, csv_signal_nums = types['signal_files'], types['signal_nums']
    excel_results = self.read_excel(csv_signal_nums)['numbers']
    noise_nums = excel_results['noise']
    magnetron_nums = excel_results['magnetron']
    if reb_noise:
        nums = noise_nums
    else:
        nums = magnetron_nums
    for num in nums:
        file_name = f'str{num}.csv'
        file_data = self.open_file(file_name, reduced=True)
        t, u, dt = file_data['time'], file_data['voltage'], file_data['time_resolution']
        pl_density = self.read_excel(file_name)['dicts'][num]['Ток плазмы, А']
        fft_results = self.fft_amplitude(t, u)
        freqs, amps = fft_results['frequency'], fft_results['amplitude']
        mean_freq = self.mean_frequency(freqs, amps)
        spectrum_mean_freq = mean_freq['mean_freq']

        fig = plt.figure(num=1, dpi=300)
        ax = fig.add_subplot(111)

        if magnetron_noise_base:
            base_inds = np.logical_or(freqs < peak_freq - peak_gate, freqs > peak_freq + peak_gate)
            base_freqs, base_amps = freqs[base_inds], amps[base_inds]
            mean_freq = self.mean_frequency(base_freqs, base_amps)
            spectrum_mean_freq = mean_freq['mean_freq']
            line, = ax.plot(base_freqs, base_amps, linewidth=0.7, color='mediumseagreen')
        else:
            line, = ax.plot(freqs, amps, linewidth=0.7, color='mediumseagreen')
        if peak:
            ax.set_xlim(left=peak_freq-30e6, right=peak_freq+30e6)
        else:
            ax.set_xlim(left=0, right=4e9)
        ax.grid(which='both', axis='both')
        ax.set_xlabel(r'$Frequency, GHz$')
        ax.set_ylabel(r'$Amplitude$')
        ax.set_title(r'$№={}, n={}$'.format(num, pl_density))
        if not peak:
            line.set_label(r'$f = {} GHz$'.format(spectrum_mean_freq))
            ax.legend()
        if peak:
            png_name = self.fft_magnetron / 'peak_{num}'
        elif reb_noise:
            png_name = self.fft_magnetron / 'reb_noise_{num}'
        elif magnetron_noise_base:
            png_name = self.fft_magnetron / 'noise_base_{num}'
        else:
            png_name = self.fft_magnetron / 'amplifier_{num}'
        fig.savefig(png_name)
