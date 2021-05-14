import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
import numpy as np
import matplotlib.gridspec as gridspec
import csv
import openpyxl as excel


def oscillograms():
    test = ProcessSignal('200925')
    csv_types = test.read_type_file()
    csv_signals = csv_types['signal_files']
    csv_signal_nums = csv_types['signal_nums']
    excel_dicts = test.read_excel(csv_signal_nums)['numbers']
    magnetron_nums = excel_dicts['magnetron']
    for signal in csv_signals:
        file = test.open_file(signal, reduced=True)
        u = file['voltage']
        t = file['time']
        dt = file['time_resolution']
        u_filt = test.bandpass_filter(t, u, 2.725e9, 2.755e9)
        num = signal[3:6]
        if int(num) < 66:
            absorbers = '№2, 3'
        elif 66 <= int(num) < 125 and num in magnetron_nums:
            absorbers = 0
        elif 125 <= int(num) <= 184 and num in magnetron_nums:
            absorbers = '№2, 3, A'
        else:
            absorbers = 0
        pl_d_nums = test.read_excel(csv_signal_nums)['dicts'].keys()
        if num in pl_d_nums:
            pl_density = test.read_excel(csv_signal_nums)['dicts'][num]['Ток плазмы, А']
            test.oscill_picture(num, t, u, u_filt, pl_density, absorbers, save=True)
#oscillograms()


def one_signal_oscillogramm(signal_num, voltage_num, exp_num, central_freq=2.71e9):
    test = ProcessSignal(str(exp_num))
    signal_file = test.open_file(signal_num, reduced=False)
    print('Signal {} data obtained'.format(signal_num[3:6]))
    pl_density = test.read_excel(signal_num)['dicts'][signal_num[3:6]]['Ток плазмы, А']
    voltage_file = test.open_file(voltage_num, reduced=False)
    print('Voltage {} data obtained'.format(voltage_num[3:6]))
    signal_u = signal_file['voltage']
    signal_t = signal_file['time'] + 50e-9
    voltage_u = voltage_file['voltage']
    voltage_t = voltage_file['time']

    ind_y_shift = signal_t < 20e-9
    y_shift = np.mean(signal_u[ind_y_shift])
    u_shift = signal_u - y_shift
    try:
        filt_freq_min, filt_freq_max = central_freq - 15e6, central_freq + 15e6
        u_filt = test.fft_filter(signal_t, signal_u, filt_freq_min, filt_freq_max)
        #u_filt_1 = test.bandpass_filter(signal_t, signal_u, 2.695e9, 2.725e9)
        print('Filtering done')
    except:
        pass
    #pl_density = test.read_excel(signal_num)['dicts'][signal_num[3:6]]['Ток плазмы, А']
    #magnetron_abs = test.read_excel(signal_num)['dicts'][signal_num[3:6]]['Поглотители в тракте магнетрона']
    fig = plt.figure(num=1, dpi=150)
    ax = fig.add_subplot(111)
    print('Creating a picture...')
    line1, = ax.plot(signal_t / 1e-9, u_shift, linewidth=0.7, color='dodgerblue')
    #line1, = ax.plot(signal_t / 1e-9, signal_u, linewidth=0.7, color='dodgerblue')
    line1.set_label('СВЧ сигнал')

    try:
        line2, = ax.plot(signal_t / 1e-9, u_filt, linewidth=0.7, color='red')
        line2.set_label('Фильтрованный сигнал (2,71 +/- 0.015) ГГц')
        
        #line4, = ax.plot(signal_t / 1e-9, u_filt_1, linewidth=0.7, color='green')
        #line4.set_label('Bandpass (2,74 +/- 0.03) ГГц')
    except:
        pass
    line3, = ax.plot(voltage_t / 1e-9, voltage_u, linewidth=0.7, color='blue')
    line3.set_label('Импульс напряжения')
    ax.set_xlabel(r'$Время, нс$', fontsize=14, fontweight='black')
    ax.set_ylabel(r'$Напряжение, В$', fontsize=14, fontweight='black')
    ax.set_ylim(bottom=-2.5, top=2.5)
    #ax.legend()
    ax.grid(which='both', axis='both')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    '''
    if magnetron_abs is not None:
        plt.title('№ {}, n = {}, Поглотители = {}'.format(signal_num[3:6], pl_density, magnetron_abs))
    else:
        plt.title('№ {}, n = {}'.format(signal_num[3:6], pl_density))
    '''
    plt.title(f'№ {signal_num[3:6]}-{voltage_num[3:6]},n = {pl_density}')
    png_name = test.signal_pics_path / f'all_ocs_{pl_density}_{signal_num[3:6]}.png'
    fig.savefig(png_name)
    plt.close(fig)


def exp_oscillogramms(exp_num, last_num=0, first_num=0):
    test = ProcessSignal(str(exp_num))
    signal_nums_0 = test.read_type_file()['signal_files']
    voltage_nums_0 = test.read_type_file()['voltage_files']
    if last_num == 0:
        signal_nums, voltage_nums = signal_nums_0, voltage_nums_0
    else:
        signal_nums = [signal_nums_0[i] for i in range(len(signal_nums_0)) if first_num < int(signal_nums_0[i][3:6]) < last_num]
        voltage_nums = [voltage_nums_0[i] for i in range(len(voltage_nums_0)) if first_num < int(voltage_nums_0[i][3:6]) < last_num]
    print('Signals:', signal_nums)
    for element in zip(signal_nums, voltage_nums):
        one_signal_oscillogramm(element[0], element[1], exp_num)


#exp_oscillogramms(210423, last_num=61, first_num=0)


def magnetron_osc(exp_num):
    test = ProcessSignal(str(exp_num))
    path = r'C:\Users\d_Nice\Documents\SignalProcessing\2021\{}'.format(str(exp_num))
    file_list = os.listdir(path)
    for file in file_list:
        if 'str' in file:
            data = test.open_file(file, reduced=False)
            t, u = data['time'], data['voltage']
            fig = plt.figure(num=1, dpi=150)
            ax = fig.add_subplot(111)
            print(f'Creating a picture {file}...')
            line1, = ax.plot(t / 1e-9, u, linewidth=0.7)
            ax.set_xlabel(r'$Время, нс$', fontsize=14, fontweight='black')
            ax.set_ylabel(r'$Напряжение, В$', fontsize=14, fontweight='black')
            ax.grid(which='both', axis='both')
            ax.set_ylim(bottom=-4, top=4)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.title(f'№ {file[3:6]}')
            png_name = test.signal_pics_path / 'ocs_{}'.format(file[3:6])
            fig.savefig(png_name)
            plt.close(fig)

#magnetron_osc(210423)


def file_nums_oscillogramms(exp_num, file_nums):
    print(file_nums)
    test = ProcessSignal(str(exp_num))
    for file_num in file_nums:
        #file_name = f'str{file_num:03d}.csv'
        file_name = f'str{file_num}.csv'
        data = test.open_file(file_name, reduced=False)
        t, u = data['time'], data['voltage']
        fig = plt.figure(num=1, dpi=150)
        ax = fig.add_subplot(111)
        print(f'Creating a picture {file_num}...')
        line1, = ax.plot(t / 1e-9, u, linewidth=0.7)
        ax.set_xlabel(r'$Время, нс$', fontsize=14, fontweight='black')
        ax.set_ylabel(r'$Напряжение, В$', fontsize=14, fontweight='black')
        ax.grid(which='both', axis='both')
        min_y, max_y = min(u) - 0.25, max(u) + 0.25
        ax.set_ylim(bottom=min_y, top=max_y)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(f'№ {file_num}')
        png_name = test.signal_pics_path / f'ocs_{file_num}.png'
        fig.savefig(png_name)
        plt.close(fig)


#file_nums_oscillogramms(210414, [f'{i:03d}' for i in range(65, 91)])


def pl_densities_from_excel(exp_num, file_nums_list):
    proc = ProcessSignal(f'{exp_num}')
    all_files_list = os.listdir(proc.exp_file_path)
    pl_d_vals = []
    for file_num in file_nums_list:
        file_name = f'str{file_num}.csv'
        try:
            pl_density = proc.read_excel(file_name)['dicts'][file_num]['Ток плазмы, А']
        except:
            pl_density = 0
        pl_d_vals.append(pl_density)
    return(pl_d_vals)

#pl_densities_from_excel(210414, [f'{i:03d}' for i in range(121, 151)])


def file_nums_oscillogramms_density_vals(exp_num, file_nums):
    print(file_nums)
    test = ProcessSignal(str(exp_num))
    density_vals = pl_densities_from_excel(exp_num, file_nums)
    density_mtrx = np.reshape(np.asarray(density_vals), (int(len(density_vals)/2), 2))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(density_vals)/2), 2))
    print(density_mtrx[0, ])
    central_freq = 2.714E9
    for j in range(nums_mtrx.shape[0]):
        gs = gridspec.GridSpec(2, 1)
        fig = plt.figure(num=1, dpi=200, figsize=[11.69, 12.27])
        for i in range(nums_mtrx.shape[1]):
            file_num = f'{nums_mtrx[j, i]}'
            file_name = f'str{file_num}.csv'
            data = test.open_file(file_name, reduced=False)
            if int(nums_mtrx[j, i]) % 2 == 0:
                t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
            else:
                t, u = data['time'], (data['voltage'] - np.mean(data['voltage']))
            filt_freq_min, filt_freq_max = central_freq - 15e6, central_freq + 15e6
            u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)

            ax = fig.add_subplot(gs[i, 0])
            print(f'Creating a picture {file_num}...')
            line1, = ax.plot(t / 1e-9, u, linewidth=0.7)
            line2, = ax.plot(t / 1e-9, u_filt, linewidth=0.7)
            ax.set_xlabel(r'$Время, нс$', fontsize=14, fontweight='black')
            ax.set_ylabel(r'$Напряжение, В$', fontsize=14, fontweight='black')
            ax.grid(which='both', axis='both')
            min_y, max_y = min(u) - 0.25, max(u) + 0.25
            ax.set_ylim(bottom=min_y, top=max_y)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            ax.set_title(f'№ {file_num}, n={density_mtrx[j, i]}')
        png_name = test.signal_pics_path / f'ocs_{density_mtrx[j,0]}_{nums_mtrx[j, 0]}_{nums_mtrx[j, 1]}.png'
        fig.savefig(png_name)
        plt.close(fig)


file_nums_oscillogramms_density_vals(210423, [f'{i:03d}' for i in range(131, 187)])


def write_2_antennas_osc_csv(exp_num, file_nums):
    print(f'Experiment {exp_num}')
    test = ProcessSignal(str(exp_num))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(file_nums) / 2), 2))
    central_freq = 2.714E9
    for j in range(nums_mtrx.shape[0]):
        antennas_data_dict = {}
        for i in range(nums_mtrx.shape[1]):
            file_num = f'{nums_mtrx[j, i]}'
            file_name = f'str{file_num}.csv'
            data = test.open_file(file_name, reduced=False)
            filt_freq_min, filt_freq_max = central_freq - 15e6, central_freq + 15e6
            if int(nums_mtrx[j, i]) % 2 == 0:
                t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
                u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                antennas_data_dict['t_main'], antennas_data_dict['u_main'] = t, u
                antennas_data_dict['u_filt_main'] = u_filt
                plt.plot(t, u)
            else:
                #t, u = data['time'], (data['voltage'] - np.mean(data['voltage'])) / 3.16
                t, u = data['time'], (data['voltage'] - np.mean(data['voltage']))
                u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                antennas_data_dict['t'], antennas_data_dict['u'] = t, u
                antennas_data_dict['u_filt'] = u_filt
                plt.plot(t + 4.347E-9, u)

        plt.show()
        '''
        print(f'Writing {nums_mtrx[j, i]} csv file...')
        print(antennas_data_dict.keys())

        dif_file_path = test.exp_file_path / '2_antennas_CSV'
        dif_file_path.mkdir(parents=True, exist_ok=True)
        file = open(str(dif_file_path / f'{nums_mtrx[j, 0]}_{nums_mtrx[j, 1]}.csv'), 'w', newline='')
        with file:
            writer = csv.writer(file)
            writer.writerow(['t_main(s)', 'V_main(V)', 'V_main_filt(V)', 't(s)', 'V(V)', 'V_filt(V)'])
            for i in range(1, max(len(antennas_data_dict['t_main']), len(antennas_data_dict['t']))):
                writer.writerow([antennas_data_dict['t_main'][i], antennas_data_dict['u_main'][i], antennas_data_dict['u_filt_main'][i],
                                 antennas_data_dict['t'][i], antennas_data_dict['u'][i], antennas_data_dict['u_filt'][i]])
        '''

#write_2_antennas_osc_csv(210423, [165, 166])


def find_zeros(time, voltage):
    time, voltage = np.asarray(time), np.asarray(voltage)
    voltage_signs = np.sign(voltage)
    l = 0
    for i in range(len(voltage_signs) - 1):
        if l < 1:
            if voltage_signs[i] != voltage_signs[i + 1] and voltage_signs[i + 1] > 0:
                start_ind = i
                l += 1
    zero_inds = [start_ind]
    m = 0
    for j in range(start_ind + 1, len(voltage_signs) - 1):
        if m < 2:
            if voltage_signs[j] != voltage_signs[j + 1] and voltage_signs[j + 1] != 0:
                zero_inds.append(j)
                m += 1
    time_part_max, volt_part_max = time[zero_inds[0]:zero_inds[1]], voltage[zero_inds[0]:zero_inds[1]]
    time_part_min, volt_part_min = time[zero_inds[1]:zero_inds[2]], voltage[zero_inds[1]:zero_inds[2]]
    #print(zero_inds)
    ind_max, ind_min = np.argmax(volt_part_max), np.argmin(volt_part_min)
    time_max, volt_max = time_part_max[ind_max], volt_part_max[ind_max]
    time_min, volt_min = time_part_min[ind_min], volt_part_min[ind_min]
    #plt.plot(time, voltage)
    #plt.plot(time_part_min, volt_part_min)
    #plt.plot(time_part_max, volt_part_max)
    #plt.show()
    max_min_dict = {'time_max': time_max,
                    'time_min': time_min,
                    'volt_max': volt_max,
                    'volt_min': volt_min}
    return max_min_dict


def two_antennas_max_table(exp_num, file_nums):
    print(f'Experiment {exp_num}')
    test = ProcessSignal(str(exp_num))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(file_nums) / 2), 2))
    print(nums_mtrx)
    central_freq = 2.714E9
    compare_pts = [75 * 1e-9, 155 * 1e-9, 200 * 1e-9]

    ex_table = excel.Workbook()
    ex_table.create_sheet(title='Integral', index=0)
    sheet = ex_table['Integral']
    sheet['A1'] = exp_num
    sheet['A2'] = 'No (центр/бок)'
    sheet['B2'] = 'n, отн.ед.'
    sheet['C1'] = 'Без фильтра'
    #sheet['С2'] = f't_1 = {np.round(compare_pts[0] / 1e-9, 0)}'
    sheet['D2'] = f't_2 = {np.round(compare_pts[1] / 1e-9, 0)} c'
    sheet['E2'] = f't_3 = {np.round(compare_pts[2] / 1e-9, 0)} c'

    sheet['F1'] = 'Фильтрованный'
    sheet['F2'] = f't_1 = {np.round(compare_pts[0] / 1e-9, 0)} c'
    sheet['G2'] = f't_2 = {np.round(compare_pts[1] / 1e-9, 0)} c'
    sheet['H2'] = f't_3 = {np.round(compare_pts[2] / 1e-9, 0)} c'

    for pt in compare_pts:
        for j in range(nums_mtrx.shape[0]):
            antennas_data_dict = {}
            for i in range(nums_mtrx.shape[1]):
                print(nums_mtrx[j, i], 'pt=', pt)
                file_num = f'{nums_mtrx[j, i]}'
                file_name = f'str{file_num}.csv'
                data = test.open_file(file_name, reduced=False)
                filt_freq_min, filt_freq_max = central_freq - 15e6, central_freq + 15e6
                if int(nums_mtrx[j, i]) % 2 == 0:
                    t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
                    u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                    pl_density = test.read_excel(file_name)['dicts'][file_num]['Ток плазмы, А']
                    ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
                    t, u = t[ind_mask], u[ind_mask]
                    peak_dict = find_zeros(t, u)
                    time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict['time_min'], peak_dict['volt_min']
                    delta_main = volt_max - volt_min

                    u_filt = u_filt[ind_mask]
                    filt_dict = find_zeros(t, u_filt)
                    time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
                    time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
                    delta_main_filt = volt_max_filt - volt_min_filt

                    cell = sheet.cell(row=j + 3, column=2)
                    cell.value = f'{pl_density}'

                    plt.plot(t, u)
                    plt.plot(time_min, volt_min,  marker='o')
                    plt.plot(time_max, volt_max,  marker='o')
                else:
                    t, u = data['time'] + 4.347E-9, (data['voltage'] - np.mean(data['voltage']))
                    u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                    ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
                    t, u = t[ind_mask], u[ind_mask]
                    peak_dict = find_zeros(t, u)
                    time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict['time_min'], peak_dict['volt_min']
                    delta_sub = volt_max - volt_min

                    u_filt = u_filt[ind_mask]
                    filt_dict = find_zeros(t, u_filt)
                    time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
                    time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
                    delta_sub_filt = volt_max_filt - volt_min_filt

                    plt.plot(t, u_filt)
                    plt.plot(time_min_filt, volt_min_filt, marker='o')
                    plt.plot(time_max_filt, volt_max_filt, marker='o')
            u_relat_filt = delta_main_filt / delta_sub_filt
            u_relat = delta_main / delta_sub

            if pt == compare_pts[0]:
                column = 3
                column_filt = 6
            elif pt == compare_pts[1]:
                column = 4
                column_filt = 7
            else:
                column = 5
                column_filt = 8
            cell = sheet.cell(row=j + 3, column=column)
            cell.value = f'{np.round(u_relat, 3)}'
            cell = sheet.cell(row=j + 3, column=column_filt)
            cell.value = f'{np.round(u_relat_filt, 3)}'

            print(u_relat_filt)
            plt.title(nums_mtrx[j, i])
            #plt.show()
            cell = sheet.cell(row=j + 3, column=1)
            cell.value = f'{int(nums_mtrx[j, i]) -1 } / {nums_mtrx[j, i]}'

    path = test.excel_folder_path / f'Ampl_{exp_num}_2_antennas.xlsx'
    ex_table.save(path)
#two_antennas_max_table(210423, [f'{i:03d}' for i in range(131, 187)])


def two_antennas_ampl_relate(exp_num, file_nums, pts):
    print(f'Experiment {exp_num}')
    test = ProcessSignal(str(exp_num))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(file_nums) / 2), 2))
    print(nums_mtrx)
    central_freq = 2.714E9
    compare_pts = [pts[i] * 1E-9for i in range(len(pts))]
    for pt in compare_pts:
        for j in range(nums_mtrx.shape[0]):
            for i in range(nums_mtrx.shape[1]):
                print(nums_mtrx[j, i], 'pt=', pt)
                file_num = f'{nums_mtrx[j, i]}'
                file_name = f'str{file_num}.csv'
                data = test.open_file(file_name, reduced=False)
                filt_freq_min, filt_freq_max = central_freq - 15e6, central_freq + 15e6
                if int(nums_mtrx[j, i]) % 2 == 0:
                    t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
                    u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                    pl_density = test.read_excel(file_name)['dicts'][file_num]['Ток плазмы, А']
                    ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
                    t, u = t[ind_mask], u[ind_mask]
                    peak_dict = find_zeros(t, u)
                    time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict['time_min'], peak_dict['volt_min']
                    delta_main = volt_max - volt_min

                    u_filt = u_filt[ind_mask]
                    filt_dict = find_zeros(t, u_filt)
                    time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
                    time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
                    delta_main_filt = volt_max_filt - volt_min_filt

                    plt.plot(t, u)
                    plt.plot(time_min, volt_min,  marker='o')
                    plt.plot(time_max, volt_max,  marker='o')
                else:
                    t, u = data['time'] + 4.347E-9, (data['voltage'] - np.mean(data['voltage']))
                    u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                    ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
                    t, u = t[ind_mask], u[ind_mask]
                    peak_dict = find_zeros(t, u)
                    time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict['time_min'], peak_dict['volt_min']
                    delta_sub = volt_max - volt_min

                    u_filt = u_filt[ind_mask]
                    filt_dict = find_zeros(t, u_filt)
                    time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
                    time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
                    delta_sub_filt = volt_max_filt - volt_min_filt

                    plt.plot(t, u_filt)
                    plt.plot(time_min_filt, volt_min_filt, marker='o')
                    plt.plot(time_max_filt, volt_max_filt, marker='o')
            u_relat_filt = delta_main_filt / delta_sub_filt
            u_relat = delta_main / delta_sub
            print('Отн_без_фильтра:', u_relat)
            print('Oтн_фильтр:', u_relat_filt)

#two_antennas_ampl_relate(210423, [139, 140], [200, 300, 400, 500, 550])


def signal_envelope(t, u, dt, time_frame=1e-8):
    num_of_pts = int(time_frame / dt)
    times = t[::num_of_pts]
    time_step = np.abs(times[0] - times[1])
    t_env = []
    u_env = []
    for time in times:
        time_lim = time + time_step
        ind_t = np.logical_and(t > time, t <= time_lim)
        t_frame = t[ind_t]
        u_frame = u[ind_t]
        ind_max = np.argmax(u_frame)
        t_env.append(t_frame[ind_max])
        u_env.append(u_frame[ind_max])
    envelope_dict = {'env_time': t_env,
                     'env_voltage': u_env}
    return envelope_dict


def envelope_max_part(self, t, u, dt, max_part):
    env = self.signal_envelope(t, u, dt)
    env_t = np.array(env['env_time'])
    env_u = np.array(env['env_voltage'])
    max_env_u = np.max(env_u)
    ind_lim = env_u > max_part * max_env_u
    t_lim = env_t[ind_lim]
    ind_stop = len(t_lim - 1)
    l = 0
    for i in range(len(t_lim) - 1):
        if l < 1:
            if t_lim[i + 1] - t_lim[i] > 90e-9:
                ind_stop = i
                l += 1
    try:
        t_bond = t_lim[:ind_stop:]
        u_zeros = np.zeros(t_bond.size)
        t_min, t_max = t_bond[0], t_bond[-1]
        ind_use = np.logical_and(t >= t_min, t <= t_max)
        t_use = t[ind_use]
        u_use = u[ind_use]
        signal_part_dict = {'signal_time': t_use,
                            'signal_voltage': u_use}
        return signal_part_dict
    except:
        pass

def file_nums_oscillogramms_10_ns(exp_num, file_nums, start_time, density_vals=False, filt=False):
    print(file_nums)
    start_time = start_time * 1e-9
    end_time = start_time + 10e-9
    test = ProcessSignal(str(exp_num))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(file_nums) / 2), 2))
    if density_vals:
        density_vals = pl_densities_from_excel(exp_num, file_nums)
        density_mtrx = np.reshape(np.asarray(density_vals), (int(len(density_vals) / 2), 2))
    for j in range(nums_mtrx.shape[0]):
        fig = plt.figure(num=1, dpi=200)
        for i in range(nums_mtrx.shape[1]):
            file_num = f'{nums_mtrx[j, i]}'
            file_name = f'str{file_num}.csv'
            data = test.open_file(file_name, reduced=False)
            if int(file_num) % 2 == 0:
                t, u = data['time'], (data['voltage'] - np.mean(data['voltage'])) / 3
                '''
                filt_freq_min, filt_freq_max = 2.714e9 - 15e6, 2.714e9 + 15e6
                u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                '''
                time_inds = np.logical_and(t >= start_time, t <= end_time)
            else:
                time_shift = 4.35 * 1e-9
                t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
                '''
                filt_freq_min, filt_freq_max = 2.714e9 - 15e6, 2.714e9 + 15e6
                u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
                '''
                time_inds = np.logical_and(t >= start_time - time_shift, t <= end_time - time_shift)
            time_part_centr, u_part_centr = t[time_inds], u[time_inds]
            ax = fig.add_subplot(111)
            print(f'Creating a picture {file_num}...')
            if int(file_num) % 2 == 0:
                line1, = ax.plot(time_part_centr / 1e-9, u_part_centr, linewidth=1.2, color='red')
                line1.set_label('Сигнал на оси')
                #min_y, max_y = min(u_part_centr) - 0.25, max(u_part_centr) + 0.25
            else:
                line1, = ax.plot((time_part_centr + time_shift) / 1e-9, u_part_centr, linewidth=1.2)
                line1.set_label('Сигнал сбоку')
                min_y, max_y = min(u_part_centr) - 0.25, max(u_part_centr) + 0.25
        ax.set_xlabel(r'$Время, нс$', fontsize=14, fontweight='black')
        ax.set_ylabel(r'$Напряжение, В$', fontsize=14, fontweight='black')
        ax.grid(which='both', axis='both')
        ax.set_ylim(bottom=min_y, top=max_y)
        ax.legend()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        if density_vals:
            ax.set_title(f'№ {nums_mtrx[j, 0]} - {nums_mtrx[j, 1]} (магнетрон), n={density_mtrx[j, i]}, {int(start_time / 1e-9)} - {(start_time /1e-9) + 10}ns')
        else:
            ax.set_title(f'№ {nums_mtrx[j, 0]} - {nums_mtrx[j, 1]}, {int(start_time / 1e-9)} - {(start_time / 1e-9) + 10} ns')
        png_name = test.signal_pics_path / f'{exp_num}_{nums_mtrx[j, 0]}_{nums_mtrx[j, 1]}_магнетрон_{int(start_time / 1e-9)}ns.png'
        fig.savefig(png_name)
        plt.close(fig)


#file_nums_oscillogramms_10_ns(210423, [f'{i:03d}' for i in range(165, 167)], start_time=120, density_vals=True)


def rename_osc_files(exp_num):
    proc = ProcessSignal(f'{exp_num}')
    all_files_list = os.listdir(proc.signal_pics_path)
    osc_files_list = [all_files_list[i] for i in range(len(all_files_list)) if 'png' in all_files_list[i] and 'all' not in all_files_list[i]]
    print('files:', osc_files_list)
    density_vals = pl_densities_from_excel(exp_num)
    for i, file in enumerate(osc_files_list):
        new_file_name = file.replace(f'{file}', f'{file[4:7]}_{density_vals[i]}.png')
        old_file_path = proc.signal_pics_path / file
        new_file_path = proc.signal_pics_path / new_file_name
        os.rename(old_file_path, new_file_path)


#rename_osc_files(210414)


