import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
import numpy as np


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
    png_name = test.signal_pics_path / 'all_ocs_{}.png'.format(pl_density)
    fig.savefig(png_name)
    plt.close(fig)


def exp_oscillogramms(exp_num):
    test = ProcessSignal(str(exp_num))
    signal_nums_0 = test.read_type_file()['signal_files']
    print(signal_nums_0)
    signal_nums = [signal_nums_0[i] for i in range(len(signal_nums_0)) if signal_nums_0[i] != 'str015.csv']
    voltage_nums = test.read_type_file()['voltage_files']
    for element in zip(signal_nums, voltage_nums):
        one_signal_oscillogramm(element[0], element[1], exp_num)


#exp_oscillogramms(210414)


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

#magnetron_osc(210326)


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


def pl_densities_from_excel(exp_num):
    proc = ProcessSignal(f'{exp_num}')
    all_files_list = os.listdir(proc.exp_file_path)
    file_nums_list = [f'{i:03d}' for i in range(65, 91)]
    pl_d_vals = []
    for file_num in file_nums_list:
        file_name = f'str{file_num}.csv'
        pl_density = proc.read_excel(file_name)['dicts'][file_num]['Ток плазмы, А']
        pl_d_vals.append(pl_density)
    return(pl_d_vals)

#pl_densities_from_excel(210414)


def file_nums_oscillogramms_density_vals(exp_num, file_nums):
    print(file_nums)
    test = ProcessSignal(str(exp_num))
    density_vals = pl_densities_from_excel(exp_num)
    for i, file_num in enumerate(file_nums):
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
        plt.title(f'№ {file_num}, n={density_vals[i]}')
        png_name = test.signal_pics_path / f'ocs_{file_num}_{density_vals[i]}.png'
        fig.savefig(png_name)
        plt.close(fig)


file_nums_oscillogramms_density_vals(210414, [f'{i:03d}' for i in range(65, 91)])


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


