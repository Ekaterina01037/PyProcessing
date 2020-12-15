import matplotlib.pyplot as plt
from ProcessClass_10 import ProcessSignal
import os
from pathlib import Path
import numpy as np
import csv
from scipy.stats import linregress
import openpyxl as xl

folder_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\phase_shift')

files = os.listdir(folder_path)
def open_file(file_name):
    file_path = folder_path / '{}'.format(file_name)
    t = []
    u = []
    with open(str(file_path)) as File:
        reader = csv.reader(File)
        for row in reader:
            t.append(float(row[3]))
    with open(str(file_path)) as File:
        reader = csv.reader(File)
        for row in reader:
            u.append(float(row[4]))
    dt = t[1] - t[0]
    t = np.array(t)
    u = np.array(u)
    file_dict = {'time': t,
                 'voltage': u,
                 'time_resolution': dt}
    return file_dict


def cut_files():
    for file in files:
        file_data = open_file(file)
        t, u = file_data['time'], file_data['voltage']
        mask_inds = np.logical_and(t >= 120e-9, t <= 650e-9)
        t_cut, u_cut = t[mask_inds], u[mask_inds]
        name_str = file[0:6] + '_cut' + file[6:]
        file_path = folder_path / name_str
        file = open(str(file_path), 'w', newline='')
        with file:
            writer = csv.writer(file)
            for i in range(0, len(t_cut)):
                writer.writerow([0, 0, 0, t_cut[i], u_cut[i]])

#cut_files()

def signal_periods():
    files = os.listdir(folder_path)
    for file in files:
        if 'cut' in file:
            file_data = open_file(file)
            time, voltage = file_data['time'], file_data['voltage']
            voltage_signs = np.sign(voltage)
            m = 0
            l = 0
            start_inds = []
            end_inds = []
            for i in range(len(voltage_signs) - 1):
                if l < 1:
                    if voltage_signs[i] != voltage_signs[i + 1] and voltage_signs[i + 1] != 0:
                        m += 1
                        if m % 2 != 0:
                           start_inds.append(i)
                           end_inds.append(i+1)
            start_volt = voltage[start_inds]
            start_time = time[start_inds]
            end_volt = voltage[end_inds]
            end_time = time[end_inds]
            zero_xs = []
            for i in range(len(start_time)):
                x = np.array([start_time[i], end_time[i]])
                y = np.array([start_volt[i], end_volt[i]])
                k, b, r_value, p_value, std_err = linregress(x, y)
                new_y = 0
                new_x = (new_y - b) / k
                zero_xs.append(new_x)
            periods = np.diff(zero_xs)
            zeros = zero_xs[1::]
            ex_table = xl.Workbook()
            ex_table.create_sheet(title='Период', index=0)
            sheet = ex_table['Период']
            sheet['A1'] = 'Время'
            sheet['B1'] = 'Период'
            for z in range(periods.size):
                cell = sheet.cell(row=z + 2, column=1)
                cell.value = zeros[z]
                cell = sheet.cell(row=z + 2, column=2)
                cell.value = periods[z]
            path = folder_path / f'{file[:6:]}.xlsx'
            ex_table.save(path)
signal_periods()