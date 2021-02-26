import numpy as np
import openpyxl as xl
from pathlib import Path
import os


def read_excel_spectra_data(exp_num):
    excel_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\{}\Excel'.format(exp_num))
    noise_base_folder_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\{}\Pictures\FFT\Спектр шумового пъедестала'.format(exp_num))
    reb_noise_folder_path = Path(r'C:\Users\d_Nice\Documents\SignalProcessing\2020\{}\Pictures\FFT\Спектр шумов РЭП'.format(exp_num))
    path_list = [noise_base_folder_path, reb_noise_folder_path]
    for path in path_list:
        file_list = os.listdir(path)
        for file in file_list:
            if 'xlsx' in file:
                wb = xl.load_workbook(path/file)
                sheet = wb['Full FFT']
                rows = sheet.max_row

                if path == noise_base_folder_path:
                    nums_base = np.asarray([int(sheet[f'A{i}'].value) for i in range(2, rows+1)])
                    n_s_base, freqs_base = np.asarray([float(sheet[f'C{i}'].value) for i in range(2, rows+1)]), np.asarray([float(sheet[f'B{i}'].value) for i in range(2, rows+1)])
                    sorted_inds = np.argsort(n_s_base)
                    nums_base_sorted, n_s_base_sorted, freqs_base_sorted = nums_base[sorted_inds], n_s_base[sorted_inds], freqs_base[sorted_inds]
                else:
                    nums_reb = np.asarray([int(sheet[f'A{i}'].value) for i in range(2, rows+1)])
                    n_s_reb, freqs_reb = np.asarray([float(sheet[f'C{i}'].value) for i in range(2, rows+1)]), np.asarray([float(sheet[f'B{i}'].value) for i in range(2, rows+1)])
                    sorted_inds = np.argsort(n_s_reb)
                    nums_reb_sorted, n_s_reb_sorted, freqs_reb_sorted = nums_reb[sorted_inds], n_s_reb[sorted_inds], freqs_reb[sorted_inds]

    new_table = xl.Workbook()
    new_table.create_sheet(title='Mean_freq_noise', index=0)
    sheet = new_table['Mean_freq_noise']
    sheet['A1'] = 'Номер(шум)'
    sheet['B1'] = 'Плотность плазмы, отн.ед.'
    sheet['C1'] = 'Средняя частота, ГГц'

    sheet['E1'] = 'Номер'
    sheet['F1'] = 'Плотность плазмы, отн.ед.'
    sheet['G1'] = 'Средняя частота, ГГц'
    for z in range(len(nums_reb_sorted)):
        cell = sheet.cell(row=z + 2, column=1)
        cell.value = int(nums_reb_sorted[z])
        cell = sheet.cell(row=z + 2, column=2)
        cell.value = n_s_reb_sorted[z]
        cell = sheet.cell(row=z + 2, column=3)
        cell.value = freqs_reb_sorted[z]

    for v in range(len(nums_base_sorted)):
        cell = sheet.cell(row=v + 2, column=5)
        cell.value = int(nums_base_sorted[v])
        cell = sheet.cell(row=v + 2, column=6)
        cell.value = n_s_base_sorted[v]
        cell = sheet.cell(row=v + 2, column=7)
        cell.value = freqs_base_sorted[v]

    new_table_path = excel_path/f'mean_freq_{exp_num}.xlsx'
    new_table.save(new_table_path)


read_excel_spectra_data(210225)