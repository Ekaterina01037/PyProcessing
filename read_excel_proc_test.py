import openpyxl as xl
from openpyxl.utils import get_column_letter
from pathlib import Path

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