import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

i = Image.open(r'C:\Users\d_Nice\Documents\SignalProcessing\2021\210621\im_4.png')
iar = np.asarray(i)
print(iar)

plt.imshow(iar)
#print(iar)
#plt.show()

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
    print(zero_inds)
    try:
        time_part_max, volt_part_max = time[zero_inds[0]:zero_inds[1]], voltage[zero_inds[0]:zero_inds[1]]
        time_part_min, volt_part_min = time[zero_inds[1]:zero_inds[2]], voltage[zero_inds[1]:zero_inds[2]]

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
    except:
        pass

x, y = [], []
for i in range(iar.shape[1]):
    for j in range(iar.shape[0]):
        point = iar[j][i]
        if point[0] != 0 and point[1] != 0 and point[2] != 0:
            if point[1] < 230 and point[2] < 235 and point[0] == 255:
                print(point)
                if 75 <= i < 560 and 70 < j:
                    x.append(i)
                    y.append(j)
plt.plot(x, y, marker='.', linestyle=' ')
plt.show()

x, y = x[:300:], y[:300:]
plt.plot(x, y)
plt.plot([x[np.argmin(y)], x[np.argmax(y)]], [np.min(y), np.max(y)], marker='.', linestyle=' ')
plt.show()
double_ampl = np.max(y) - np.min(y)

print(double_ampl)


def two_antennas_max_table_new(exp_num, file_nums, start_time=100, end_time=150, time_step=25):
    print(f'Experiment {exp_num}')
    test = ProcessSignal(str(exp_num))
    nums_mtrx = np.reshape(np.asarray(file_nums), (int(len(file_nums) / 2), 2))
    print(nums_mtrx)
    central_freq = 2.714E9
    pts = np.arange(start_time, end_time + time_step, time_step)
    compare_pts = [pts[i] * 1E-9 for i in range(len(pts))]
    col_nums = [i for i in range(3, 3 + len(compare_pts))]
    comment = test.read_excel(f'str{nums_mtrx[0, 0]}.csv')['dicts'][f'{nums_mtrx[0, 0]}']['Комментарий']

    ex_table = excel.Workbook()
    ex_table.create_sheet(title='Integral', index=0)
    sheet = ex_table['Integral']
    sheet['A1'] = exp_num
    sheet['B1'] = comment
    sheet['A2'] = 'Без фильтра'
    sheet['A3'] = 'No (центр/бок)'
    sheet['B3'] = 'n, отн.ед.'
    filt_row = 5 + nums_mtrx.shape[0]
    sheet[f'A{filt_row}'] = 'Фильтрованный'
    sheet[f'A{filt_row + 1}'] = 'No (центр/бок)'
    sheet[f'B{filt_row + 1}'] = 'n, отн.ед.'

    for n in range(len(col_nums)):
        letter = get_column_letter(col_nums[n])
        sheet[f'{letter}3'] = f't_{n + 1} = {np.round(compare_pts[n] / 1e-9, 0)}'
        sheet[f'{letter}{filt_row}'] = f't_{n + 1} = {np.round(compare_pts[n] / 1e-9, 0)} c'

    for j in range(nums_mtrx.shape[0]):
        vals_for_mean = []
        k = 0
        for pt in compare_pts:
            row_ind = 0
            print(nums_mtrx[j, 0], 'pt=', pt)
            file_num = f'{nums_mtrx[j, 0]}'
            file_name = f'str{file_num}.csv'
            data = test.open_file(file_name, reduced=False)
            filt_freq_min, filt_freq_max = central_freq - 30e6, central_freq + 30e6
            t, u = data['time'], data['voltage'] - np.mean(data['voltage'])
            u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
            pl_density = test.read_excel(file_name)['dicts'][file_num]['Ток плазмы, А']
            comment = test.read_excel(file_name)['dicts'][file_num]['Комментарий']
            ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
            t, u = t[ind_mask], u[ind_mask]
            peak_dict = find_zeros(t, u)
            time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict[
                'time_min'], peak_dict['volt_min']
            delta_main = volt_max - volt_min

            u_filt = u_filt[ind_mask]
            filt_dict = find_zeros(t, u_filt)
            time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
            time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
            delta_main_filt = volt_max_filt - volt_min_filt

            cell = sheet.cell(row=j + 4, column=2)
            cell.value = f'{pl_density}'

            cell = sheet.cell(row=j + 1 + filt_row, column=2)
            cell.value = f'{pl_density}'

            else:
            t_shift = 0
            t, u = data['time'] + t_shift, (data['voltage'] - np.mean(data['voltage']))
            u_filt = test.fft_filter(t, u, filt_freq_min, filt_freq_max)
            ind_mask = np.logical_and(t >= pt, t <= pt + 1E-9)
            t, u = t[ind_mask], u[ind_mask]
            peak_dict = find_zeros(t, u)

            time_max, volt_max, time_min, volt_min = peak_dict['time_max'], peak_dict['volt_max'], peak_dict[
                'time_min'], peak_dict['volt_min']
            delta_sub = volt_max - volt_min

            u_filt = u_filt[ind_mask]
            filt_dict = find_zeros(t, u_filt)
            time_max_filt, volt_max_filt = filt_dict['time_max'], filt_dict['volt_max']
            time_min_filt, volt_min_filt = filt_dict['time_min'], filt_dict['volt_min']
            delta_sub_filt = volt_max_filt - volt_min_filt


if int(file_num) < 156:
    u_relat_filt = delta_main_filt / delta_sub_filt
    u_relat = delta_main / delta_sub
else:
    u_relat_filt = delta_sub_filt / delta_main_filt
    u_relat = delta_sub / delta_main

ell = sheet.cell(row=j + 4, column=column)
cell.value = f'{np.round(u_relat, 3)}'
cell = sheet.cell(row=j + 1 + filt_row, column=column)
cell.value = f'{np.round(u_relat_filt, 3)}'

sheet[f'A{j + 4}'] = f'{int(nums_mtrx[j, i]) - 1} / {nums_mtrx[j, i]}'
sheet[f'A{filt_row + j + 1}'] = f'{int(nums_mtrx[j, i]) - 1} / {nums_mtrx[j, i]}'

if pt <= 325:
    vals_for_mean.append(u_relat_filt)

column = col_nums[k]
k = k + 1

mean_filt = np.mean(vals_for_mean)
stand_dev_filt = tstd(vals_for_mean)
cell = sheet.cell(row=j + 1 + filt_row, column=len(pts) + 2)
cell.value = f'{np.round(mean_filt, 3)}'
cell = sheet.cell(row=j + 1 + filt_row, column=len(pts) + 3)
cell.value = f'{np.round(stand_dev_filt, 3)}'

path = test.excel_folder_path / f'{exp_num}_{nums_mtrx[j, i]}_{comment[:6]}_30MHz.xlsx'
ex_table.save(path)
exception_list = [139, 140, 155]
# exception_list = []
# two_antennas_max_table(210707, [f'{i:03d}' for i in range(55, 188) if i not in exception_list])

# file_nums = [f'{i:03d}' for i in range(158, 188) if i not in exception_list]