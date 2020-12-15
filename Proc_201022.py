import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from ProcessClass_10 import ProcessSignal
import matplotlib.pyplot as plt
import os

test = ProcessSignal('201026')
csv_signals = os.listdir(test.exp_file_path)


def file_classification(signals):
    magnetron_files = []
    signal_files = []
    voltage_files = []
    for signal in signals:
        if 'csv' in signal:
            file = test.open_file(signal, reduced=False)
            t, u = file['time'], file['voltage']
            dt = np.abs(t[1] - t[0])
            env_curv = test.signal_envelope(t, u, dt)
            u_env_curv = env_curv['env_voltage']
            part_env_data = test.useful_part(t, u, dt, max_part=0.45)
            try:
                u_env, t_env = np.asarray(part_env_data['signal_voltage']), np.asarray(part_env_data['signal_time'])
                delta_time = t_env[-1] - t_env[0]
                if delta_time > 400e-9 and np.mean(u_env_curv) > 1:
                    magnetron_files.append(signal)
                    plt.plot(t, u)
                    plt.title(signal[3:6])
                    plt.show()
                elif np.mean(u) < 0:
                    voltage_files.append(signal)
                else:
                    signal_files.append(signal)
            except:
                pass
    classif_dict = {'magnetron': magnetron_files,
                    'signals': signal_files,
                    'voltage': voltage_files}
    return classif_dict

file_classification(csv_signals)

def integrals_26():
    signal_files = file_classification(csv_signals)['signals']
    for file in signal_files:
        num = file[3:6]
        if int(num) > 27:
            pl_density = test.read_excel(file)['dicts'][num]['Ток плазмы, А']
            data = test.open_file(file, reduced=False)
            t, u = data['time'], data['voltage']
            dt = np.abs(t[1] - t[0])

            env = test.signal_envelope(t, u, dt)
            env_t = np.array(env['env_time'])
            env_u = np.array(env['env_voltage'])
            max_env_u = np.max(env_u)
            ind_lim = env_u > 0.2 * max_env_u
            t_lim_start = env_t[ind_lim][0]
            print(t_lim_start)
            t_lim_stop = t_lim_start + 262e-9

            time_inds = np.logical_and(t > t_lim_start, t < t_lim_stop)
            t_calc = t[time_inds]
            u_calc = u[time_inds]
            plt.plot(t, u)
            plt.plot(t_calc, u_calc)
            plt.show()

#integrals_26()


def magnetron_spectra(magnetron_files, signal_files):
    for m_file, s_file in zip(magnetron_files, signal_files):
        num = s_file[3:6]
        m_num = m_file[3:6]
        pl_density = test.read_excel(s_file)['dicts'][num]['Ток плазмы, А']
        s_data = test.open_file(s_file, reduced=False)
        t, u = s_data['time'], s_data['voltage']
        dt = np.abs(t[1] - t[0])

        cut_data = test.useful_part(t, u, dt, max_part=0.2)
        t_cut, u_cut = cut_data['signal_time'], cut_data['signal_voltage']



        m_data = test.open_file(m_file, reduced=False)
        m_t, m_u = m_data['time'], m_data['voltage']
        m_ind_cut = np.logical_and(m_t > t_cut[0], m_t < t_cut[-1])
        t_m_cut = m_t[m_ind_cut]
        u_m_cut = m_u[m_ind_cut]

        #calculating_the_spectra
        signal_fft = test.fft_amplitude(t_cut, u_cut)
        s_amps, s_freqs = signal_fft['amplitude'], signal_fft['frequency']

        magnetron_fft = test.fft_amplitude(t_m_cut, u_m_cut)
        m_amps, m_freqs = magnetron_fft['amplitude'], magnetron_fft['frequency']

        #plotting

        fig = plt.figure(num=1, dpi=150)
        ax = fig.add_subplot(111)
        print('Creating a picture...')
        line1, = ax.plot(s_freqs, s_amps, linewidth=0.7, color='mediumseagreen')
        #line2, = ax.plot(m_freqs, m_amps, linewidth=1.5, color='blue')
        line1.set_label('Сигнал на выходе')
        #line2.set_label('Сигнал на входе')
        ax.set_xlabel(r'$Частота, ГГц$', fontsize=14, fontweight='black')
        ax.set_ylabel(r'$Амплитуда$', fontsize=14, fontweight='black')
        ax.set_xlim(left=0, right=4e9)
        ax.set_ylim(bottom=0)
        ax.legend()
        ax.grid(which='both', axis='both')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('№ {}, n = {}'.format(s_file[3:6], pl_density))
        png_name = test.fft_pics_path / 'spectrum_{}'.format(num)
        fig.savefig(png_name)
        plt.close(fig)
        '''
        m_fig = plt.figure(num=1, dpi=150)
        m_ax = m_fig.add_subplot(111)
        print('Creating a picture...')
        m_line1, = m_ax.plot(m_freqs, m_amps, linewidth=1.5, color='blue')
        m_line1.set_label('Сигнал на входе')
        m_ax.set_xlabel(r'$Частота, ГГц$', fontsize=14, fontweight='black')
        m_ax.set_ylabel(r'$Амплитуда$', fontsize=14, fontweight='black')
        m_ax.set_xlim(left=2.68e9, right=2.74e9)
        m_ax.set_ylim(bottom=0)
        m_ax.grid(which='both', axis='both')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        m_ax.legend()
        plt.title('№ {}, n = {}'.format(m_file[3:6], pl_density))
        png_name = test.fft_pics_path / 'magnetron_peak_{}'.format(m_num)
        m_fig.savefig(png_name)
        plt.close(m_fig)
        #plt.plot(m_t, m_u)
        #plt.plot(t_m_cut, u_m_cut)
        #plt.show()
        '''
#magnetron_spectra(magnetron_files, signal_files)