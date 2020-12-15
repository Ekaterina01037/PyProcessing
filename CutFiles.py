from ProcessClass import ProcessSignal

#Классификация файлов сигнал / напряжение
def write_file():
    fft_test = ProcessSignal('191106')
    fft_test.files_classification()
#write_file()

#Чтение классифицированных сигналов из файла
def load_text():
    fft_test = ProcessSignal('190925')
    csv_types = fft_test.read_type_file()
    print(csv_types)
#load_text()

#Обрезка файлов
def write_csv():
    test = ProcessSignal('191106')
    test.reduce_files()
write_csv()

#График обрезанного сигнала
def load_red_csv():
    test = ProcessSignal('191001')
    types = test.read_type_file()
    csv_signals = types['signal_files']
    for signal in csv_signals:
        open_dict = test.open_file(signal, reduced=True)
        t = open_dict['time']
        u = open_dict['voltage']
        plt.plot(t, u)
        plt.show()
#load_red_csv()