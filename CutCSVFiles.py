from ProcessClass import ProcessSignal

#классифицирует файлы сигнал/напряжение и записывает их в таблицу:
def write_file():
    fft_test = ProcessSignal('191126')
    fft_test.files_classification()
#write_file()

#загружает номера классифицированных файлов из таблицы:
def load_text():
    fft_test = ProcessSignal('190925')
    txt_path = fft_test.exp_file_path / 'types.txt'
    csv_types = fft_test.read_type_file()
    print(csv_types)
#load_text()

#режет исходные CSV
def write_csv():
    test = ProcessSignal('191120')
    test.reduce_files()
#write_csv()