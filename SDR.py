from multiprocessing import Pool, cpu_count
import os
import seq_split as sq
import gen_seq
import os
import time
import shutil
import argparse

def mult_iteration(j, N = 3, repeat_set = 6, len_subseq = 410, E_val = 1, number_set_of_subseq = 8,
                    s_name_txt_orig = 'iter4/example_new.txt', name_folder_root = 'SDR1/1/', N1 = 40, number_of_iteration = 10):
    #j - номер начального набора
    #N - количество различных подпоследовательностей в начальном наборе
    #repeat_set - количество повторов подпоследовательностей в начальном наборе 
    #len_subseq - длина подпоследовательностей в начальном наборе
    #E-val - E-value
    #number_set_of_subseq - количество начальных наборов
    #s_name_txt_orig - название файла с исходной последовательностью
    #name_folder_root - название папки, где сохранять результат
    #N1 - количество подпоследовательностей, которые нужно брать на шаге улучшения
    #number_of_iteration - количество итераций улучшения

    #создание папок под процесс
    name_process = os.getpid()
    
    name_folder = os.path.join(name_folder_root, str(name_process))

    if not os.path.isdir(name_folder):
        os.mkdir(name_folder)

    #создание исходной подпоследовательности в папке
    s_local_txt = 's.txt'
    s_name_txt = os.path.join(name_folder, s_local_txt)
    if not os.path.exists(s_name_txt):
        shutil.copy(s_name_txt_orig, s_name_txt)

    #создание исходной подпоследовательности в папке в формате FASTA
    s_local_fa = 's.fa'
    s_name_fasta = os.path.join(name_folder, s_local_fa)
    if not os.path.exists(s_name_fasta):
        sq.txt_into_fasta(s_name_txt, s_name_fasta)

    #создание файла для записи лучшего по процессу результата
    best_res_txt_local = 'best_result.txt'
    best_res_txt = os.path.join(name_folder, best_res_txt_local)
    if not os.path.exists(best_res_txt):
        f = open(best_res_txt, 'w')
        f.write('')
        f.close()

    #создание файла для записи позиций лучшего по процессу результата
    best_nhmm_posision_local = 'best_posision.fa'
    best_nhmm_posision = os.path.join(name_folder, best_nhmm_posision_local)
    if not os.path.exists(best_nhmm_posision):
        f = open(best_nhmm_posision, 'w')
        f.write('')
        f.close()

    #создание файла для записи текущего набора подпоследовательностей
    name_set_subseq = os.path.join(name_folder, 'set_subseq.fa')
    #создание файла для записи выравненного текущего набора подпоследовательностей
    name_set_subseq_msa = os.path.join(name_folder, 'set_subseq_msa.fa')
    #создание файла для записи текущего набора подпоследовательностей в формате Stockholm
    name_set_subseq_msa_stockh = os.path.join(name_folder, 'set_subseq_msa.sto')
    #создание файла для записи текущей информации о работе nHMMER
    name_result_nhmmer_info = os.path.join(name_folder, 'result_nhmmer_info.fa')
    #создание файла для записи текущей информации о позициях, найденных nHMMER 
    name_result_nhmmer_posision = os.path.join(name_folder, 'result_nhmmer_posision.fa')
    #создание файла для записи текущих найденных подпоследовательностей в формате txt
    name_set_res = os.path.join(name_folder, 'resulst.txt')

    #создание начального набора в формате FASTA
    sq.creating_fasta_set(s_name_txt, N, len_subseq, repeat_set, name_set_subseq)

    #запись количества лучших по потоку и текущих найденных наборов
    number_of_subseq_i = 0
    number_of_subseq_i_best = 0
    
    #запуск итерационной процедуры улучшения
    for i in range(number_of_iteration):
        #получение множественного выравнивания набора
        sq.muscle_msa(name_set_subseq, name_set_subseq_msa)
        #перевод его в формат Stockholm
        sq.fasta_into_stockholm(name_set_subseq_msa, name_set_subseq_msa_stockh)

        #запуск поиска наиболее схожих с текущем набором подпоследовательностей исходной строки
        sq.nhmmer_on_stock_msa(name_set_subseq_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 20)
        #создание нового набора в формате FASTA
        sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N1, name_set_subseq)
        #получение количества найденных подпоследовательностей
        number_of_subseq_i = sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_set_res)
        
        #критерий прекращения итерации улучшения
        if number_of_subseq_i <= N:
            break

    #получение лучшего числа найденных по потоку подпоследовательностей 
    number_of_subseq_i_best = sq.creating_txt_res_nhmmer_set(s_name_txt, best_nhmm_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = best_res_txt)

    #измениние лучшего по потоку значения в случае большем найденном количестве подпоследовательностей
    if number_of_subseq_i > number_of_subseq_i_best:
        #изменение лучшего набора
        f = open(name_set_res, 'r')
        l = f.read()
        f.close()

        f = open(best_res_txt, 'w')
        f.write(l)
        f.close()

        #изменение лучших позиций
        f = open(name_result_nhmmer_posision, 'r')
        l = f.read()
        f.close()

        f = open(best_nhmm_posision, 'w')
        f.write(l)
        f.close()

    #возвращение названия папки по потоку
    return name_folder


def increasing_best(name_folder, mean_predict_len_subseq = 400, working = True):
    #применяет операцию увеличения длин до определенного среднего для лучшего по потоку значения
    
    #получение файла с исходной последовательностью
    s_local_txt = 's.txt'
    s_name_txt = os.path.join(name_folder, s_local_txt)
    if not os.path.exists(s_name_txt):
        f = open(s_name_txt, 'w')
        f.write('')
        f.close()

    #получение файла с лучшими позициями по потоку
    best_nhmm_posision_local = 'best_posision.fa'
    best_nhmm_posision = os.path.join(name_folder, best_nhmm_posision_local)
    if not os.path.exists(best_nhmm_posision):
        f = open(best_nhmm_posision, 'w')
        f.write('')
        f.close()

    #получение файла с лучшим результатом по потоку
    best_res_txt_local = 'best_result.txt'
    best_res_txt = os.path.join(name_folder, best_res_txt_local)
    if not os.path.exists(best_res_txt):
        f = open(best_res_txt, 'w')
        f.write('')
        f.close()

    #создание файля для новых лучший позиций
    best_nhmm_posision_local_new = 'best_posision_new.fa'
    best_nhmm_posision_new = os.path.join(name_folder, best_nhmm_posision_local_new)
    if not os.path.exists(best_nhmm_posision_new):
        f = open(best_nhmm_posision_new, 'w')
        f.write('')
        f.close()

    #получение лучшего количества найденных повторов по потоку
    number_of_subseq =  sq.creating_txt_res_nhmmer_set(s_name_txt, best_nhmm_posision, min_len_subseq = 0, 
                                                       max_len_subseq = 1000, N =  10 ** 6, name_res_txt = os.path.join(name_folder, 'best_without_increase_result.txt'))
    #запуск итерации улучшения
    if number_of_subseq >= 6 and working:
        number_of_subseq = sq.increase(s_name_txt, best_nhmm_posision, 
                name_result_nhmmer_posision_new=best_nhmm_posision_new,
                name_folder=name_folder,
                best_name_res_txt = os.path.join(name_folder, 'best_with_increase_result.txt'),
                mean_predict=mean_predict_len_subseq, N=70, diff_variants=5)
    
    #возвращает новое улучшенное количество значений
    return number_of_subseq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Ping script")

    parser.add_argument('-f', '--filename', dest='s_name', required=True, help='Name of file with sequence')
    parser.add_argument('--ml', '--mean_len', dest='mean_len', default=400, type=int, required=False, help='Mean len of repeat sequence')
    parser.add_argument('--ns', '--number_subseq', dest='number_subseq', default=100, type=int, required=False, help='Number of different search')
    parser.add_argument('-N', '--number_diff_subseq', dest='N', default=3, type=int, required=False, help='Number of different initial subsequences')
    parser.add_argument('--rp', '--repeat_subseq', dest='repeat_subseq', default=6, type=int, required=False, help='Number of repeat initial subsequences')
    parser.add_argument('-l', '--len_subseq', dest='len_subseq', default=410, type=int, required=False,  help='Len initial subsequence')
    parser.add_argument('-E', '--E_value', dest='E', default=0.01, type=float, required=False,  help='E-value')
    parser.add_argument('--ni', '--number_iteration', dest='number_iteration', default=10, type=int, required=False,  help='Number of iteration in one proccess')
    parser.add_argument('--N_search', dest='N_search', default=40, type=int, required=False,  help='Number of selected sequences')
    parser.add_argument('--increase_len', dest='increase_len', default=True, type=bool, required=False,  help='Increase len, bool')
    parser.add_argument('--nf', '--name_folder', dest='name_folder_root', default='SDR', type=str, required=False,  help='Name folder for files')
    parser.add_argument('-o', '--name_file_out', dest='name_file_out', default='-', type=str, required=False,  help='Name file for result')
    
    args = parser.parse_args()

    #запуск таймера времни
    time_start = time.time()

    #создание общей папки
    name_folder_root = args.name_folder_root

    if not os.path.isdir(name_folder_root):
        os.mkdir(name_folder_root)

    #создание общей папки, завязанной на время запуска
    name_folder = os.path.join(name_folder_root, str(int(time_start)))

    if not os.path.isdir(name_folder):
        os.mkdir(name_folder)
  
    #функция, необходимая для запуска основной программы по средствам Pool.map, принимающей функцию с 1 аргументом
    #запуск основной программы с переданными начальными аргументами
    def iteration_run(x):
        return mult_iteration(x, N = args.N, repeat_set = args.repeat_subseq, len_subseq = args.len_subseq,
                               E_val = args.E, number_set_of_subseq = args.number_subseq, 
                               s_name_txt_orig = args.s_name, 
                               name_folder_root = name_folder, N1 = args.N_search, number_of_iteration = args.number_iteration)
    
    #функция, необходимая для запуска функции увеличения длин по средствам Pool.map, принимающей функцию с 1 аргументом
    def iteration_run_increase(x):
        return increasing_best(x, mean_predict_len_subseq = args.mean_len, working = args.increase_len)
    
    #создание процессов по количеству доступных ядер и запуск основной программы
    with Pool(processes=cpu_count()) as pool:
        values = [j + 1 for j in range(args.number_subseq)]
        #results_folder - названия папок, созданных для отдельных процессов
        results_folder = pool.map(iteration_run, values)
    
    #отсчет времени работы основной программ поиска
    time_finish_search = time.time()
    
    #создание массива с названиями папок без повтороений
    results_folder_set = set(results_folder)
    results_folder_diff = []
    for res_folder in results_folder_set:
        results_folder_diff.append(res_folder)

    #запуск по процессам функций увеличения длин подпоследовательностей
    with Pool(processes=cpu_count()) as pool2:
        results_N = pool2.map(iteration_run_increase, results_folder_diff)

    #получение максимального количества найденных подпоследовательностей и определение процесса, в котором он найден
    max_N = max(results_N)
    max_N_folder = results_folder_diff[results_N.index(max_N)]

    #получение времени работы алгоритма поиска и алгоритма увеличения длин
    time_finish = time.time()
    time_res_seach = time_finish_search - time_start
    time_res_increase = time_finish - time_finish_search

    #print('N:', *results_N)
    
    #запись в файл или в основной поток вывода результата о времени, количества найденных подпоследовательностей и названия папки, где лежит ответ
    s = 'N = {}, max_N_folder = {}, time_search = {}, time_increase = {}'.format(max_N, max_N_folder, round(time_res_seach / 60, 2),  round(time_res_increase / 60, 2))
    if args.name_file_out == '-':
        print(s)
    else:
        f = open(args.name_file_out)
        f.write(s)
        f.close()


