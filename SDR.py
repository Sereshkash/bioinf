from multiprocessing import Pool, cpu_count
import os
import seq_split as sq
import gen_seq
import os
import time
import shutil


def mult_iteration(j, N = 3, repeat_set = 6, len_subseq = 410, E_val = 1, number_set_of_subseq = 8,
                    s_name_txt_orig = 'iter4/example_new.txt', name_folder_root = 'SDR1/1/', N1 = 40, number_of_iteration = 10):
    
    name_process = os.getpid()
    #print(j, name_process)
    
    name_folder = os.path.join(name_folder_root, str(name_process))

    if not os.path.isdir(name_folder):
        os.mkdir(name_folder)

    s_local_txt = 's.txt'
    s_name_txt = os.path.join(name_folder, s_local_txt)
    if not os.path.exists(s_name_txt):
        shutil.copy(s_name_txt_orig, s_name_txt)

    s_local_fa = 's.fa'
    s_name_fasta = os.path.join(name_folder, s_local_fa)
    if not os.path.exists(s_name_fasta):
        sq.txt_into_fasta(s_name_txt, s_name_fasta)

    best_res_txt_local = 'best_result.txt'
    best_res_txt = os.path.join(name_folder, best_res_txt_local)
    if not os.path.exists(best_res_txt):
        f = open(best_res_txt, 'w')
        f.write('')
        f.close()

    best_nhmm_posision_local = 'best_posision.fa'
    best_nhmm_posision = os.path.join(name_folder, best_nhmm_posision_local)
    if not os.path.exists(best_nhmm_posision):
        f = open(best_nhmm_posision, 'w')
        f.write('')
        f.close()


    name_set_subseq = os.path.join(name_folder, 'set_subseq.fa')
    name_set_subseq_msa = os.path.join(name_folder, 'set_subseq_msa.fa')
    name_set_subseq_msa_stockh = os.path.join(name_folder, 'set_subseq_msa.sto')
    name_result_nhmmer_info = os.path.join(name_folder, 'result_nhmmer_info.fa')
    name_result_nhmmer_posision = os.path.join(name_folder, 'result_nhmmer_posision.fa')
    name_set_res = os.path.join(name_folder, 'resulst.txt')

    sq.creating_fasta_set(s_name_txt, N, len_subseq, repeat_set, name_set_subseq)
    number_of_subseq_i = 0
    number_of_subseq_i_best = 0
    
    for i in range(number_of_iteration):
        sq.muscle_msa(name_set_subseq, name_set_subseq_msa)
        sq.fasta_into_stockholm(name_set_subseq_msa, name_set_subseq_msa_stockh)

        sq.nhmmer_on_stock_msa(name_set_subseq_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 20)
        sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N1, name_set_subseq)
        number_of_subseq_i = sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_set_res)
        
        if number_of_subseq_i <= N:
            break

    #print(number_of_subseq_i)

    number_of_subseq_i_best = sq.creating_txt_res_nhmmer_set(s_name_txt, best_nhmm_posision, min_len_subseq = 0, max_len_subseq = 1000, N =  10 ** 6, name_res_txt = name_set_res)

    if number_of_subseq_i > number_of_subseq_i_best:
        f = open(name_set_res, 'r')
        l = f.read()
        f.close()

        f = open(best_res_txt, 'w')
        f.write(l)
        f.close()

        f = open(name_result_nhmmer_posision, 'r')
        l = f.read()
        f.close()

        f = open(best_nhmm_posision, 'w')
        f.write(l)
        f.close()

    return name_folder


def increasing_best(name_folder, mean_predict_len_subseq = 400):
    s_local_txt = 's.txt'
    s_name_txt = os.path.join(name_folder, s_local_txt)
    if not os.path.exists(s_name_txt):
        f = open(s_name_txt, 'w')
        f.write('')
        f.close()

    best_nhmm_posision_local = 'best_posision.fa'
    best_nhmm_posision = os.path.join(name_folder, best_nhmm_posision_local)
    if not os.path.exists(best_nhmm_posision):
        f = open(best_nhmm_posision, 'w')
        f.write('')
        f.close()

    best_res_txt_local = 'best_result.txt'
    best_res_txt = os.path.join(name_folder, best_res_txt_local)
    if not os.path.exists(best_res_txt):
        f = open(best_res_txt, 'w')
        f.write('')
        f.close()


    best_nhmm_posision_local_new = 'best_posision_new.fa'
    best_nhmm_posision_new = os.path.join(name_folder, best_nhmm_posision_local_new)
    if not os.path.exists(best_nhmm_posision_new):
        f = open(best_nhmm_posision_new, 'w')
        f.write('')
        f.close()

    number_of_subseq =  sq.creating_txt_res_nhmmer_set(s_name_txt, best_nhmm_posision, min_len_subseq = 0, 
                                                       max_len_subseq = 1000, N =  10 ** 6, name_res_txt = 'best_without_increase_result.txt')
    if number_of_subseq > 20:
        number_of_subseq = sq.increase(s_name_txt, best_nhmm_posision, 
                name_result_nhmmer_posision_new=best_nhmm_posision_new,
                name_folder=name_folder,
                best_name_res_txt = os.path.join(name_folder, 'best_with_increase_result.txt'),
                mean_predict=mean_predict_len_subseq, N=70, diff_variants=5)
    
    return number_of_subseq

if __name__ == '__main__':
    s_name_txt_orig = 'example_1_2.txt'
    mean_predict_len_subseq = 400
    number_set_of_subseq = 50
    N = 3
    repeat_set = 6
    len_subseq = 410
    E_val = 20
    N1 = 40
    number_of_iteration = 10

    time_start = time.time()

    name_folder_root = 'SDR1'

    if not os.path.isdir(name_folder_root):
        os.mkdir(name_folder_root)

    name_folder = os.path.join(name_folder_root, str(int(time_start)))

    if not os.path.isdir(name_folder):
        os.mkdir(name_folder)

    print('S')

    def iteration_run(x):
        return mult_iteration(x, N = 3, repeat_set = 6, len_subseq = 410,
                               E_val = 20, number_set_of_subseq = number_set_of_subseq, 
                               s_name_txt_orig = s_name_txt_orig, 
                               name_folder_root = name_folder, N1 = 40, number_of_iteration = 10)
    
    
    def iteration_run_increase(x):
        return increasing_best(x, mean_predict_len_subseq = mean_predict_len_subseq)
    
    with Pool(processes=cpu_count()) as pool:
        values = [j + 1 for j in range(number_set_of_subseq)]
        results_folder = pool.map(iteration_run, values)
    time_finish_search = time.time()
    
    results_folder_set = set(results_folder)
    results_folder_diff = []
    for res_folder in results_folder_set:
        results_folder_diff.append(res_folder)

    print(round((time_finish_search - time_start) / 60, 2))

    with Pool(processes=cpu_count()) as pool2:
        results_N = pool2.map(iteration_run_increase, results_folder_diff)

    max_N = max(results_N)
    max_N_folder = results_folder[results_N.index(max_N)]

    print('N:', *results_N)

    time_finish = time.time()
    time_res_seach = time_finish_search - time_start
    time_res_increase = time_finish - time_finish_search
    s = 'N = {}, max_N_folder = {}, time_search = {}, time_increase = {}'.format(max_N, max_N_folder, round(time_res_seach / 60, 2),  round(time_res_increase / 60, 2))
    print(s)


