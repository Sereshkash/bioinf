from multiprocessing import Pool, cpu_count
import os
import seq_split as sq
import gen_seq
import os
import time


def mult_iteration(j, N = 3, repeat_set = 6, len_subseq = 410, E_val = 1, number_set_of_subseq = 8):
    if not os.path.isdir('iter/'):
        os.mkdir('iter/')

    name_folder = 'iter/' + 'it' + str(j) + '/'

    if not os.path.isdir(name_folder):
        os.mkdir(name_folder)

    name_set_subseq0 = name_folder + 'example_subseq.fa'
    name_set_subseq0_msa = name_folder + 'example_subseq_msa.fa'
    name_set_subseq0_msa_stockh = name_folder + 'example_subseq_msa.sto'
    name_result_nhmmer_info_0 = name_folder + 'name_result_nhmmer_info_0.fa'
    name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info.fa'
    
    name_result_nhmmer_posision = name_folder + 'name_result_nhmmer_posision.fa'
    name_set_res = name_folder + 'name_set_res.txt'

    N1 = 40

    sq.creating_fasta_set(s_name_txt, N, len_subseq, repeat_set, name_set_subseq0)
    sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info_0, name_result_nhmmer_posision, 50)
    sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N1, name_set_subseq0)

    for i in range(9):
        #sq.creating_fasta_set(s_name_txt, N, len_subseq, name_set_subseq0)
        name_set_res = name_folder + 'name_set_res_' + str(i+1) + '.txt'
        name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info_' + str(i+1) + '.fa'
        sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
        sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
        sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 50)
        sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N1, name_set_subseq0)
        sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 500, N =  10 ** 6, name_res_txt = name_set_res)
        f1 = open(name_set_res)
        s = f1.readlines()
        f1.close()
        if len(s) // 2 >= 1000 or len(s) // 2 <= N:
            break
        
    i = 10
    name_set_res = name_folder + 'name_set_res.txt'
    name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info_' + str(i+1) + '.fa'
    sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 20)
    sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 500, N =  10 ** 6, name_res_txt = name_set_res)
    f1 = open(name_set_res)
    s = f1.readlines()
    f1.close()
    if len(s) // 2 >= 10:
        sq.increase_len_subseq(s_name_txt, name_result_nhmmer_posision, name_result_nhmmer_posision, left_step = 0, right_step = 0, min_len_subseq =  20, max_len_subseq = 1000)
        sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, N1, name_set_subseq0)
        sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 500, N =  10 ** 6, name_res_txt = name_set_res)
        sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
        sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
        sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 50)
        sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 700, N =  10 ** 6, name_res_txt = name_set_res)
    if len(s) // 2 >= 400:
        print('len(s) = {}, iter = {}'.format(len(s) // 2, j))
            
        
if __name__ == '__main__':
    time_start = time.time()
    x_l = 1.2
    gen_seq_S = gen_seq.generation_seq(10**7, 10**3, 400, x_l, alphabet = 'GTAC')
    s_name_txt = 'iter3/example_new.txt'
    f = open(s_name_txt, 'w')
    f.write(gen_seq_S)
    f.close()
    s_name_fasta = 'iter3/example_new.fa'

    sq.txt_into_fasta(s_name_txt, s_name_fasta)

    N = 3
    repeat_set = 6
    len_subseq = 410
    E_val = 1
    number_set_of_subseq = 8

    #for i in range(nubmer_set_of_subseq):

    j = 1


    with Pool(processes=cpu_count()) as pool:
        values = [j + 1 for j in range(number_set_of_subseq)]
        results = pool.map(mult_iteration, values)

    max_N = 0
    max_N_index = 0
    for j in range(1, number_set_of_subseq + 1):
        name_folder = 'iter/' + 'it' + str(j) + '/'
        name_set_res = name_folder + 'name_set_res.txt'
        f = open(name_set_res)
        s = f.readlines()
        f.close()
        N_j = len(s) // 2
        print(N_j)
        if N_j > max_N:
            max_N = N_j
            max_N_index = j


    time_finish = time.time()
    time_res = time_finish - time_start
    s = 'x = {}, N = {}, max_N_index = {}, time = {}'.format(x_l, max_N, max_N_index, round(time_res / 60, 2))
    print(s)


