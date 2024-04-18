from multiprocessing import Pool, cpu_count
import os
import seq_split as sq
import gen_seq
import os
import time


def mult_iteration(j):
    
    #print(j)

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

    for i in range(10):
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
        if len(s) // 2 >= 1000 or len(s) // 2 <= 0:
            break
        
    i = 10
    name_set_res = name_folder + 'name_set_res.txt'
    name_result_nhmmer_info = name_folder + 'name_result_nhmmer_info_' + str(i+1) + '.fa'
    sq.muscle_msa(name_set_subseq0, name_set_subseq0_msa)
    sq.fasta_into_stockholm(name_set_subseq0_msa, name_set_subseq0_msa_stockh)
    sq.nhmmer_on_stock_msa(name_set_subseq0_msa_stockh, s_name_fasta, name_result_nhmmer_info, name_result_nhmmer_posision, 20)
    sq.creating_fasta_new_nhmmer_set(s_name_txt, name_result_nhmmer_posision, 10 ** 6, name_set_res)
    sq.creating_txt_res_nhmmer_set(s_name_txt, name_result_nhmmer_posision, min_len_subseq = 0, max_len_subseq = 500, N =  10 ** 6, name_res_txt = name_set_res)

if __name__ == '__main__':
    for l in range(10, 22):
        x_l = l / 20
        time_start = time.time()

        gen_seq_S = gen_seq.generation_seq(10**7, 10**3, 400, x = x_l, alphabet = 'GTAC')
        s_name_txt = 'example_new.txt'
        f = open(s_name_txt, 'w')
        f.write(gen_seq_S)
        f.close()
        s_name_fasta = 'example_new.fa'

        sq.txt_into_fasta(s_name_txt, s_name_fasta)

        N = 3
        repeat_set = 6
        len_subseq = 410
        E_val = 50
        number_set_of_subseq = 400

        #for i in range(nubmer_set_of_subseq):

        j = 1

        time_start = time.time()
        with Pool(processes=cpu_count()) as pool:
            values = [j + 1 for j in range(number_set_of_subseq)]
            results = pool.map(mult_iteration, values)

        f = open('data_res_comp_with_repeat1.txt', 'a')
        f.write('x = {}:\n'.format(x_l))
        f.close()

        max_N = 0
        max_N_index = 0
        set_N = ''
        for j in range(1, number_set_of_subseq):
            name_folder = 'iter/' + 'it' + str(j) + '/'
            name_set_res = name_folder + 'name_set_res.txt'
            f = open(name_set_res, 'r')
            s = f.readlines()
            f.close()
            N_j = len(s) // 2
            set_N += str(N_j) + ' '
            if N_j > max_N:
                max_N = N_j
                max_N_index = j
        
        f = open('data_res_comp_with_repeat1.txt', 'a')
        f.write('{}\n'.format(set_N))
        f.close()


        time_finish = time.time()
        time_res = time_finish - time_start
        s = 'x = {}, N = {}, time = {}\n'.format(x_l, max_N, round(time_res / 60, 2))
        f = open('res_comp_with_repeat1.txt', 'a')
        f.write(s)
        f.close()
        print('x_i = {}'.format(x_l))


