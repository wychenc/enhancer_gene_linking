
# coding: utf-8

# In[6]:

def run_correlation_method(output_path, rna_mat_file, rna_row_labels, rna_col_labels, dnase_mat_file, dnase_row_labels, dnase_col_labels, chr_name, visualize):         
    import numpy as np
    import pybedtools
    from scipy.stats.stats import pearsonr
    from scipy.stats.stats import spearmanr
    from random import randint
    from matplotlib import pyplot as plt
    from numpy import trapz
    import codecs

    rna_d_new = {}
    rna_row_label = []
    with open(rna_row_labels,'r') as rna_r_f:
        for l in rna_r_f:
            items = l.strip().split('\t')
            rna_row_label.append(items[0] + ',' + items[1] + ',' + items[2])
    rna_r_f.close()

    rna_col_label = []
    with codecs.open(rna_col_labels,'r',encoding = 'ISO-8859-1') as rna_c_f:
        for l in rna_c_f:
            rna_col_label.append(l.strip())
    rna_c_f.close()
    
    # make rna_d_new where keys are tss_loc and cell type, and value is counts
    rna_mat = np.loadtxt(rna_mat_file, dtype=float)
    for r in np.arange(rna_mat.shape[0]):
        for c in np.arange(rna_mat.shape[1]):
            rna_d_new[(rna_row_label[r], rna_col_label[c])] = rna_mat[r,c]

    dna_d_new = {}
    dna_row_label = []
    with open(dnase_row_labels,'r') as dna_r_f:
        for l in dna_r_f:
            items = l.strip().split('\t')
            dna_row_label.append(items[0] + ',' + items[1] + ',' + items[2])
    dna_r_f.close()
    
    dna_col_label = []
    with codecs.open(dnase_col_labels,'r',encoding = 'ISO-8859-1') as dna_c_f:
        for l in dna_c_f:
            dna_col_label.append(l.strip())
    dna_c_f.close()

    # make rna_d_new where keys are tss_loc and cell type, and value is counts
    dna_mat = np.loadtxt(dnase_mat_file, dtype=float)
    for r in np.arange(dna_mat.shape[0]):
        for c in np.arange(dna_mat.shape[1]):
            dna_d_new[(dna_row_label[r], dna_col_label[c])] = dna_mat[r,c]

    # identifies dnase peaks <1MB of each tss
    a = pybedtools.example_bedtool(rna_row_labels)
    b = pybedtools.example_bedtool(dnase_row_labels)
    dnase_within_1MB_of_rna_tss = a.window(b, w=1000000)
    np.savetxt(output_path+'dnase_within_1MB_of_rna_tss.bed',dnase_within_1MB_of_rna_tss, fmt="%s")

    loc_d = {} # dictionary where key is rna loc and value is list of dnase locs < 1 MB
    with open(output_path+'dnase_within_1MB_of_rna_tss.bed','r') as f4: 
        for l4 in f4:
            items = l4.strip().split('\t')
            if tuple(items[0:3]) in loc_d:
                loc_d[tuple(items[0:3])].append(items[3:6]) 
            else:
                loc_d[tuple(items[0:3])] = [items[3:6]]

    # get correlation scores for each rna-dnase pair in loc_list
    score_list_pearson = []
    score_list_spearman = []
    rand_score_list_pearson = []
    rand_score_list_spearman = []
    loc_list = []

    if chr_name == 'all':
        loc_d_sub = loc_d
    else:
        loc_d_sub = {k:v for (k,v) in loc_d.items() if str(chr_name) in k}  

    for rna_loc in loc_d_sub: # for each rna loc 
        dna_loc_list = loc_d_sub[rna_loc] # get list of dna loc <1MB 
        rna_ct_all_samples = [] # list of rna counts across 122 samples for 1 specific rna loc

        for sample in rna_col_label: # for each of 122 samples
            rna_loc_joint = rna_loc[0]+','+rna_loc[1]+','+rna_loc[2]
            rna_ct_all_samples.append(rna_d_new.get((rna_loc_joint, sample)))# got rna counts vector for 1 rna loc

        for dna_loc in dna_loc_list: # for each <1MB dna loc  
            dna_ct_all_samples = [] # list of dna counts across 122 samples for 1 specific dnase loc

            for sample in rna_col_label:
                dna_loc_joint = dna_loc[0]+','+dna_loc[1]+','+dna_loc[2]
                dna_ct_all_samples.append(dna_d_new.get((dna_loc_joint, sample))) # got dnase counts vector for 1 <1MB dnase

            rna_vec = np.asarray(rna_ct_all_samples).astype(np.float)
            dna_vec = np.asarray(dna_ct_all_samples).astype(np.float)

            if np.std(rna_vec)==0 or np.std(dna_vec)==0: # to avoid nan's            
                continue

            loc_list.append((rna_loc, dna_loc)) 
            score_list_pearson.append(pearsonr(rna_vec, dna_vec))
            score_list_spearman.append(spearmanr(rna_vec, dna_vec)) 

            random.shuffle(dna_vec)
            rand_score_list_pearson.append(pearsonr(rna_vec, dna_vec))
            rand_score_list_spearman.append(spearmanr(rna_vec, dna_vec))

    np.savetxt(output_path+chr_name+'_rna_dnase_pair_locations.txt', loc_list, fmt="%s")
    np.savetxt(output_path+chr_name+'_pearson_correlation.txt', score_list_pearson, fmt="%s")
    np.savetxt(output_path+chr_name+'_spearman_correlation.txt', score_list_spearman, fmt="%s")

    if visualize:
        val_pearson = [i[0] for i in score_list_pearson]
        val_spearman = [i[0] for i in score_list_spearman]
        val_rand_pearson = [i[0] for i in rand_score_list_pearson]
        val_rand_spearman = [i[0] for i in rand_score_list_spearman]

        freq_pearson = []
        freq_rand_pearson = []
        freq_spearman = []
        freq_rand_spearman = []
        tick_marks = np.linspace(-1, 1, 101)

        left_tick = -1
        right_tick = -0.98
        for frame in np.arange(100):
            freq_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_pearson))
            freq_rand_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_pearson))
            freq_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_spearman))
            freq_rand_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_spearman))

            left_tick += 0.02
            right_tick += 0.02

        x = tick_marks[0:100]
        y1 = freq_pearson
        y2 = freq_rand_pearson
        y3 = freq_spearman
        y4 = freq_rand_spearman

        # plot frequencies of correlation values
        plt.figure(figsize=(10, 10))
        plt.plot(x,y1,'b',label='pearson real')
        plt.plot(x,y2,'--b',label='pearson null')
        plt.plot(x,y3,'k',label='spearman real')
        plt.plot(x,y4,'--k',label='spearman null')
        plt.legend(loc='upper right')
        plt.title('frequencies of correlation values: real against null')
        plt.xlabel('correlation values')
        plt.ylabel('frequency')
        plt.savefig(output_path+'frequencies_of_correlation_values.png')

        # plot frequencies of absolute pearson correlation values
        freq_pearson = []
        freq_rand_pearson = []
        freq_spearman = []
        freq_rand_spearman = []
        
        tick_marks = np.linspace(-1, 1, 101)

        val_pearson_abs = [abs(number) for number in val_pearson]
        val_rand_pearson_abs = [abs(number) for number in val_rand_pearson]
        val_spearman_abs = [abs(number) for number in val_spearman]
        val_rand_spearman_abs = [abs(number) for number in val_rand_spearman]

        left_tick = -1
        right_tick = -0.98
        for frame in np.arange(100):
            freq_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_pearson_abs))
            freq_rand_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_pearson_abs))
            freq_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_spearman_abs))
            freq_rand_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_spearman_abs))

            left_tick += 0.02
            right_tick += 0.02

        x = tick_marks[0:100]
        y1 = freq_pearson
        y2 = freq_rand_pearson
        y3 = freq_spearman
        y4 = freq_rand_spearman

        plt.figure(figsize=(10, 10))
        plt.plot(x,y1,'b',label='pearson real')
        plt.plot(x,y2,'--b',label='pearson null')
        plt.plot(x,y3,'k',label='spearman real')
        plt.plot(x,y4,'--k',label='spearman null')
        plt.legend(loc='upper right')
        plt.title('frequencies of absolute correlation values: real against null')
        plt.xlabel('correlation values')
        plt.ylabel('frequency')
        plt.savefig(output_path+'frequencies_of_absolute_correlation_values.png')

        # plot area under frequency curve
        freq_pearson = []
        freq_rand_pearson = []
        auc_pearson = []
        auc_rand_pearson = []
        
        freq_spearman = []
        freq_rand_spearman = []
        auc_spearman = []
        auc_rand_spearman = []
        
        tick_marks = np.linspace(-1, 1, 101)
        left_tick = -1
        right_tick = -0.98
        for frame in np.arange(100):
            freq_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_pearson))
            auc_pearson.append(trapz(freq_pearson, dx=0.01))
            freq_rand_pearson.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_pearson))
            auc_rand_pearson.append(trapz(freq_rand_pearson, dx=0.01))
            
            freq_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_spearman))
            auc_spearman.append(trapz(freq_spearman, dx=0.01))
            freq_rand_spearman.append(sum(((i > left_tick) & (i < right_tick))for i in val_rand_spearman))
            auc_rand_spearman.append(trapz(freq_rand_spearman, dx=0.01))
            
            left_tick += 0.02
            right_tick += 0.02

        x = tick_marks[0:100]
        y1 = auc_pearson
        y2 = auc_rand_pearson
        y3 = auc_spearman
        y4 = auc_rand_spearman

        plt.figure(figsize=(10, 10))
        plt.plot(x,y1,'b',label='pearson real')
        plt.plot(x,y2,'--b',label='pearson shuffled')
        plt.plot(x,y3,'k',label='spearman real')
        plt.plot(x,y4,'--k',label='spearman null')
        plt.legend(loc='upper right')
        plt.title('area under frequency curve of correlation values: real against null')
        plt.xlabel('correlation values')
        plt.ylabel('area under frequency curve')
        plt.savefig(output_path+'area_under_frequency_curve_of_correlation_values.png')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='computes pearson and spearman correlations of rna and dnase peaks across different cell samples')
    parser.add_argument('output_path', help='desired path to the output files, ending with a forward slash.')
    parser.add_argument('rna_matrix',help='file tab delimited with rows of rna locations and columns of cell samples')
    parser.add_argument('rna_row_labels',help='bed file of rna tss down each row of rna_matrix')
    parser.add_argument('rna_col_labels',help='list of cell samples across each column of rna_matrix')
    parser.add_argument('dnase_matrix',help='file tab delimited with rows of dnase peak locations and columns of cell samples')  
    parser.add_argument('dnase_row_labels',help='bed file of dnase peaks down each row of dnase_matrix')
    parser.add_argument('dnase_col_labels',help='list of cell samples across each column of dnase_matrix')
    parser.add_argument('-s','--subsample',help='set to "all" for no subsampling, or set to chr+number, e.g. "chr21", to compute correlations for one specific chromosome') 
    parser.add_argument('-v','--visualize',help='visualize various parameter plots of obtained correlation values',action='store_true') 
    args = parser.parse_args()
                        
    #if args.method=='correlation':
    run_correlation_method(args.output_path,args.rna_matrix,args.rna_row_labels,args.rna_col_labels,args.dnase_matrix,args.dnase_row_labels,args.dnase_col_labels,args.subsample,args.visualize)
    
main()

