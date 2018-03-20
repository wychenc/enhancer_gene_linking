
# coding: utf-8

# In[28]:

# uses gencode.v19.annotation.gtf 
def make_rna_matrix(output_path):
    import numpy as np
    from sklearn.decomposition import PCA
    from matplotlib import pyplot as plt

    d = {}
    biosample_matched_dnase_set=set()
    gene_id_set=set()
    bios_cell_d={}
    bios_color_d={}
    with open('/mnt/data/integrative/metadata/rna_metadata_2016-12-05.txt','r') as f1: 
        for l1 in f1:
            if ("rsem_quant_file" in l1) or ("UW_Stam" not in l1): # skip header line & only keep data from Stam lab
                continue
            items = l1.strip().split('\t')   
            if items[7] == "" or items[7] == "ENCSR000EII.Primary_CD20+_B_cells_from_peripheral_blood.Duke_Crawford.DNase-seq":
                continue
            biosample_matched_dnase = items[7]
            rsem_quant_file = items[6]
            biosample_matched_dnase_set.add(biosample_matched_dnase)
            bios_cell_d[biosample_matched_dnase] = items[10]
            bios_color_d[biosample_matched_dnase] = items[12]

            with open(rsem_quant_file,'r') as f2: # open rsem link
                for l2 in f2:
                    if "count" in l2:
                        continue
                    link = l2.split('\t')
                    gene_id = link[0] # did not take transcript_id
                    expected_count = link[4] 
                    gene_id_set.add(gene_id)
                    d[(gene_id, biosample_matched_dnase)] = expected_count
    f1.close()
    biosample_matched_dnase_list=list(biosample_matched_dnase_set) # what i'll use to correlate rna and dnase in same sample
    gene_id_list=list(gene_id_set)

    len_bio = len(biosample_matched_dnase_list)
    len_gen = len(gene_id_list)
    rna_matrix = np.empty([len_gen, len_bio], dtype=list)

    for gene_idx in range(len_gen):
        current_gene = gene_id_list[gene_idx]
        for bios_idx in range(len_bio):
            current_bios = biosample_matched_dnase_list[bios_idx]
            for match_idx in range(len_bio):
                current_match = biosample_matched_dnase_list[match_idx]
                current_counts = d.get((current_gene, current_bios))
                rna_matrix[gene_idx, bios_idx] = current_counts

    np.savetxt(output_path+'rna_counts_gene-by-cellline_unsorted_unlabeled.txt', rna_matrix,fmt="%s",delimiter='\t')
    np.savetxt(output_path+'rna_counts_gene-by-cellline_unsorted_column_labels.txt', biosample_matched_dnase_list, fmt="%s",delimiter='\t')
    np.savetxt(output_path+'rna_counts_gene-by-cellline_unsorted_row_labels_genes.txt', gene_id_list, fmt="%s",delimiter='\t')

    # change gene loc to tss using gtf file
    orig_d = {}
    gtf_d = {}
    with open('/srv/scratch/wychen66/get_correlation/rna_outputs/gencode.v19.annotation.gtf', 'r') as f:
        for l in f:
            all_items = l.strip().split('\t')
            if len(all_items)>1 and all_items[2] == 'gene':
                items = l.split(';')[0].split('\t')  
                ids = items[8]
                orig_loc = items[0]+'\t'+items[3]+'\t'+items[4]
                if items[6] == '+': 
                    locs = items[0]+'\t'+items[3]+'\t'+str(int(items[3])+1)
                elif items[6] == '-':
                    locs = items[0]+'\t'+items[4]+'\t'+str(int(items[4])+1)

                orig_d[ids] = orig_loc
                gtf_d[ids] = locs
    
    # load list of rna loc
    rna_tss_loc = []
    rna_set = set()
    with open(output_path+'rna_counts_gene-by-cellline_unsorted_row_labels_genes.txt','r') as file:
        for line in file:
            items = line.strip() 
            rna_set.add(items)
            for key, val in orig_d.items():
                if val == items:
                    rna_tss_loc.append(gtf_d.get(key))
    np.savetxt(output_path+'rna_counts_gene-by-cellline_unsorted_row_labels_tss.txt', rna_tss_loc, fmt="%s")
    
    # these lists are for labeling pca plots
    color_list = []
    tissue_list = []
    for i in np.arange(len(biosample_matched_dnase_list)):
        color = bios_color_d.get(biosample_matched_dnase_list[i])
        color_list.append(color)
        tissue = bios_cell_d.get(biosample_matched_dnase_list[i])
        tissue_list.append(tissue)

    # add top row with eid_matched_danse of each col, and add left column with gene_id of each row
    rna_matrix_labeled = np.insert(rna_matrix, 0, gene_id_list, axis=1) 
    bios_matched_dnase_list_title = [''] + biosample_matched_dnase_list
    rna_matrix_labeled = np.insert(rna_matrix_labeled, 0, [bios_matched_dnase_list_title], axis=0) 
    np.savetxt(output_path + 'rna_counts_gene-by-cellline_unsorted_labeled.txt', rna_matrix_labeled,fmt="%s",delimiter='\t')

    '''
    # first sort rna matrix so that rows of gene_id are ordered according to chromosomal location in bash, then
    sorted_rna_matrix = np.zeros([rna_matrix.shape[0],rna_matrix.shape[1]])
    gene_labels=[]
    with open('/srv/scratch/wychen66/get_correlation/rna_outputs/rna_gene_location_list_with_row_number_sorted.bedGraph', 'r') as f:
        row_count = 0
        for l in f:
            idx = int(l.strip().split('\t')[3])
            gene_labels.append(l.strip().split('\t')[0:3])
            sorted_rna_matrix[row_count,:] = rna_matrix[idx,:]
            row_count += 1
    f.close()
    np.savetxt(output_path+'rna_counts_gene-by-cellline_sorted_unlabeled.txt', sorted_rna_matrix,fmt="%s",delimiter='\t')
    np.savetxt(output_path+'rna_counts_gene-by-cellline_sorted_unlabeled_row_label.txt', gene_labels,fmt="%s",delimiter='\t')
    np.savetxt(output_path+'rna_counts_gene-by-cellline_sorted_column_label.txt', biosample_matched_dnase_list,fmt="%s",delimiter='\t')
    '''
    # below prepares for pca analysis
    sorted_rna_matrix = rna_matrix
    rna_transpose = np.transpose(sorted_rna_matrix)

    # normalize over sequencing depth
    col_sum = np.sum(rna_transpose,axis=1)
    col_sum = col_sum.reshape(col_sum.shape[0],-1)
    rna_transpose_depth_normalized = np.arcsinh(10000000* np.divide(rna_transpose,col_sum))

    # plot PCA of rna matrix
    pca = PCA(n_components=2)
    pca.fit(rna_transpose_depth_normalized) # fit model to data
    rna_pca = pca.transform(rna_transpose_depth_normalized) # reduce dimensions
    print('reduced dimensions:', rna_pca.shape)
    cell_n = np.arange(len(biosample_matched_dnase_list))
    print('pca variance ratio:', pca.explained_variance_ratio_)

    plt.figure(figsize=(10, 10))
    for i, c, cell_label in zip(cell_n, color_list, tissue_list): 
        plt.scatter(rna_pca[i, 0], rna_pca[i, 1], color=c, label=cell_label)
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('rna PCA')
    plt.savefig(output_path+'pca_of_rna_matrix.png')
    '''
    # https://matplotlib.org/examples/statistics/violinplot_demo.html
    pos = np.arange(sorted_rna_matrix.shape[1]) + 1
    data = sorted_rna_matrix
    plt.plot(nrows=5, ncols=5,figsize=(10,10))
    plt.violinplot(data, pos, points=20, widths=1, showmeans=True, showextrema=True, showmedians=True)
    plt.title('rna violin plot')
    plt.show()

    # density plot across ~ 100 samples overlaid on eachother
    x = np.arange(sorted_rna_matrix.shape[0])
    y = []
    for i in np.arange(sorted_rna_matrix.shape[1]):
        y.append(sorted_rna_matrix[:,i])
    y = np.transpose(y)
    plt.figure(figsize=(10, 10))
    plt.plot(x,y,'k',linewidth=0.1)
    plt.title('rna density plot across 123 samples')
    plt.xlabel('57820 genes')
    plt.ylabel('counts')
    plt.show()
    '''
    # make frequency plot of each column of dnase_matrix taking on various values
    freq_matrix = np.zeros([123,100])
    tick_marks = np.linspace(-1, 22889435, 101)

    for sample in np.arange(123):
        left_tick = 0
        right_tick = 1
        for frame in np.arange(100):
            freq_matrix[sample,frame] = ((left_tick < sorted_rna_matrix[:,sample]) & (sorted_rna_matrix[:,sample] < right_tick)).sum()
            left_tick += 1
            right_tick += 1
    x = tick_marks[0:100]
    y = np.transpose(freq_matrix)

    plt.figure(figsize=(10, 10))
    plt.plot(x,y,'k',linewidth=0.1)
    plt.title('frequency plot of each column of rna_matrix taking on various values across 123 samples')
    plt.xlabel('count values in rna_matrix')
    plt.ylabel('frequency of count value')
    plt.savefig(output_path+'frequency_plot_of_rna_matrix_values.png')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='makes ')
    parser.add_argument('output_path', help='the desired path to the output files, ending with a forward slash')
    args = parser.parse_args()
    make_rna_matrix(args.output_path)
    print('Output path: '+args.output_path)
    
main()

