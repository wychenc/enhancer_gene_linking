
# coding: utf-8

# In[2]:

def make_dnase_matrix_post_shell_script(input_path, output_path)
    import numpy as np
    from sklearn.decomposition import PCA
    from matplotlib import pyplot as plt

    dnase_matrix = np.loadtxt(input_path+'final_dnase_matrix.txt',dtype=float)
    dnase_transpose = np.transpose(dnase_matrix)

    dnase_bios_labels=[]
    with open(input_path+'folder_name.txt','r') as f1: 
        for l1 in f1:
            items = l1.strip().split('\t')
            dnase_bios_labels.append(items[0])
    f1.close()
    np.savetxt(output_path+'final_dnase_matrix_column_labels.txt', dnase_bios_labels, fmt="%s")    

    # pca analysis, normalized over sequencing depth
    col_sum = np.sum(dnase_transpose,axis=1)
    col_sum = col_sum.reshape(col_sum.shape[0],-1) 
    dnase_transpose_depth_normalized = np.arcsinh(1000000* np.divide(dnase_transpose,col_sum))

    pca = PCA(n_components=2)
    pca.fit(dnase_transpose_depth_normalized) # fit model to data
    dnase_pca = pca.transform(dnase_transpose_depth_normalized) # reduce dimensions
    print('reduced dimensions:', dnase_pca.shape)

    # make vector of labels
    pca_labels =[]
    with open(input_path+'pca_labels.txt','r') as f1: 
        for l1 in f1:
            l1=l1.strip().replace('\n','')
            pca_labels.append(l1)
    f1.close()

    # make vector of colors
    color_list=[]
    with open(input_path+'pca_color_codes.txt','r') as f2: 
        for l2 in f2:
            l2=l2.strip().replace('\n','')
            color_list.append(l2)
    f2.close()

    # plot PCA output
    cell_n = np.arange(len(pca_labels))
    print('pca variance ratio:', pca.explained_variance_ratio_)

    plt.figure(figsize=(10, 10))
    for i, c, cell_label in zip(cell_n, color_list, pca_labels): 
        plt.scatter(dnase_pca[i, 0], dnase_pca[i, 1], color=c, label=cell_label)

    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('dnase PCA')
    plt.savefig(output_path+'pca_of_dnase_matrix.png')

    '''
    # viollin plot: https://matplotlib.org/examples/statistics/violinplot_demo.html
    pos = np.arange(dnase_matrix.shape[1]) + 1
    data = dnase_matrix
    plt.plot(nrows=5, ncols=5,figsize=(10,10))
    plt.violinplot(data, pos, points=20, widths=1, showmeans=True, showextrema=True, showmedians=True)
    plt.title('dnase violin plot')
    plt.show()

    # density plot across ~ 100 samples overlaid on eachother
    x = np.arange(dnase_matrix.shape[0])
    y = []
    for i in np.arange(dnase_matrix.shape[1]):
        y.append(dnase_matrix[:,i])
    y = np.transpose(y)
    plt.figure(figsize=(15, 15))
    plt.plot(x,y,'k',linewidth=0.1)
    plt.title('dnase density plot across 122 samples')
    plt.xlabel('762379 dnase peaks')
    plt.ylabel('counts')
    plt.show()
    '''
    # make frequency plot of each column of dnase_matrix taking on various values
    freq_matrix = np.zeros([122,100])
    tick_marks = np.linspace(-1, 2800, 101)

    for sample in np.arange(122):
        left_tick = 0
        right_tick = 1
        for frame in np.arange(100):
            freq_matrix[sample,frame] = ((left_tick < dnase_matrix[:,sample]) & (dnase_matrix[:,sample] < right_tick)).sum()
            left_tick += 1
            right_tick += 1
    x = tick_marks[0:100]
    y = np.transpose(freq_matrix)

    plt.figure(figsize=(10, 10))
    plt.plot(x,y,'k',linewidth=0.1)
    plt.title('frequency plot of each column of dnase_matrix taking on various values across cell samples')
    plt.xlabel('count values in dnase_matrix')
    plt.ylabel('frequency of count value')
    plt.savefig(output_path+'frequency_plot_of_dnase_matrix_values.png')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='the path to the input files, ending with a forward slash', type=str)
    parser.add_argument('output_path', help='the desired path to the output files, ending with a forward slash', type=str)
    args = parser.parse_args()
    make_dnase_matrix_post_shell_script(args.input_path,args.output_path)
    print('Output path: ' + str(args.output_path))
    
main()

