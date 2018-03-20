
# coding: utf-8

# In[4]:

def make_dnase_matrix_pre_shell_script(output_path):
    # make 5 useful dnase file columns into lists
    import numpy as np
    lab_list=[]
    biosample_matched_rna_list=[]
    eid_matched_rna_list=[]
    pooled_tagAlign_list=[]
    idr_file_list=[]
    cell_list=[]
    color_list=[]
    dnase_folder_name_list=[]

    with open('/mnt/data/integrative/metadata/dnase_metadata_2016-12-05.txt','r') as f_txt: 
        for line in f_txt:
            if 'idr_file' in line: # skip over header line
                continue
            items = line.split('\t')
            lab_list.append(items[23])
            color_list.append(items[22])
            biosample_matched_rna_list.append(items[17])
            dnase_folder_name_list.append(items[0])
            eid_matched_rna_list.append(items[18])
            pooled_tagAlign_list.append(items[6])
            idr_file_list.append(items[9])
            cell_list.append(items[20])

    biosample_matched_rna = np.asarray(biosample_matched_rna_list)
    eid_matched_rna = np.asarray(eid_matched_rna_list)
    pooled_tagAlign = np.asarray(pooled_tagAlign_list)
    idr_file = np.asarray(idr_file_list)
    color = np.asarray(color_list)
    dnase_folder_name = np.asarray(dnase_folder_name_list)

    bios_eid_tag_idr_color_folder_all = np.stack((biosample_matched_rna,eid_matched_rna,pooled_tagAlign,idr_file,color,dnase_folder_name))
    bios_eid_tag_idr_color_folder_all = bios_eid_tag_idr_color_folder_all.transpose()
    f_txt.close()

    # restrict lab to "UW_Stam" and biosample_matched_rna to not blank, obtain reduced matrix
    # from prior computations: lab: 351, bios: 123, idx: 122 (to keep)
    pca_labels_final = []
    idx=[] # get indices that fit above restrictions
    idx_counter = 0
    for i in biosample_matched_rna_list:
        idx_counter += 1
        if i is not "":
            thelab = lab_list[idx_counter] 
            if "Stam" in thelab:
                idx.append(idx_counter)

    bios_eid_tag_idr_color_folder=np.asarray([0,0,0,0,0,0]) # add first row of zeros
    for i in idx:
        row = bios_eid_tag_idr_color_folder_all[i,:]
        bios_eid_tag_idr_color_folder = np.vstack((bios_eid_tag_idr_color_folder, row))
        pca_labels_final.append(cell_list[i])

    bios_eid_tag_idr_color_folder = np.delete(bios_eid_tag_idr_color_folder, (0), axis=0) # delete first row of zeros

    idr_files = bios_eid_tag_idr_color_folder[:,3] # make bed file of the idr_files to manipulation in shell script
    tag_files = bios_eid_tag_idr_color_folder[:,2]
    bios_eid_list = bios_eid_tag_idr_color_folder[:,0:2]
    pca_color = bios_eid_tag_idr_color_folder[:,4]
    dnase_folder_name = bios_eid_tag_idr_color_folder[:,5]

    np.savetxt(output_path + 'idr_files.bed', idr_files, fmt="%s", delimiter='\t')
    np.savetxt(output_path + 'tag_files.bed', tag_files, fmt="%s", delimiter='\t') # this is the column labels
    np.savetxt(output_path + 'bios_eid_for_final_dnase_matrix.txt', bios_eid_list, fmt="%s", delimiter='\t')
    np.savetxt(output_path + 'pca_labels.txt', pca_labels_final, fmt="%s", delimiter='\t')
    np.savetxt(output_path + 'pca_color_codes.txt', list(pca_color), fmt="%s", delimiter='\t')
    np.savetxt(output_path + 'folder_name.txt', list(dnase_folder_name), fmt="%s", delimiter='\t')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_path', help='the desired path to the output files, ending with a forward slash', type=str)
    args = parser.parse_args()
    make_dnase_matrix_pre_shell_script(args.output_path)
    print('Output path: ' + str(args.output_path) + 'idr_files.bed')
    
main()

