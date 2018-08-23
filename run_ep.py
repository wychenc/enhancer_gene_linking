
from __future__ import print_function
import argparse
import subprocess as subp
import gzip
import os
import fnmatch
import re
from scipy.stats.stats import pearsonr
import numpy as np
import base64

def main():
    
    parser = argparse.ArgumentParser(description='')
    
    outdir_parser=argparse.ArgumentParser(add_help=False)
    outdir_parser.add_argument('--outdir',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/outdir')

    enh_mat_parser=argparse.ArgumentParser(add_help=False)
    enh_mat_parser.add_argument('--enh_mat',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/enh_mat.gz')

    gene_mat_parser=argparse.ArgumentParser(add_help=False)
    gene_mat_parser.add_argument('--gene_mat',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/gene_mat.gz')

    enh_pos_parser=argparse.ArgumentParser(add_help=False)
    enh_pos_parser.add_argument('--enh_pos',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/enh_pos.bed')

    gene_pos_parser=argparse.ArgumentParser(add_help=False)
    gene_pos_parser.add_argument('--gene_pos',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/gene_pos.bed')

    dist_parser=argparse.ArgumentParser(add_help=False)
    dist_parser.add_argument('--dist',type=int,default=1000000)

    bashrc_parser=argparse.ArgumentParser(add_help=False)
    bashrc_parser.add_argument('--bashrc',default='/mnt/lab_data/kundaje/users/oursu/code/forcici/enhancer_gene_linking/ep_bashrc')

    methods_parser=argparse.ArgumentParser(add_help=False)
    methods_parser.add_argument('--methods',default='correlation')
    
    subparsers = parser.add_subparsers(help='', dest='command')
    subparsers.required = True #http://bugs.python.org/issue9253#msg186387
    
    preprocess_parser=subparsers.add_parser('preprocess',parents=[outdir_parser,enh_mat_parser,gene_mat_parser,enh_pos_parser,gene_pos_parser,dist_parser,bashrc_parser],help='')
    ep_parser=subparsers.add_parser('ep',parents=[outdir_parser,enh_mat_parser,gene_mat_parser,methods_parser],help='')
    
    args = vars(parser.parse_args())
    command = args.pop("command", None)

    command_methods = {'preprocess': preprocess,
                         'ep': ep}
    command_methods[command](**args)

def preprocess(outdir,enh_mat,gene_mat,enh_pos,gene_pos,dist,bashrc):
    print('------------------ checking input files')
    enh_cells=gzip.open(enh_mat,'r').readline().decode('utf8').strip().split('\t')[4:]
    gene_cells=gzip.open(gene_mat,'r').readline().decode('utf8').strip().split('\t')[4:]
    enh_pos=gzip.open(enh_mat,'r').readline().decode('utf8').strip().split('\t')[:5]
    gene_pos=gzip.open(gene_mat,'r').readline().decode('utf8').strip().split('\t')[:5]

    if len(set(enh_cells))!=len(enh_cells):
        print('Repeated column in the enhancer matrix')
        exit
    if len(set(gene_cells))!=len(gene_cells):
        print('Repeated column in the gene matrix')
        exit
    common_cells=set(enh_cells).intersection(gene_cells)
    if len(common_cells)<len(set(enh_cells)):
        print('WARNING: cell lines in enhancer matrix not found in gene matrix: '+','.join(list(set(enh_cells).setdiff(gene_cells))))
    if len(common_cells)<len(set(gene_cells)):
        print('WARNING: cell lines in gene matrix not found in gene matrix: '+','.join(list(set(gene_cells).setdiff(enh_cells))))

    #TODO: check all elements in enh and gene matrix are in the --enh_pos and --gene_pos files
    #enh_pos and gene_pos matrices are extracted from enh_mat and gene_mat    

    print('------------------ assigning enhancer-gene pairs to investigate using a distance threshold')
    #write a script to run bedtools window. the script will first load the bashrc, then run the command
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/scripts'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/data'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/results'])
    sname=outdir+'/scripts/connect_enh2gene_pos.sh'
    s=open(sname,'w')
    s.write('. '+bashrc+'\n')
    s.write('bedtools window -w '+str(dist)+' -a '+enh_pos+' -b '+gene_pos+' | gzip > '+outdir+'/data/enh2gene.pos.gz'+'\n')
    s.close()
    run_script(sname)
    print('Result: '+outdir+'/data/enh2gene.pos.gz')
    
    print('------------------ splitting files (enh2gene, enh_matrix, gene_matrix) by chromosome')
    sname=outdir+'/scripts/split_by_chromo.sh'
    s=open(sname,'w')
    s.write('. '+bashrc+'\n')
    for f in [outdir+'/data/enh2gene.pos.gz']:
        s.write('zcat -f '+f+' | awk -v outfile='+outdir+'/data/'+os.path.basename(f)+' \'{print $0>outfile"."$1}\''+'\n')
    s.close()
    run_script(sname)
    
def run_script(sname):
    subp.check_output(['bash','-c','chmod 755 '+sname])
    subp.check_output(['bash','-c',sname])
    
def ep(outdir,enh_mat,gene_mat,methods):
    methods_list=methods.split(',')
    #go through each chromosome
    #extract the enh2gene pairs to compute
    #compute (yay!)

    for chromo_f in fnmatch.filter(os.listdir(outdir+'/data/'), 'enh2gene.pos.gz.*'):
        chromo=re.sub('enh2gene.pos.gz.','',chromo_f)
        print(chromo)

        files_dict={}
        for method in methods_list:
            files_dict[method]=gzip.open(outdir+'/results/'+method+'.'+chromo+'.gz','w')
        
        enh_data={}
        gene_data={}
        for line in gzip.open(enh_mat,'r').readlines()[1:]:
            items=line.decode('utf8').strip().split('\t')
            if items[0]!=chromo:
                continue
            enh=items[0]+':'+items[1]+'-'+items[2]
            values=np.array([float(x) for x in items[4:]])
            #TODO: make sure cell types are in the correct order
            if enh in enh_data:
                print('Enhancer '+enh+' repeated')
                exit
            enh_data[enh]=values

        for line in gzip.open(gene_mat,'r').readlines()[1:]:
            items=line.decode('utf8').strip().split('\t')
            if items[0]!=chromo:
                continue
            gene=items[0]+':'+items[1]+'-'+items[2]
            values=np.array([float(x) for x in items[4:]])
            #TODO: make sure cell types are in the correct order
            if gene in gene_data:
                print('Gene '+gene+' repeated')
                exit
            gene_data[gene]=values

        #TODO: index by gene names rather than the position in the genome
        
        #extract enh2gene pairs
        for line in open(outdir+'/data/enh2gene.pos.gz.'+chromo,'r').readlines():
            items=line.strip().split('\t')
            e_chr,e_start,e_end,g_chr,g_start,g_end=items[0],items[1],items[2],items[3],items[4],items[5]
            ename=e_chr+':'+e_start+'-'+e_end
            gname=g_chr+':'+g_start+'-'+g_end
            if ename in enh_data.keys() and gname in gene_data.keys():
                e_values=enh_data[ename]
                g_values=gene_data[gname]

                if 'correlation' in methods_list:
                    #TODO: deal with nans
		    if np.var(e_values) is not 0:
		    	enh2gene_value=get_corr(e_values,g_values)
                    else:	
		    	enh2gene_value='nan'

                    #write to file
                    files_dict['correlation'].write(str(e_chr+'\t'+e_start+'\t'+e_end+'\t'+g_chr+'\t'+g_start+'\t'+g_end+'\t'+str(enh2gene_value)+'\n').encode())

                #add your own methods

def get_corr(a,b):
    return pearsonr(a,b)[0]
                

main()



