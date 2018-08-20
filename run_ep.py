
from __future__ import print_function
import argparse

def main():
    
    parser = argparse.ArgumentParser(description='')
    
    outdir_parser=argparse.ArgumentParser(add_help=False)
    outdir_parser.add_argument('--outdir')

    enh_mat_parser=argparse.ArgumentParser(add_help=False)
    enh_mat_parser.add_argument('--enh_mat')

    gene_mat_parser=argparse.ArgumentParser(add_help=False)
    gene_mat_parser.add_argument('--gene_mat')

    subparsers = parser.add_subparsers(help='', dest='command')
    subparsers.required = True #http://bugs.python.org/issue9253#msg186387
    
    preprocess_parser=subparsers.add_parser('preprocess',parents=[outdir_parser,enh_mat_parser,gene_mat_parser],help='')
    ep_parser=subparsers.add_parser('ep',parents=[outdir_parser],help='')
    
    args = vars(parser.parse_args())
    command = args.pop("command", None)

    command_methods = {'preprocess': preprocess,
                         'ep': ep}
    command_methods[command](**args)

def preprocess(outdir,enh_mat,gene_mat):
    print('build a list of which genes correspond to which ')

def ep(outdir):
    print('ep')

main()



