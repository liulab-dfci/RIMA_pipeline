#!/usr/bin/env python3
import os
import argparse
import logging
import pandas as pd
import numpy as np
import multiprocessing as mp

__doc__="""
-------------------------------------------------------- Microsatellite Instability Score ------------------------------------
This is the script to calculate the MSI score for your input samples based on expression profile.
User needs to provide a gene expression file and the cancer type. 
        Note:
        [Gene Expression file]
            1.  If it's possible, please input a normalized expression file follows the instruction:
                The gene expression value should be normalized toward a control sample which could be either normal tissues related 
                with a cancer type or mixture sample from diverse tumor samples. The log2(RPKM+1) values from a RNA-seq experiment 
                may not be meaningful unless a good reference control is available to adjust the batch effect and cancer type difference. 
                In our study, we used the all sample average in each study as the normalization control. 
                Otherwise, please enable --normalize flag to run the code. We'll do the normalization for you by:
                1. Do the log2(x+1) transformation
                2. Subtract the MEAN from your samples. 
                    - TCGA with normal samples: the MEAN is the average across normal samples
                    - Other datasets: the MEAN is the average across all samples.

            2.  If it's possible, please convert your gene identifier into Hugo Symbol based on your annotation GTF files. 
                Otherwise, we will use our annotation GTF to do the conversion, which is gencode v27.

        Please make sure the following files are under the script folder:
        .
        ├── main.py
        ├── msi_TCGA-STAD.lambda.min.coef 
        ├── gene.map
        ├── name.map
        
        Output example:
                MSI
Pt1     0.5372564094454566
Pt10    0.2577301118792513
Pt12    0.4296686626272194
"""

Scripts_Folder = os.path.dirname(os.path.realpath(__file__))

def main():
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o','--out_dir',type=str,required=True,
    help='Path to the folder you want to store the MSI score and file will be named the same as the input expression file')
    parser.add_argument('-e', '--exprsn_folder', type=str,
                        default=os.path.join(os.getenv('HOME'),'work/DATA/TCGA/hg38/RNASeq/fpkm/origin/'),
                        help='Path to the folder that only stored your expression profiles, or path of asingle expression file.')
    parser.add_argument('-m', '--msi_model', type=str,
                        default=os.path.join(
                            Scripts_Folder, 'msi_TCGA-STAD.lambda.min.coef'),
                        help='Path to the MSI model coefs.')
    parser.add_argument('--normalize', dest='normalize',action='store_true',
                        help='Whether to do normalization on input gene expression file. Default is False')
    parser.add_argument('--log2_transform', dest='log2_transform', action='store_true',
                        help='Whether to do log2(1+x) transform on input gene expression file. Default is False')
    parser.add_argument('--q_norm', dest='q_norm', action='store_true',
                        help='Whether to do quantile normalization on input gene expression file. Default is False')
    parser.add_argument('--log', type=str,
                        default='./stdout.log',
                        help='Path to store the logging information.')

    

    args = parser.parse_args()
    logging.basicConfig(filename=args.log,
                        filemode='w', level=logging.DEBUG)

    MSI_SIGNATURE=pd.read_csv(args.msi_model, index_col = 0,sep='\t')
    MSI_SIGNATURE.columns=['MSI']
    MSI_SIGNATURE['MSS']=- MSI_SIGNATURE['MSI']
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    logging.info(args)

    if os.path.isfile(args.exprsn_folder):
       getMSI('.', args.exprsn_folder, args.out_dir,
              MSI_SIGNATURE, args.normalize, args.log2_transform, args.q_norm)
    else:
        exprsn_list = [ f for f in os.listdir(args.exprsn_folder) if not f.startswith('.')]
        pool = mp.Pool(mp.cpu_count())
        pool.starmap(
            getMSI, [(args.exprsn_folder, f, args.out_dir, MSI_SIGNATURE, args.normalize, args.log2_transform, args.q_norm) for f in exprsn_list])
        pool.close()
    
    logging.info('All results are stored in {}'.format(args.out_dir))


def toEntrez(expression, BASE_DIR=Scripts_Folder):
    ''' Convert expression matrix with Ensemble ID or Gene symbol as index to expression matrix with Entrez ID
        GTF Gencode version 27
    Parameters
    ----------
    expression : pandas.DataFrame
        Expression matrix with gene identifier as index
    BASE_DIR: str
        Path to folders contains gene.map and name.map files
    Returns
    -------
    None / pandas.DataFrame
        Expression matrix with Entrez ID as index
    None / str
        Error information
    '''

    # translate the number of expression
    name_map = pd.read_csv(os.path.join(
        BASE_DIR, 'name.map'), index_col=0,sep='\t')['Entrez']
    ENS_TO_ENTREZ = pd.read_csv(os.path.join(BASE_DIR, 'gene.map'),sep='\t')[
        ['EntrezID', 'EnsembleID']].dropna(axis=0, how='any').set_index(['EnsembleID'])['EntrezID']
    ENS_TO_ENTREZ = ENS_TO_ENTREZ[~ ENS_TO_ENTREZ.index.duplicated()]

    cnt_nonint = sum(
        map(lambda v: not isinstance(v, np.number), expression.index))

    if cnt_nonint > 0:
        ind_g = name_map[expression.index.map(lambda x:x.upper())]
        ind_e = ENS_TO_ENTREZ[expression.index.map(lambda x:x.upper())]

        map_count = len(ind_g.dropna()) - len(ind_e.dropna())

        expression.index = ind_g if map_count > 0 else ind_e
        flag = (~ pd.isnull(expression.index))

        # too few genes
        if sum(flag) < expression.shape[0]/2 or sum(flag) < 10:
            expression = None
            errors = 'Only ' + str(sum(flag)) + ' out of ' + str(
                expression.shape[0]) + ' gene names are found. We only support: Ensembl ID, EntrezID, and Gene symbol'
        else:
            expression = expression.loc[flag]
            errors = None
    else:
        # Input matrix is index by Entrez ID already
        errors = None
    return expression, errors

def quantile_norm(df):
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df


def getTN(df, s):
    DIC = {'tumor': '0[0-9]$',
           'normal': '1[0-9]$'}
    return df.loc[:, df.columns.str.contains(DIC[s])]

def estMsi(exprs,msi_model):
    ol_gene = exprs.columns.intersection(msi_model.index)
    miss_gene = 1-(float(len(ol_gene)) / msi_model.shape[0])
    if miss_gene == 1:
        errors = 'No matched genes'
    elif miss_gene > 0.0:
        errors = '%.2f%% MSI signature genes are missing on input expression profile.' % (
                miss_gene*100)
    else:
        errors = None
    msi_score = np.exp(exprs[ol_gene].dot(msi_model.loc[ol_gene]))
    msi_score = msi_score.apply(lambda v: v/msi_score.sum(axis=1), axis=0)
    # in case there is duplicated samples
    msi_score = msi_score.groupby(level=0).mean()
    return msi_score['MSI'],errors


def getMSI(base_dir, filename, out_dir, msi_model, normalize, log2_transform, q_norm):
    obj = pd.read_csv(os.path.join(base_dir, filename),index_col=0,sep="\t")
    # Determining hwo to transfer the identifier
    # if the input expression profile with identifier:
    ## - Hugo symbol: Do NOT do any convertion
    ## - Entrez ID: Transfer MSI model to EtrezID indentifier
    ## - Ensemble ID: Transfer both of MSI model and input expression profile to EntrezID

    cnt_nonint = sum(
        map(lambda v: not isinstance(v, np.int), obj.index))
    if cnt_nonint >0:
        cnt_nonEsembel = sum(map(lambda v: not v.startswith('ENSG'), obj.index))
        if cnt_nonEsembel <= 0:
                logs = '-'*30 + \
                    '\n[Identifer Conversion] Input expression with ENSEMBLE ID as gene identifier.' + \
                    '\n[Identifer Conversion] Coverting indentifier of both input expression and MSI signature to EntrezID'
                obj,errors = toEntrez(obj)
                if not errors is None:
                    raise ValueError(errors)
                msi_model = toEntrez(obj)
        else:
                logs = '-'*30 + \
                    '\n[Identifer Conversion] Input expression with Hugo Symbol as gene identifier.' + \
                    '\n[Identifer Conversion] Stay the same identifier.'
                
    else: 
        logs = '-'*30 + \
            '\n[Identifer Conversion] Input expression with EntrezID as gene identifier.' + \
            '\n[Identifer Conversion] Coverting indentifier of MSI signature to EntrezID'
        msi_model, _ = toEntrez(msi_model)

    obj = obj.groupby(level=0).mean()
    if log2_transform or normalize:
        logs = logs + \
             '\n[RUNNING] Log2(x+1) transform on samples {}'.format(filename)
        obj = np.log(1+obj)
       
    if q_norm or normalize:
        logs = logs + \
             '\n[RUNNING] Quantile normalization on samples {}'.format(
                 filename)
        obj_expr = quantile_norm(obj)

    if normalize:
        # for TCGA data, normalize by normal samples
        if all(obj_expr.columns.str.contains('TCGA')):
            logs = logs+'\n[RUNNING] Detected all samples are from TCGA'
            obj_tumor = getTN(obj_expr, 'tumor')
            norm = obj_tumor.mean(axis=1)
            if getTN(obj_expr, 'normal').shape[1] > 1:
                logs = logs+'\n[RUNNING] Subtract the mean across normal samples'
                norm = getTN(obj_expr, 'normal').mean(axis=1)
            else:
                logs = logs+'\n[RUNNING] Subtract out the mean across all sample'
            obj_expr = obj_tumor.subtract(norm, axis=0)
            obj_expr.columns = obj_expr.columns.map(lambda x: x[:-3])
            obj_expr = obj_expr.T.groupby(level=0).mean()
        else:
            logs = logs+'\n[RUNNING] Detected all samples are NOT from TCGA'
            logs = logs+'\n[RUNNING] Subtract out the mean across all samples'
            obj_expr = obj_expr.subtract(obj_expr.mean(axis=1), axis=0).T
    else:
        obj_expr = obj.T
        logs = logs + '\n[RUNNING] Does not do normalization on input expression profile'

    msi,errors = estMsi(obj_expr, msi_model=msi_model)
    msi.name = 'MSI'
    if not errors is None:
        logs = logs + \
            '\n[WARNING]' + errors
    msi.to_frame().to_csv(os.path.join(out_dir,"msi_est_score.txt"),sep='\t')
    logs = logs+'\n[DONE] {}\n'.format(filename) + '-'*30

    logging.info(logs)


if __name__ == "__main__":
    main()