#!/usr/bin/env python3
import pandas as pd 
import numpy as np
import os
import argparse

__doc__="""
-------------------------------------------------------- Tumor Immune Dysfunction and Exclusion (TIDE) Score ------------------------------------
This is the script to calculate the TIDE score and immunesuppresive metrics for your gene expression profile.
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
                2. Subtract the average across your samples.

            2.  If it's possible, please convert your gene identifier into Entrez ID based on your annotation GTF files. 
                Otherwise, we will use our annotation GTF to do the conversion, which is gencode v27.
        [Cancer Type]
            We validated TIDE performance on predicting anti-PD1 and anti-CTLA4 response across several melanoma datasets and a limited dataset 
            of non-small cell lung cancer (NSCLC). TIDE may not work on cancer types other than melanoma and NSCLC (e.g., glioblastoma, or renal 
            cell carcinoma) and therapies other than anti-PD1 and anti-CTLA4 (e.g., anti-PDL1, or Car T). 

        Please make sure the following files are under the script folder:
        .
        ├── gene.map
        ├── main.py
        ├── name.map
        ├── signature
        ├── signature.std
        ├── signature.exclusion
        
        Output example:
                               Dysfunction     Exclusion       TIDE    CTL.flag        MDSC    CAF     TAM M2
TCGA-OR-A5J1-01 -0.45519311143117797    0.8290763802485979      -0.45519311143117797    True    0.03426937822731616     0.07496213114432222     0.006110807445246081
TCGA-OR-A5J2-01 -0.26669077391968443    0.8175876656653063      -0.26669077391968443    True    0.02503334313982659     0.08875374887296407     0.001955749219099661
TCGA-OR-A5J3-01 -0.17363257328368162    0.7055251453020146      0.7055251453020146      False   0.02972893945878613     0.07374222584556789     -0.003297292809209788
TCGA-OR-A5J5-01 -0.15379716859521772    0.8428418293971209      -0.15379716859521772    True    0.03746802552639951     0.0897389437620849      -0.007868256457221412
TCGA-OR-A5J6-01 -0.08426969004213913    0.300811121066252       -0.08426969004213913    True    -0.009514339999834164   0.04673702522676173     0.004314658757363608
TCGA-OR-A5J7-01 -0.23180055341324793    0.5320869807612111      -0.23180055341324793    True    0.022960453836467765    0.043525796777474285    0.005957144298245646
----------------------------------------------------------------------------------------------------------------------------------------------
"""
def toEntrez(expression,BASE_DIR):
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
    name_map = pd.read_table(os.path.join(
        BASE_DIR,'name.map'), index_col=0)['Entrez']
    ENS_TO_ENTREZ = pd.read_table(os.path.join(BASE_DIR,'gene.map'))[
        ['EntrezID', 'EnsembleID']].dropna(axis=0, how='any').set_index(['EnsembleID'])['EntrezID']
    ENS_TO_ENTREZ = ENS_TO_ENTREZ[~ ENS_TO_ENTREZ.index.duplicated()]

    cnt_nonint = sum(
        map(lambda v: not isinstance(v, np.int), expression.index))

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

def get_expression_flag(expression):
    """ Distinguish hot tumor from cold tumor"""
    set_CTL = [925, 926, 3001, 3002, 5551]
    CTL = expression.loc[set_CTL].min() > 0
    CTL.name = 'CTL.flag'
    return CTL


def estTIDEscore(expression, cancer, BASE_DIR, normalize=False):
    signature_file = '/'.join([BASE_DIR,'signature'])
    signature = pd.read_table(signature_file, index_col=0)
    signature_sd = pd.read_table(signature_file + '.std', index_col=0)
    signature_exclusion = pd.read_table(
        signature_file + '.exclusion', index_col=0)
    # Only consider studied cancer types
    if cancer == 'Melanoma':
        signature_sd = signature_sd.loc['SKCM.RNASeq']
    elif cancer == 'NSCLC':
        signature_sd = (signature_sd.loc['LUSC.RNASeq.norm_subtract'] +
                        signature_sd.loc['LUAD.RNASeq.norm_subtract'])/2
    elif cancer == 'Other':
        # use the melanoma rule for approximation
        signature_sd = signature_sd.loc['SKCM.RNASeq']
    else:
        return None, 'Cannot recognize ' + cancer

    # merging ratio for different signatures
    expression = pd.read_table(expression, index_col=0)

    # translate the number of expression
    expression, errors = toEntrez(expression, BASE_DIR)

    if errors != None:
        return None, errors
   
    expression = expression.groupby(level=0).mean()
    expression = expression.dropna(axis=0)
    # do normalization
    if normalize:
        expression = np.log2(expression + 1)
        expression = expression.apply(lambda v: v-v.mean(),axis=1)

    flag = get_expression_flag(expression)
    correlation = signature.apply(lambda v: expression.corrwith(v))
    correlation = correlation.divide(signature_sd)

    correlation['TIDE'] = correlation['Exclusion']
    correlation.loc[flag, 'TIDE'] = correlation.loc[flag, 'Dysfunction']
    # add in T-cell exclusion signatures
    correlation_exclusion = signature_exclusion.apply(
        lambda v: expression.corrwith(v))

    result = pd.concat(
        [correlation, flag,correlation_exclusion], axis=1)
    return result, None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-e', '--expression', required=True,type=str,
                        help="Path to gene expression file")
    parser.add_argument('-o','--output',required=True,type=str,
                        help="Path to store TIDE score")
    parser.add_argument('-c', '--cancer', required=True, type=str,
                        help="The cancer type of your cohorts. Options: 'Melanoma','NSCLC','Other'") 
    parser.add_argument('--base_path',required=False,type=str,default=None,
                        help='Path of this script location, where has gene.map, name.map, and singature* files. Default is None')
    parser.add_argument('--normalize', dest='normalize',action='store_true',
                        help='Whether to do normalization on input gene expression file. Default is False')
    
    scripts_folder = os.path.dirname(os.path.realpath(__file__))
    args = parser.parse_args()
    print(args)
    if args.base_path is None:
        args.base_path = scripts_folder 
    result,error = estTIDEscore(expression=args.expression,cancer=args.cancer,BASE_DIR=args.base_path,normalize=args.normalize)
    if error is None:
        out_dir = os.path.dirname(args.output)
        print(out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        result.to_csv(args.output,sep='\t')
        print('[Success] The TIDE score has been store on {}'.format(os.path.abspath(args.output)))
    else:
        print(error)
