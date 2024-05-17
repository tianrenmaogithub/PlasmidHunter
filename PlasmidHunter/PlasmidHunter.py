import argparse
import os
import pkg_resources
import shutil
import random
import json
import pickle
from .functions import *

parser=argparse.ArgumentParser(description='This program predicts the location (chromosome or plasmid) of contig seqeuences from an input fasta file based on gene content profile using machine learning. It outputs the prediction and probability of each contig.')
parser.add_argument('-i','--infile',dest='inp',required=True,help='Input fasta file. Only contigs > 1Kbp will be includded in the prediction.')
parser.add_argument('-c','--cpu',dest='cpu',required=False,default=10, type=int, help='Number of processes to launch for the Prodigal gene prediction and Diamond alignment. [Default: 10]')
parser.add_argument('-o','--outdir',dest='out',required=False,default='plasmid_prediction',help='Output directory of the prediction. The result file, predictions.tsv, includes the prediction results and corresponding probabilities [Default: plasmid_prediction]')
args=parser.parse_args()

def main():
    infile = args.inp
    model = pkg_resources.resource_filename('PlasmidHunter', 'model/guassiannb_model_params.json.gz')
    feature_genes = pkg_resources.resource_filename('PlasmidHunter', 'model/feature_genes.gz')
    database = f'{os.path.dirname(os.path.dirname(model))}/database/database.dmnd'
    if not os.path.exists(database):
        print(f'{database}\nDownloading database. This will run only once. Please wait some minutes ...')
        os.makedirs(os.path.dirname(database), exist_ok=True)
        download_file('https://zenodo.org/records/10431696/files/database.dmnd.gz?download=1', f'{database}.gz')
        unzip(f'{database}.gz')

    cpu = args.cpu
    outdir = args.out

    # filter short reads
    print('Filtering contigs < 1 Kbp out ...')
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    seq_num = filter_fasta(infile, outdir+'/input.fasta', 1000)
    cpu = min(seq_num, cpu)

    # generate gene profile
    dict1 = gene_content_profile(outdir+'/input.fasta', prodigal, database, cpu=cpu)

    # predicting
    #feature_genes = pickle.load(open(feature_genes, 'rb'))
    feature_genes = json.loads(gzip.open(feature_genes, "rt", encoding="utf-8").read())
    #clf = pickle.loads(pickle.load(open(model, 'rb')))
    clf = load_model_params(model)
    pred = predict(dict1, feature_genes, clf)
    pred.to_csv(outdir+'/predictions.tsv', sep='\t')
    os.remove(outdir+'/input.fasta')
    print('Prediction completed. See '+outdir+'/predictions.tsv for the prediction outputs and probabilities.\nPlease cite us if you find it helpful. Thank you!')

if __name__ == '__main__':
    main()

