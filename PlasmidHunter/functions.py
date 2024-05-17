import shutil
import sys
import os
import subprocess as subp
from multiprocessing import Pool
import pandas as pd
import numpy as np
from glob import glob
import json
import pickle
from collections import Counter
from Bio import SeqIO
import requests
import gzip
from sklearn.naive_bayes import GaussianNB

def unzip(infile):
    with gzip.open(infile, 'rb') as f_in:
        with open(infile[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            os.remove(infile)

def download_file(url, outfile):
    response = requests.get(url)
    with open(outfile, 'wb') as file:
        file.write(response.content)

def call(cmd='',out='',err=''):
    if subp.call('set -ex; '+cmd,shell=True)==0:
        print('\n'+out+'\n')
    else:
        print('\n'+err+'\n')
        sys.exit(1)

def filter_fasta(infile, outfile, length=1000):
    outfile = open(outfile, 'w')
    n = 0
    for rec in SeqIO.parse(infile, 'fasta'):
        if len(rec.seq) >= length:
            SeqIO.write(rec, outfile, 'fasta')
            n += 1
    outfile.close()

    return(n)

def split_fasta(infile, outdir, n=10):
    print('Now splitting input fasta file ...')
    dict1 = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
    ids = dict1.keys()
    split = np.array_split(list(ids), n)
    for i in range(len(split)):
        outfile = open('%s/%d.fasta' % (outdir, i), 'w')
        for j in split[i]:
            SeqIO.write(dict1[j], outfile, 'fasta')
        outfile.close()

def prodigal(indir, f):
    call('prodigal -a %s/%s.faa -o %s/%s.out -i %s/%s -p meta -q' % (indir, f, indir, f, indir, f), 'Prodigal prediction for %s completed.' % f, 'Prodigal prediction for %s failed.' % f)
 
def prodigal_predict(indir, prodigal_func):
    print('Now predicting protein sequences using Prodigal ...')
   
    files = [i for i in os.listdir(indir) if i.split('.')[-1] == 'fasta']
    pool = Pool(len(files))
    for f in files:
        pool.apply_async(prodigal_func, (indir, f,))
    pool.close()
    pool.join()

def diamond_alignment(query, database, cpu=10): 
    print('Now aligning protein sequences against database using Diamond ...')
    call('diamond blastp --query '+query+' --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 --id 30 --query-cover 70 --db '+database+' --out '+query+'.daa --outfmt 6 --threads '+str(cpu), 'Diamond alignment completed.', 'Diamond alignment failed.')

def daa_to_hits(daa):
    print('Now parsing Diamond alignment results ...')
    dict1 = {}
    df = pd.read_csv(daa, sep = '\t', header = None, index_col = 0)
    df = df.loc[:,1]
    df.index = ['_'.join(i.split('_')[:-1]) for i in df.index]
    dict1 = df.groupby(level=0).apply(list).to_dict()

    return(dict1)

def list_to_df(li1, name):
    dict1 = Counter(li1)
    df = pd.DataFrame({name: dict1})
    df[df>0] = 1
    
    return(df)

def mergefiles(indir, outfile):
    out = []
    for i in glob(f'{indir}/*.faa'):
        out.append(open(i, 'r').read())
    open(outfile, 'w').write(''.join(out))

def gene_content_profile(dna_file, prodigal_func, database, cpu=10): # get gene content profile of sequences using a diamond output
    print('Now generating gene content profile ...')
    temp_dir = dna_file+'_temp'
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
    split_fasta(dna_file, temp_dir, n=cpu)
    prodigal_predict(temp_dir, prodigal_func)
    mergefiles(temp_dir, f'{temp_dir}/proteins.faa')
    diamond_alignment(temp_dir+'/proteins.faa', database, cpu=cpu)
    dict1 = daa_to_hits(temp_dir+'/proteins.faa.daa')
    shutil.rmtree(temp_dir)
    
    return(dict1)

def load_model_params(param_file): #
    with gzip.open(param_file, 'rt', encoding='utf-8') as file:
        json_str = file.read()
    params = json.loads(json_str)

    model = GaussianNB()
    model.class_prior_ = np.array(params['class_prior_'])
    model.theta_ = np.array(params['theta_'])
    model.var_ = np.array(params['var_'])
    model.classes_ = np.array(params['classes_'])

    return model

def predict(dict1, feature_genes, clf):
    print('Predicting using a model ...')
    dict1 = list(dict1.items())
    l = [dict1[i:i+50] for i in list(range(0, len(dict1), 50))]
    l = [dict(i) for i in l]
    out = pd.DataFrame()

    for d in l:
        df = pd.concat([list_to_df(d[i], i) for i in d], axis=1)
        df = df.fillna(0)
        df2 = pd.DataFrame(index=feature_genes)
        df2 = pd.concat([df2, df], axis=1)
        df2 = df2.transpose()
        df2 = df2.loc[:,feature_genes]
        df2 = df2.fillna(0)
        pred_proba = clf.predict_proba(np.array(df2))
        pred = np.array([np.argmax(i) for i in pred_proba])
        df3 = pd.DataFrame(np.hstack((pred.reshape((-1,1)), pred_proba)), index=df2.index, columns=['Prediction (0: chromosome, 1: plasmid)', 'Probability of 0', 'Probability of 1'])
        out = pd.concat([out, df3])

    return(out)

