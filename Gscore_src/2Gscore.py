import numpy as np
from scipy.stats import hypergeom
import statsmodels.stats
import statsmodels.stats.multitest
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-n", type=str , default='./sample_data_output/example_pcc.txt' , help='The Pearson correlation network.')
parser.add_argument("-g", type=str , default='./sample_data_input/example_allDEG.txt' , help='The DEG list file, the DEG IDs are seperated by "\n".')
parser.add_argument("-s", type=str , default='./sample_data_input/example_GeneSet_KEGG_v102.txt' , help='The gene set collection file.')
parser.add_argument("-q", type=str , default='./sample_data_input/example_query_DEG.txt' , help='The query DEG file, the DEG IDs are seperated by "\n".')
parser.add_argument("-t", type=float , default=0.7 , help='The threshold for Pearson correlation coefficient.')
parser.add_argument("-io", type=str , default='./sample_data_output/example_individual_DEG_result.txt' , help='The output file name for individual DEG result.')
parser.add_argument("-Lo", type=str , default='./sample_data_output/example_DEG_list_result.txt' , help='The output file name for DEG list result.')

args = parser.parse_args()
input_network,input_DEG_list,input_geneset,input_query,input_cutoff, output_path_single, output_path_com = args.n, args.g, args.s, args.q, args.t, args.io, args.Lo

pcc_cutoff = input_cutoff

#read DEG list
thresh={}
fp=open(input_DEG_list,'r')
DEGs=fp.read().strip().split('\n')
fp.close()
DEGs = set(DEGs)
for DEG in DEGs:
    thresh[DEG+'-'+DEG]=1

#read pcc data
fp=open(input_network,'r') #PCC
lines=fp.readlines()
fp.close()
for line in lines[1:]: 
    temp=line.strip().split('\t') 
    if (np.abs(float(temp[2]))>=pcc_cutoff):
        thresh[temp[0]+'-'+temp[1]]=1
        thresh[temp[1]+'-'+temp[0]]=1
    else:
        thresh[temp[0]+'-'+temp[1]]=0
        thresh[temp[1]+'-'+temp[0]]=0

#read gene set data
fp=open(input_geneset,'r')
Pathway={}
AllPathway={}
lines=fp.readlines()
fp.close()
for line in lines:
    temp=line.strip().split('\t')
    pathway_gene=temp[2].split(',')
    pathway_gene_DEG=[]
    for j in pathway_gene:
        if j in DEGs:
            pathway_gene_DEG.append(j)
            AllPathway[j]='1'
    if len(pathway_gene_DEG) != 0:
        Pathway[temp[0]]=pathway_gene_DEG

#%% calculate individual DEG information
single_dic = {}
DEG_list_input=[]

for j in DEGs:
    (N,M)=(len(AllPathway),0)
    for k in AllPathway.keys():
        if j+'-'+k in thresh:
            M += thresh[j+'-'+k]
    for l in Pathway.keys():
        (n,m)=(len(Pathway[l]),0)
        for p in Pathway[l]:
            #get involved gene
            if j+'-'+p in thresh:
                m += thresh[j+'-'+p]

        P_value=(hypergeom.sf(m-1,N,M,n)) #p-val of hyper
        
        if ~(M==0 or N==0 or n==0) and (j not in single_dic):
            single_dic[j] = {'GeneSet_ID':[],'gene_ID':[],'N':[],'M':[],'n':[],'m':[],'p_value':[],'FDR_q_value':[]}
        if m!=0:
            single_dic[j]['gene_ID'].append(j)
            single_dic[j]['GeneSet_ID'].append(l)
            single_dic[j]['p_value'].append(P_value)
            single_dic[j]['N'].append(N)
            single_dic[j]['M'].append(M)
            single_dic[j]['n'].append(n)
            single_dic[j]['m'].append(m)

        DEG_list_input.append(l+'\t'+j+'\t'+str(N)+'\t'+str(M)+'\t'+str(n)+'\t'+str(m)+'\t'+str(P_value)) #write p j

#individual fdrq correction
for g in single_dic:
    fdrQ = list(statsmodels.stats.multitest.fdrcorrection(single_dic[g]['p_value'], alpha=0.05, method='indep', is_sorted=False)[1])
    single_dic[g]['FDR_q_value'] = fdrQ

#output individual result
output_dic = {'GeneSet_ID':[],'gene_ID':[],'N':[],'M':[],'n':[],'m':[],'p_value':[],'FDR_q_value':[]}
for g in single_dic:
    for index in single_dic[g]:
        output_dic[index]+=single_dic[g][index]
del single_dic

op=open(output_path_single,'w')
output_df = pd.DataFrame(output_dic)
output_df.to_csv(output_path_single,sep='\t', index=False)
op.close()

#%% calculate DEG list information
Pathway_nm={}
Pathway={}
for line in DEG_list_input:
    temp=line.split('\t')
    Pathway[temp[0]]='1'
    Pathway_nm[temp[0]+'|'+temp[1]]=[temp[2],temp[3],temp[4],temp[5]]

#input interest DEG. Here, we input top 100 DEG as example
fp=open(input_query,'r')
interest_DEG_list=fp.read().strip().split('\n')
fp.close()

DEG_list_output_dic = {'GeneSet_ID':[],'N':[],'M':[],'n':[],'m':[],'p_value':[],'FDR_q_value':[]}
for l in Pathway.keys():
    (N,M)=(0,0)
    (n,m)=(0,0)
    left = []
    right = []
    for gene in interest_DEG_list:
        if (l+'|'+gene) in Pathway_nm.keys():
            N+=int(Pathway_nm[l+'|'+gene][0])
            M+=int(Pathway_nm[l+'|'+gene][1])
            n+=int(Pathway_nm[l+'|'+gene][2])
            m+=int(Pathway_nm[l+'|'+gene][3])
            
    if (N!=0):
        DEG_list_output_dic['GeneSet_ID'].append(l)
        DEG_list_output_dic['p_value'].append(hypergeom.sf(m-1,N,M,n))
        DEG_list_output_dic['N'].append(N)
        DEG_list_output_dic['M'].append(M)
        DEG_list_output_dic['n'].append(n)
        DEG_list_output_dic['m'].append(m)

fdrQ = list(statsmodels.stats.multitest.fdrcorrection(DEG_list_output_dic['p_value'], alpha=0.05, method='indep', is_sorted=False)[1])
DEG_list_output_dic['FDR_q_value'] = fdrQ

op=open(output_path_com,'w')
DEG_list_output_df = pd.DataFrame(DEG_list_output_dic)
DEG_list_output_df.to_csv(output_path_com,sep='\t', index=False)
op.close()
