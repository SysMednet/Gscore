import pandas as pd
import sys
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score
from scipy import stats

method=sys.argv[1]

path="/data1/LanYun/GScore/for_figs/fig5D/output/TCGA_347/TCGA_individual_analysis_16cancertype_metaz_for_volcano_OR_merge.txt"
df=pd.read_csv(path, sep='\t')
df_1=df.sort_values(by=['gene_survival_metaz'], ascending=False)
df_1['Gscore_association_metaz'] = df_1['Gscore_association_metaz'].astype(float)
df_1['NEA_association_metaz'] = df_1['NEA_association_metaz'].astype(float)

gene=df_1['gene'].tolist()
gene_survival_metaz_list=df_1['gene_survival_metaz'].tolist()
Gscore_association_metaz_list=df_1[f'{method}_association_metaz'].tolist()

# NEA_association_metaz_list=df_1['NEA_association_metaz'].tolist()

##Gscore
dic_gene_sur_metaz={}
dic_Gscore_prog_right={}
dic_Gscore_prog_left={}
dic_Gscore_prog_unpass={}
#right
Gscore_gene_survival_metaz_list_right=[]
Gscore_association_metaz_list_right=[]
#left
Gscore_gene_survival_metaz_list_left=[]
Gscore_association_metaz_list_left=[]
#unpass
Gscore_gene_survival_metaz_list_unpass=[]
Gscore_association_metaz_list_unpass=[]

for i,j,k in zip(gene,gene_survival_metaz_list,Gscore_association_metaz_list):
    if  k>=1.64 and (j>=1.96) : #k>1.64 and  and (j>=1.96) (j<=-1.96)   (j<1.96 and j>-1.96)
        #right
        dic_Gscore_prog_right.setdefault(i,[]).append(k)
        Gscore_gene_survival_metaz_list_right.append(j)
        Gscore_association_metaz_list_right.append(k)
    elif k>=1.64 and (j<=-1.96) :
        #left
        dic_Gscore_prog_left.setdefault(i,[]).append(k)
        Gscore_gene_survival_metaz_list_left.append(j)
        Gscore_association_metaz_list_left.append(k)
    else:
        #unpass
        dic_Gscore_prog_unpass.setdefault(i,[]).append(k)
        Gscore_gene_survival_metaz_list_unpass.append(j)
        Gscore_association_metaz_list_unpass.append(k)
    dic_gene_sur_metaz[i]=j
    
#right
Gscore_metaz_list_right=[dic_gene_sur_metaz[i] for i in sorted(dic_Gscore_prog_right.keys()) ] 
Gscore_qval_median_list_right=[np.median(dic_Gscore_prog_right[i]) for i in sorted(dic_Gscore_prog_right.keys())]
Gscore_qval_mean_list_right=[np.mean(dic_Gscore_prog_right[i]) for i in sorted(dic_Gscore_prog_right.keys())]
#left
Gscore_metaz_list_left=[dic_gene_sur_metaz[i] for i in sorted(dic_Gscore_prog_left.keys()) ] 
Gscore_qval_median_list_left=[np.median(dic_Gscore_prog_left[i]) for i in sorted(dic_Gscore_prog_left.keys())]
Gscore_qval_mean_list_left=[np.mean(dic_Gscore_prog_left[i]) for i in sorted(dic_Gscore_prog_left.keys())]


plt.figure(dpi=500)
plt.axhline(y=1.64, color='black', linestyle='--',lw=0.5)
plt.axvline(x=+1.96, color='black', linestyle='--',lw=0.5)
plt.axvline(x=-1.96, color='black', linestyle='--',lw=0.5)
plt.axvline(x=0, color='black', linestyle='--',lw=0.5)
plt.scatter(Gscore_gene_survival_metaz_list_left, Gscore_association_metaz_list_left,color='#00b050', marker='o',s=5,label="Gscore")
plt.scatter(Gscore_gene_survival_metaz_list_right, Gscore_association_metaz_list_right,color='r', marker='o',s=5,label="Gscore")
plt.scatter(Gscore_gene_survival_metaz_list_unpass, Gscore_association_metaz_list_unpass,color='gray', marker='o',s=5,label="Gscore")
x=[-10,-8,-6,-4,-2,0,2,4,6,8,10]
# y=[0,20,40,60,80,100,120,140,160]
y=[5,10,15,20,25,30,35,40]
plt.xticks(x)
plt.yticks(y)
plt.ylim(ymin=0)

plt.xlabel('Meta-z(prognostic significance)')
plt.ylabel('Meta-z(association qvalue)')
# plt.show()
path_save=f"/data1/LanYun/GScore/for_figs/fig5D/output/TCGA_347/5D_TCGA_OR_volcano_{method}.tiff"
plt.savefig(path_save)
