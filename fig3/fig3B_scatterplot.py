import os
import scipy as sp
import scipy.stats as st
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf

output_dir = 'fig3B_output_figure/'
methods = ['NEA','ROntoTools','SPIA','PADOG','GSA','GSEA','ORA']

file_path = open('D:/experiment/gscore/KEGG_api_347_hsa_pathway_20220518_v102.txt','r')
lines_path = file_path.readlines()
file_path.close()
ID_to_name = {}
for path in lines_path:
    tmp = path.strip().split('\t')
    ID_to_name[tmp[0]] = tmp[1]

data_df = pd.read_csv('fig3B_input.txt', delimiter = "\t")
sns.set(font_scale=3)

'''file_top = open(output_dir+'meta-z_sum_all.txt','w')
top_dic = {'pathway_ID':[],'pathway_name':[],'method1':[],'methods1_metaz':[],'method2':[],'methods2_metaz':[],'sum_of_metaz':[]}
top_df = pd.DataFrame(top_dic)'''

for m in methods:
    sns.set_style("ticks")
    ax2 = plt.figure(dpi = 300, figsize=(15,7))
    
    enrich_value1=abs(st.norm.ppf(0.05))
    if m=='GSEA':
        enrich_value2=abs(st.norm.ppf(0.05)) #0.25
    elif m=='GSA':
        enrich_value2=abs(st.norm.ppf(0.05)) #0.1
    else:
        enrich_value2=abs(st.norm.ppf(0.05))
    
    def annotate(data, **kws):
        tau, p = st.kendalltau(data[m],data['Gscore_fdrq'])
        
        ax = plt.gca()
        ax.text(0.5, 1.05, '\u03C4={:.2f}, p={:.2g}'.format(tau, p),
                transform=ax.transAxes,size=45,fontweight ='bold')
    ###top sum####
    top_tmp_dic = {'pathway_ID':[],'pathway_name':[],'method1':[],'methods1_metaz':[],'method2':[],'methods2_metaz':[],'sum_of_metaz':[]}
    top_tmp_dic['pathway_ID'] = list(data_df['Unnamed: 0'])
    top_tmp_dic['pathway_name'] = [ID_to_name[x] for x in list(data_df['Unnamed: 0'])]
    top_tmp_dic['method1'] = [m]*len(list(data_df['Unnamed: 0']))
    top_tmp_dic['methods1_metaz'] = list(data_df[m])
    top_tmp_dic['method2'] = ['Gscore_fdrq']*len(list(data_df['Unnamed: 0']))
    top_tmp_dic['methods2_metaz'] = list(data_df['Gscore_fdrq'])
    top_tmp_dic['sum_of_metaz'] = list(data_df[m]+data_df['Gscore_fdrq'])
    top_tmp_df = pd.DataFrame(top_tmp_dic)
    top_tmp_df = top_tmp_df.sort_values(by=['sum_of_metaz'], ascending=False)
    #top_df = top_df.append(top_tmp_df,ignore_index=True) #.head(3)
    
    tau, p = st.kendalltau(data_df[m],data_df['Gscore_fdrq'])
    if tau>0:
        color='red'
    else:
        color='green'
    
    ax2 = sns.lmplot(x=m, y='Gscore_fdrq', data=data_df, scatter_kws={"s": 100,'color':'black'}, line_kws={'color':color,'ls':'--','lw':4}, height=10, aspect=1.1,ci=0)
    ax2.map_dataframe(annotate)
    plt.axhline(enrich_value1, ls='--')
    plt.axvline(enrich_value2, ls='--')
    ax2.set(ylim=(-(np.max(data_df['Gscore_fdrq'])/100*1.5), None),xlim=(-(np.max(data_df[m])/100*1.5), None))
    
    plt.ylabel('Meta-z (Gscore)', fontweight='bold', size=50)
    plt.xlabel('Meta-z'+' ('+m+')', fontweight='bold', size=50)

    '''if np.max(data_df[m])>enrich_value: 
        plt.fill_between([enrich_value, np.max(data_df[m])+20], -10, 1.64, alpha=0.3, color='green')
    plt.fill_between([-10, enrich_value], 1.64, np.max(data_df[g]+20), alpha=0.3, color='blue')'''
    
    plt.savefig(output_dir+'fig3B_'+m+'_vs_Gscore.png',bbox_inches='tight')
    plt.show()

#top_df.to_csv(output_dir+'meta-z_sum_all.txt',sep='\t',index=False)
#file_top.close()

