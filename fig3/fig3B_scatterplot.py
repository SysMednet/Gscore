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

data_df = pd.read_csv('fig3B_input.txt', delimiter = "\t")
sns.set(font_scale=3)

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


