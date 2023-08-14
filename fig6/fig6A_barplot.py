from scipy.stats import fisher_exact
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

sns.set(context='notebook', style='ticks', font_scale=1.5)

data_sets = ['GSE157103','GSE150316']
methods = ['Gscore_fdrq','NEA_merged6']
tar_gene = ['top_5%','top_1%','896']
label = ['Top 5%','Top 1%','DisGeNET &\nUniProtKB']
alt_type = ['allDEG']
methods_label = ['Gscore','NEA']*len(data_sets)*len(label)

order_list = []
for d in data_sets:
    for alt in alt_type:
        order_list.append(d+'_'+alt)
        
xlabel = ' '*5
for d in data_sets:
    xlabel+=d+(' '*8)

#palette
color1 = sns.color_palette('tab10')
color2 = sns.color_palette('pastel')
#color = [color1[1],color2[1],color1[0],color2[0]]+color1[2:4]+color1[6:]
color = [color1[1],color1[0]]

figure_df = pd.DataFrame()
for tar,l in zip(tar_gene,label):
    
    for alt in alt_type:
        data_df = pd.read_csv('single/COVID_'+alt+'_'+tar+'.txt', delimiter = "\t")
        
        if figure_df.empty:
            figure_df = data_df[((data_df['methods']=='Gscore_fdrq') | (data_df['methods']=='NEA_merged6')) & ((data_df['DataSet']=='GSE157103') | (data_df['DataSet']=='GSE150316'))]
            figure_df['gene'] = [l,l,l,l]
        else:
            data_tmp = data_df[((data_df['methods']=='Gscore_fdrq') | (data_df['methods']=='NEA_merged6')) & ((data_df['DataSet']=='GSE157103') | (data_df['DataSet']=='GSE150316'))]
            data_tmp['gene'] = [l,l,l,l]
            figure_df = figure_df.append(data_tmp,ignore_index = True)
    
    figure_df['x'] = list(figure_df['methods']+'_'+figure_df['gene'])
    #figure_df['methods_label'] = methods_label

for d in data_sets:
    figure_df2 = figure_df[figure_df.DataSet==d]
    
    ax1 = plt.figure(dpi = 300, figsize=(7,3))
    ax1 = sns.barplot(data=figure_df2, x="gene", y="-log10P", hue="methods",edgecolor='black',palette=color)
    ax1.legend([],[], frameon=False)

    ax1.set_ylabel('')
    ax1.set_xlabel('COVID-19-related genes',fontweight ='bold',size=23)
    #plt.legend(bbox_to_anchor=(1.5,1),loc='upper right', prop={'size': 25})
    
    plt.yticks(size=20)
    plt.xticks(size=20)
    ax1.tick_params(bottom=False)
    ax1.set_title(d,fontweight ='bold',size=25)
    
    ax2 = ax1.twiny()
    ax2.axhline(y=1.301, ls='--', color='black')
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(top=False)
    
    plt.savefig('single/'+d+'_allDEG.png',bbox_inches='tight')
    plt.show()
        
