from scipy.stats import fisher_exact
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math


sns.set(context='notebook', style='ticks', font_scale=1.6) #2 sets
data_sets = ['GSE157103','GSE150316']
    
methods = ['Gscore_fdrq','NEA_merged6','ROntoTools','SPIA','PADOG','GSA','GSEA','ORA']
alt = 'greater'

#pastel
color1 = sns.color_palette('tab10')
color2 = sns.color_palette('pastel')
color = [color1[1],color1[0]]+color1[2:4]+color1[6:] #,color2[0],color2[1]

data_df = pd.read_csv('fig5A_input.txt', delimiter = "\t")
data_df2 = data_df
data_df2['set_method'] = list(data_df['DataSet']+'_'+data_df['methods'])

tmp_df1 = data_df2[data_df2['DataSet']=='GSE150316']
tmp_df2 = data_df2[data_df2['DataSet']=='GSE157103']
data_df2 = tmp_df2.append(tmp_df1)

ax1 = plt.figure(dpi = 300, figsize=(10,4.2)) #2 sets

#ax1 = sns.barplot(data=data_df2, x="DataSet", y="(-log10 P)",order=data_sets, hue="methods",edgecolor='black',palette=color)
ax1 = sns.barplot(data=data_df2, x="set_method", y="(-log10 P)",edgecolor='black',palette=color)
handles, labels = ax1.get_legend_handles_labels()
labels = ['Gscore','NEA']+labels[4:]
#plt.legend(handles,labels,bbox_to_anchor=(1.4,1.045),loc='upper right')

plt.xticks([],rotation=90,size=25)
plt.yticks([0,1,2,3],size=30)
ax1.set(xlim=(-0.9,2*8-0.1))

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.set_ylabel('Statistical\nsignificance\n $(-log_{10}p)$',fontweight ='bold',size=35) 
#ax1.set(title=dataset)

line_x = 7.5
for l in range(len(data_sets)-1):
    ax3 = ax1.twinx()
    ax3.axvline(x=line_x, ls='-', color='gray')
    ax3.axes.get_yaxis().set_visible(False)
    ax3.tick_params(top=False)
    line_x+=8

ax2 = ax1.twiny()
ax2.axhline(y=1.301, ls='--', color='black')
ax2.axes.get_xaxis().set_visible(False)
ax2.tick_params(top=False)
  
plt.savefig('fig5A.png',bbox_inches='tight')
plt.show()

       
        