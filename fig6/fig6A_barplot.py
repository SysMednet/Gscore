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

#palette
color1 = sns.color_palette('tab10')
color2 = sns.color_palette('pastel')
color = [color1[1],color1[0]]

figure_df = pd.read_csv('fig6A_input.txt', delimiter = "\t")

for d in data_sets:
    figure_df2 = figure_df[figure_df.DataSet==d]
    
    ax1 = plt.figure(dpi = 300, figsize=(7,3))
    ax1 = sns.barplot(data=figure_df2, x="gene", y="-log10P", hue="methods",edgecolor='black',palette=color)
    ax1.legend([],[], frameon=False)

    ax1.set_ylabel('')
    ax1.set_xlabel('COVID-19-related genes',fontweight ='bold',size=23)
    
    plt.yticks(size=20)
    plt.xticks(size=20)
    ax1.tick_params(bottom=False)
    ax1.set_title(d,fontweight ='bold',size=25)
    
    ax2 = ax1.twiny()
    ax2.axhline(y=1.301, ls='--', color='black')
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(top=False)
    
    plt.savefig('fig6A_'+d+'.png',bbox_inches='tight')
    plt.show()
        
