import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

file_name = {'fig2B_input_GEO_rank.txt':['fig2B_GEO_rank.png',1,'GEO Rank'],
             'fig2B_input_TCGA_rank.txt':['fig2B_TCGA_rank.png',1,'TCGA Rank'],
             }

color = sns.color_palette('tab10')
tmp = color[0]
color[0] = color[1]
color[1] = tmp

sns.set(context='notebook', style='ticks', font_scale=2.7)
text_size = 30

#color1[0] = (1.0, 1.0, 0)

for f_name in file_name:
    v_list = file_name[f_name]
    
    data_df = pd.read_csv(f_name, delimiter = "\t")
    data_melt = pd.melt(data_df)
    medians = data_melt.groupby('variable', sort=False)['value'].median()
    
    ax = plt.figure(dpi = 300, figsize=(20,7))
    sns.set_style("whitegrid")
    ax = sns.violinplot(data=data_melt, x="variable", y="value",scale = "width",width=0.4,inner=None,cut=0,palette=color)
    
    plt.scatter(x=range(len(medians)),y=medians,c="k")
    
    for i in range(len(medians)):
        left = i-0.3
        right = i+0.3
    
        median_df = pd.DataFrame({'x_value':[left,right],'y_value':[medians[i],medians[i]]})
        sns.lineplot(data=median_df,x='x_value', y='y_value',color='black')
        
        if v_list[1]==1:
            ax.text(i, medians[i]+3, medians[i],size=text_size,ha='center',weight='bold')
        else:
            if medians[i]<0.01:
                ax.text(i, medians[i]+0.03, "%.1e"%medians[i],size=text_size,ha='center',weight='bold')
            else:
                if medians[i]>=0.97:
                    ax.text(i, medians[i]-0.07, round(medians[i],2),size=text_size,ha='center',weight='bold')
                else:
                    ax.text(i, medians[i]+0.03, round(medians[i],2),size=text_size,ha='center',weight='bold')
    
    if v_list[1]==1:
        ax.set(xlabel='',ylabel='Ranks of target pathways',ylim=(0,186)) #,title=v_list[2]
    else:
        ax.set(xlabel='',ylabel='P-values of target pathways',ylim=(-0.03,1.03)) #,title=v_list[2]
    ax.set_ylabel(ax.get_ylabel(), fontdict={'weight': 'bold'})
    counter=0
    for tc in ax.get_xticklabels():
        if counter<6:
            tc.set_color("black") #red
        else:
            tc.set_color("black") #blue
        counter+=1 
    ax.set_xticklabels(ax.get_xticklabels())
        
    plt.xticks(rotation=30)
    
    '''ax2 = ax.twinx()
    ax2.axvline(x=5.5, ls='-', color='grey')
    ax2.axes.get_yaxis().set_visible(False)
    ax2.tick_params(top=False)
    
    line_x = 0.5
    for l in range(9):
        ax2 = ax.twinx()
        ax2.axvline(x=line_x, ls='--', color='grey')
        ax2.axes.get_yaxis().set_visible(False)
        ax2.tick_params(top=False)
        line_x+=1'''
    
    
    plt.savefig(v_list[0],bbox_inches='tight')
    plt.show()
