import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

m = 'Gscore_fdrq'
data_sets = ['GSE150316','GSE157103']

#color for gene in database heatmap
cmap_dict = {1: '#EE7700', 3: '#003C9D', 2: '#66009D'}
cmap = ListedColormap([cmap_dict[i] for i in range(1,4)])

heat_color = sns.color_palette('Reds')[1:]

for d in data_sets:
    data_df = pd.read_csv('fig6B_input/COVID_single_z_'+m+'_'+d+'.txt', delimiter = "\t",index_col=0)
    data_df = data_df.loc[:,(data_df != 0).any(axis=0)]
    data_df.replace(0, np.nan, inplace=True)
    count_df = pd.read_csv('fig6B_input/COVID_single_count_'+m+'_'+d+'.txt', delimiter = "\t",index_col=0)
    exist_df = pd.read_csv('fig6B_input/COVID_single_exist_'+m+'_'+d+'.txt', delimiter = "\t",index_col=0)
    exist_df.index = list(exist_df['gene'])
    exist_df = exist_df.drop(['gene'], axis=1)
    exist_df.replace(0, np.nan, inplace=True)

    sns.set(font_scale=3, style='white')
    edge_color = 'gray'
    
    #heatmap
    ax2 = plt.figure(dpi = 300,figsize=(25,len(data_df)*0.8)) #len(data_df)*0.3
    
    ax2 = sns.heatmap(data_df, cmap=heat_color, cbar=False)
    ax2.axhline(y=0, color=edge_color,linewidth=6)
    ax2.axhline(y=data_df.shape[0], color=edge_color,linewidth=6)
    ax2.axvline(x=0, color=edge_color,linewidth=6)
    ax2.axvline(x=data_df.shape[1], color=edge_color,linewidth=6)
    plt.yticks(fontweight='bold',size=50)
    plt.xticks(rotation=90)

    ax2.set_title('')
    plt.savefig('fig6B_output/COVID_single_z_'+m+'_'+d+'.png',bbox_inches='tight')
    plt.show()
    
    #count barplot
    ax1 = plt.figure(dpi = 300,figsize=(3,len(count_df)*0.8))
    
    ax1 = sns.barplot(data=count_df, x="count", y="gene" , orient = 'h',color='gray')
    plt.xticks([0,70,140],fontweight='bold',rotation=90,size=50)
    ax1.set_yticklabels([])
    plt.ylabel('')
    plt.xlabel('')
    plt.savefig('fig6B_output/COVID_single_count_'+m+'_'+d+'.png',bbox_inches='tight')
    plt.show()
    
    #gene in database heatmap
    ax3 = plt.figure(dpi = 300,figsize=(3,len(count_df)*0.8))
        
    ax3 = sns.heatmap(exist_df,cbar=False,cmap=cmap,linewidths=10, linecolor='white')
    ax3.axhline(y=0, color=edge_color,linewidth=5)
    ax3.axhline(y=exist_df.shape[0], color=edge_color,linewidth=5)
    ax3.axvline(x=0, color=edge_color,linewidth=5)
    ax3.axvline(x=exist_df.shape[1], color=edge_color,linewidth=5)
    
    ax3.set(yticklabels=[],ylabel='',xticklabels=[],xlabel='')
    plt.savefig('fig6B_output/COVID_single_exist_'+m+'_'+d+'.png',bbox_inches='tight')
    plt.show()
#%%
#create heatmap legend
import matplotlib as mpl

fig = plt.figure(dpi=300)
ax = fig.add_axes([0.05, 0.80, 0.6, 0.07]) #(x,y,寬度,高度)
cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal',
                       ticks=[0,1,2,3,4,5,6],
                       cmap='Reds',
                       norm=mpl.colors.Normalize(0, 6.5),  # vmax and vmin
                       ticklocation='bottom')
cb.ax.tick_params(labelsize=25)
cb.set_label(label='Association\nSignificance (z-score)',size=25, weight='bold')
plt.savefig('fig6B_output/colorbar.png', bbox_inches='tight')
plt.show()
