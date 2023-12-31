import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import powerlaw
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data_sets = {'GSE157103':' '*12,'GSE150316':''}
xtick = ['Gscore','NEA','ROntoTools','SPIA','PADOG','GSA','GSEA','ORA']*len(data_sets)
xlabel = ' '*8

for d in data_sets:
    xlabel+=d+data_sets[d]

size = []
size_dic = {}
scale_size = {1:850,5:780,10:700,20:550,30:400,50:250,100:150,200:100,300:70,347:50}
num_list = range(1,347)
    
size = [850]
rank_list = list(scale_size.keys())
for i in range(1,len(rank_list)):
    s = rank_list[i]-rank_list[i-1]
    size+=[scale_size[rank_list[i]]]*s

for i in range(len(size)):
    size_dic[i+1] = size[i]*2.5

sns.set(context='notebook', style='ticks', font_scale=4)

data_df = pd.read_csv('fig5B_input.txt', delimiter = "\t")
num_of_path = data_df['path_name'].describe().unique()[1]

ax_com = plt.figure(dpi = 300, figsize=(20,num_of_path+2))
mycol0 = sns.color_palette("coolwarm")
mycol = []
for c in mycol0:
    mycol.append(c)
mycol.reverse()
ax_com = sns.scatterplot(data=data_df, x="methods_dataset", y="path_name", hue="adjp",palette="coolwarm_r",size="Rank",sizes=size_dic)
ax_com.set_xticklabels(xtick)

ax_com.set_xlabel(xlabel,fontweight ='bold',fontsize=60,loc="left")
ax_com.set(xlim=(-0.9,2*8-0.1))

ax_com.set(ylabel='')
plt.xticks(rotation=90,fontsize=55)

handles, labels = ax_com.get_legend_handles_labels()
handles = handles[labels.index('Rank'):]
labels = labels[labels.index('Rank'):]
scale = [1,10,20,25,30,len(handles)-1]
handles2 = [handles[0]]
labels2 = [labels[0]]

for hl in scale:
    handles2.append(handles[hl])
    labels2.append(int(float(labels[hl])))

plt.legend(handles=handles2,labels=labels2 ,bbox_to_anchor=(1.27,1.02),loc='upper right',handleheight=2) #17

max_num = max(data_df['adjp'])
adjp_min = [x for x in list(data_df['adjp']) if str(x)!='nan']
adjp_min.sort()
norm = plt.Normalize(adjp_min[0],max_num)
sm = plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm)

##(x,y,寬度,高度)
axins = ax_com.inset_axes((1.03,-0.25,0.03,0.45)) #2 sets
clb = plt.colorbar(sm,cax=axins, orientation="vertical")
clb_text = {'fontsize':52}
clb.ax.set_title('FDR $\it{q}$ value/\nNominal $\it{p}$ value', loc='left',**clb_text)

line_y = 0.5
for l in range(num_of_path-1):
    ax2 = ax_com.twiny()
    ax2.axhline(y=line_y, ls='--', color='gray')
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(top=False)
    line_y+=1

line_x = 7.5
for l in range(len(data_sets)-1):
    ax2 = ax_com.twinx()
    ax2.axvline(x=line_x, ls='-', color='Black')
    ax2.axes.get_yaxis().set_visible(False)
    ax2.tick_params(top=False)
    line_x+=8

plt.savefig('fig5B.png',bbox_inches='tight')
plt.show()