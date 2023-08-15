import scipy.stats as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

xtick = ['Gscore','NEA','ROntoTools','SPIA','PADOG','GSA','GSEA','ORA']

sns.set(context='notebook', style='ticks',font_scale=4.5)
#palette
color1 = sns.color_palette('tab10')
color2 = sns.color_palette('pastel')
color = [color1[1],color1[0]]+color1[2:4]+color1[6:] #,color2[1]

op_dir = 'fig3A.png'

data_df = pd.read_csv('fig3A_input.txt', delimiter = "\t",index_col=0)
methods = list(data_df.keys())
total_tar = len(data_df)

result_dic = {'methods':[],'counter':[]}

for m in methods:
    if m=='GSEA':
        enrich_value=abs(st.norm.ppf(0.05)) #0.25
    elif m=='GSA':
        enrich_value=abs(st.norm.ppf(0.05)) #0.1
    else:
        enrich_value=abs(st.norm.ppf(0.05))
    
    result_dic['methods'].append(m)
    counter = 0
    
    for d in data_df[m]:
        if d>=enrich_value:
            counter+=1
    result_dic['counter'].append((counter/total_tar)*100) 
result_df = pd.DataFrame.from_dict(result_dic)

counter_list = list(result_df['counter'])

ax1 = plt.figure(dpi = 300, figsize=(10,17))
ax1 = sns.barplot(data=result_df , x="counter", y="methods" ,edgecolor='black' ,palette=color ,orient='h')

############改bar寬度(橫的)#############
for patch in ax1.patches:
    current_width = patch.get_height()
    diff = current_width - 0.5
    patch.set_height(0.5)
    patch.set_y(patch.get_y() + diff * .5)
#######################################

'''ax1.set_xticklabels(xtick)
plt.xticks(rotation=20)
ax1.set(title=o,ylabel='',xlabel='') #Enriched pan-cancer-related pathways (%)'''

ax1.set_yticklabels(xtick)
plt.xlabel('')
plt.ylabel('')
plt.xticks([0,25,50,75,100])
#ax1.set(title=o,ylabel='',xlabel='') #Enriched pan-cancer-related pathways (%)
        
'''for c in range(len(counter_list)):
    if counter_list[c]<15:
        ax1.text(counter_list[c]+8, c+0.1, str(round(counter_list[c],2)),size=37,ha='center',weight='bold')
    else:
        ax1.text(counter_list[c]-11, c+0.1, str(round(counter_list[c],2)),size=37,ha='center',weight='bold')'''

plt.savefig(op_dir,bbox_inches='tight')
plt.show()