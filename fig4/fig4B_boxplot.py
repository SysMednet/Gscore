import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy
from statannotations.Annotator import Annotator

table = pd.read_csv('fig4B_input.txt','\t',index_col=0)

pre_df = {'value':[],
          'method':[],
          'score_type':[]}

for value in ['precision','accuracy', 'FPR', 'recall']:
    values={"Gscore":[], "NEA":[]}
    print(value)
    for i in range(len(table['cutoff'])):
        if table['cutoff'][i] == 0.05:
            if table[value][i] != '-': #Gscore
                pre_df['value'].extend([float(table[value][i])])
                pre_df['method'].extend(['Gscore'])
                values["Gscore"].append(float(table[value][i]))
                if value != 'FPR':
                    score_type = value.capitalize()
                    pre_df['score_type'].extend([score_type])
                else:
                    score_type = 'Fale positive\nrate'
                    pre_df['score_type'].extend([score_type])
                
                    
                
            if table[value+'.1'][i] != '-': #NEA
                pre_df['value'].extend([float(table[value+'.1'][i])])
                pre_df['method'].extend(['NEA'])
                values["NEA"].append(float(table[value+'.1'][i]))
                if value != 'FPR':
                    score_type = value.capitalize()
                    pre_df['score_type'].extend([score_type])
                else:
                    score_type = 'Fale positive\nrate'
                    pre_df['score_type'].extend([score_type])
df = pd.DataFrame(pre_df)

#%%
plt.figure(dpi=300,figsize=(8,6))
sns.set(context='notebook', style='ticks', font_scale=1.9)
order = ['Accuracy','Precision', 'Recall', 'Fale positive\nrate']
color1 = sns.color_palette("tab10")
color_list = []
color_list.append(color1[1])
color_list.append(color1[0])
fig = sns.boxplot(data=df, x='score_type', y = 'value',hue='method',order=order, showfliers = False, palette=color_list, linewidth=2)
#fig.set_title(dataset_type+'_'+path_num+'_Wilcoxon')
plt.xticks(rotation=0)
plt.yticks([0,0.25,0.5,0.75,1.0])
plt.ylim(0.0, 1.4)
fig = sns.swarmplot(data=df,x='score_type', y = 'value',hue='method',order=order,dodge = True, color='.2', size=5)
handles, tmp = fig.get_legend_handles_labels()
plt.legend(handles=handles[:-2],labels=tmp[:-2])
fig.set_ylabel('Percentage (*100%)')
fig.set_xlabel('')

#%%
x = 'score_type'
y = 'value'

pairs=[(('Precision','Gscore'), ('Precision','NEA')), (('Accuracy','Gscore'), ('Accuracy','NEA')), (('Fale positive\nrate','Gscore'), ('Fale positive\nrate','NEA')), (('Recall','Gscore'), ('Recall','NEA'))]

annotator = Annotator(fig, pairs, data=df, x=x, y=y, order=order,hue='method')
annotator.configure(test='Wilcoxon', text_format='star', loc='inside',pvalue_thresholds=[[1e-3, "***"], [1e-2, "**"], [0.05, "*"], [1, "ns"]])
annotator.apply_and_annotate()

plt.savefig('fig4B.png', bbox_inches='tight')
plt.show()