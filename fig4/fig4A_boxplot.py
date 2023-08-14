import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

path = ''
for dataset_type in ['16TCGA']:
    for pw_num in ['347']:
        table = pd.read_csv(path+'5A_Allto'+pw_num+'_'+dataset_type+'_FDRQ_delm0.txt','\t',index_col=0)
        
        pre_df = {'value':[],
                  'label':[]}
        for i in range(len(table['Gscore_enriched_%'])):
            pre_df['value'].extend([float(table['Gscore_enriched_%'][i])])
            pre_df['label'].extend(['Gscore'])
            pre_df['value'].extend([float(table['NEA_enriched_%'][i])])
            pre_df['label'].extend(['NEA'])

        df = pd.DataFrame(pre_df)   # 將資料組合成名稱為df的data frame
        
        sns.set(context='notebook', style='ticks')
        #sns.set_context("poster", font_scale = 2, rc={"grid.linewidth": 1})
        plt.figure(figsize=(3, 5))
        order = ['Gscore', 'NEA']
        color1 = sns.color_palette("tab10")
        color_list = []
        for i in [1, 0]:
            color_list.append(color1[i])
        for i in range(2,len(color1)):
            color_list.append(color1[i])
        fig = sns.boxplot(data=df, x='label', y = 'value', order=order, showfliers = False, palette=color_list,linewidth=1.5)
        #fig.set_title(dataset_type+'_'+pw_num, size=52)
        fig.set_title('')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        sns.swarmplot(data=df,x='label', y = 'value', color='.2', size=5)
        fig.set_ylabel('') #Percentage (*100%)
        fig.set_xlabel('')
                
        x = 'label'
        y = 'value'
        pairs=[('Gscore', 'NEA')]
        annotator = Annotator(fig, pairs, data=df, x=x, y=y, order=order)
        annotator.configure(test='Wilcoxon', text_format='star', loc='inside', fontsize=25, line_width=1.5, pvalue_thresholds=[[1e-3, "***"], [1e-2, "**"], [0.05, "*"], [1, "ns"]])
        annotator.apply_and_annotate()
        plt.savefig('new/'+dataset_type+'_'+pw_num+'_FDRQ_delm0.png', dpi=300, bbox_inches='tight')
        plt.show()
        