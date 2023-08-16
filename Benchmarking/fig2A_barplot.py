import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
plot_d = {'75':{},'TCGA':{}}
dtypes = ['75','TCGA']

non_enrich = '#B0C4DE'
enrich = '#FFBF00'
order_l = ['Gscore','NEA','ROntoTools','SPIA','CePaORA','CePaGSA','PADOG','GSA','GSEA','ORA']
methods = {'Gscore':0,'NEA':1,'ROntoTools':2,'SPIA':3,'cepaORA':4,'CepaGSA':5,'PADOG':6,'GSA':7,'GSEA':8,'ORA':9}

for k in plot_d:
  with open ('fig2A_input.txt') as f:
    lines =f.readlines()
    for line in lines[1:]:
      data = line.strip().split('\t')
      if data[2] not in plot_d[k]:
        plot_d[k][data[2]] = [0,0,0,0,0,0,0,0,0,0]
      for m in methods:
        if data[0] == m and data[1]==k:
          plot_d[k][data[2]][methods[m]] = float(data[3])*100

  plot_df = pd.DataFrame(plot_d[k],index=order_l)
  ax = plot_df.plot(kind='bar',stacked=True,color=[enrich,non_enrich,'gray'],edgecolor='black',fontsize=15)
  plt.title('')
  plt.xlabel('Methods',fontsize=15)
  plt.ylabel('Percentage (%)',fontsize=15)
  ax.legend(bbox_to_anchor=(1.35,0.6),fontsize=10)
  plt.savefig('fig2A.png',bbox_inches='tight',dpi=300)
  plt.show()
    
