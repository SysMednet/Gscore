import numpy as np
from os import listdir
import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

filt = []
d = {}
d_pw = {}
pw = ['hsa03320','hsa04020','hsa04024','hsa04110'
      ,'hsa04115','hsa05206','hsa05207','hsa04010',
      'hsa04510','hsa04512']

#%%
with open ('D:\\Job\\handover\\Pan-cancer\\all_z_metaz_OR.txt') as f:
  lines = f.readlines()
  for line in lines:
    data = line.strip().split('\t')
    if data[0] not in d_pw and data[0] in pw:
      d_pw[data[0]] = data[1]
    if data[0] in pw and data[2] !='Gscore_v1':
      method = data[2].split('_')[0]
      if method not in d:
        d[method] = {}
      filt.append(line)
#%%
for f in filt:
  data = f.strip().split('\t')
  method = data[2].split('_')[0]
  if data[0] not in d[method]:
      d[method][d_pw[data[0]]] = float(data[-1])
      
#order = ['Calcium signaling pathway','MicroRNAs in cancer','cAMP signaling pathway','Chemical carcinogenesis - receptor activation','Cell cycle','p53 signaling pathway','PPAR signaling pathway','ECM-receptor interaction','Focal adhesion','MAPK signaling pathway']
order=['PPAR signaling pathway','Calcium signaling pathway','cAMP signaling pathway','Cell cycle','p53 signaling pathway','MicroRNAs in cancer','Chemical carcinogenesis - receptor activation','MAPK signaling pathway','Focal adhesion','ECM-receptor interaction']


df = pd.DataFrame(d)
df = df.reindex(order)
#df = df.sort_index()
#df = df.reindex(order_list)
df = df.replace(0, np.nan)
array = np.array(df)
r,c = array.shape
label = np.zeros((r,c),dtype=str)

for i in range(r):
  for j in range (c):
    if array[i][j] > 1.64:
      label[i][j] = 'O'
    else:
      label[i][j] = ''
      
sns.set(style='white')
ax = sns.heatmap(data=df,cmap='Reds',annot=label, fmt='',vmin=0,vmax=120,linewidths=1,linecolor='black',clip_on=False,annot_kws={"fontsize":15})
plt.xticks(fontsize=18, rotation=90)
plt.yticks(fontsize=18, rotation=0)
#plt.savefig('D:\\Job\\handover\\Pan-cancer\\Figs\\meta.png',dpi=300,bbox_inches='tight')
plt.show()