res = []
with open ('D:\\Job\\handover\\Fig2\\data_for_fig\\percentage_all_methods\\association_count.txt') as f:
  lines = f.readlines()
  for line in lines[1:]:
    data = line.strip().split('\t')
    if data[1] == '75':
      dv = 75*347
    elif data[1] == 'TCGA':
      dv = 16*347    
    per = int(data[3])/dv
    data[-1] = str(per)
    data.append(str(1)+'\n')
    res.append('\t'.join(data))

with open ('D:\\Job\\handover\\Fig2\\data_for_fig\\percentage_all_methods\\association_per.txt','w') as op:
  op.write('method\ttype\tcount\tpercentage\tall\n')
  for i in res:
    op.write(i)