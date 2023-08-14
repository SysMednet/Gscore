import os

file = 'For_paper_figure_347'
sets = 'TCGA'#75 or TCGA
cutoff_l = ['all']
#cutoff_l = ['400','250','500','1000']
for cutoff in cutoff_l:
    paths = [
            os.path.join(file, fname)
            for fname in os.listdir(file)
            if sets in fname
        ]
    
    methods = []
    for i in paths:
        tmp = i.split('\\')[-1].split('_')[0]
        methods.append(tmp)
        
    file1 = open('D:\Job\handover\impli_expli\KEGGpw\hsa_347_pathway.txt')
    file1_r = file1.read()
    file1.close()
    
    if sets=='75':
        file2 = open('GEO75_chosen_set.txt')
        file2_r = file2.read()
        file2.close()
        
        file3 = open('75_related+target_pathways.txt')
        file3_r = file3.read()
        file3.close()
    else:
        file2 = open('TCGA_list.txt')
        file2_r = file2.read()
        file2.close()
        
        file3 = open('TCGA_related+target_pathways.txt')
        file3_r = file3.read()
        file3.close()
    
    pathways = file1_r.split('\n')[:-1]
    data_sets = file2_r.split('\n')[:-1]
    data3 = file3_r.split('\n')[:-1]
    
    related = {}
    for i in data3:
        tmp = i.split('\t')
        related[tmp[0]] = tmp[1:]
    
    result_dic = {}
    
    for i in data_sets:
        result_dic[i] = {}
        for j in pathways:
            if i not in related:
                np = 'NA'
            elif j in related[i]:
                np = 'P'
            else:
                np = 'N'
            result_dic[i][j] = {}
            result_dic[i][j]['NP'] = np
            
            for k in methods:
                result_dic[i][j][k] = '-'
        
    for i in paths:
        me = i.split('\\')[-1].split('_')[0]
        
        if 'gscore' in i:
            paths2 = [os.path.join(i, fname)
                      for fname in os.listdir(i)]
            
            for j in paths2:
                data_set = j.split('\\')[-1].split('_KEGG_')[0]
                
                filed = open(j)
                filed_r = filed.read()
                filed.close()
                
                data = filed_r.split('\n')[1:-1]
                for k in data:
                    if 'all_DEG' in k and 'up' not in k and 'down' not in k:
                    #if '_'+cutoff in k:
                        tmp = k.split('\t')
                        ID = tmp[0]
                        pv = tmp[-1]
                        
                        result_dic[data_set][ID]['gscore'] = pv
                        
        elif 'ORA' in i and 'cepa' not in i:
            paths2 = [os.path.join(i, fname)
                      for fname in os.listdir(i)]
            
            for j in paths2:
                data_set = j.split('\\')[-1].split('_ORA_')[0]
                
                filed = open(j)
                filed_r = filed.read()
                filed.close()
                
                data = filed_r.split('\n')[1:-1]
                for k in data:
                    if 'all_DEG' in k and 'up' not in k and 'down' not in k:
                    #if '_'+cutoff in k:
                        tmp = k.split('\t')
                        ID = tmp[0]
                        pv = tmp[-1]
                        
                        result_dic[data_set][ID]['ORA'] = pv
                        
        elif 'SPIA' in i or 'ronto' in i or 'NEA' in i:
            paths2 = [os.path.join(i, fname)
                      for fname in os.listdir(i)
                      if cutoff in fname]
            
            for j in paths2:
                data_set = j.split('\\')[-1].split('_output_')[0]
                
                filed = open(j)
                filed_r = filed.read()
                filed.close()
                
                data = filed_r.split('\n')[1:-1]
                for k in data:
                    #if 'NA' not in k:           
                        tmp = k.split('\t')
                        ID = tmp[0]
                        pv = tmp[-1]
                        
                        result_dic[data_set][ID][me] = pv
                        
        elif 'CePa' not in i and ('GSEA' in i or 'padog' in i or 'GSA' in i):
            paths2 = [os.path.join(i, fname)
                      for fname in os.listdir(i)]
            
            for j in paths2:
                data_set = j.split('\\')[-1].split('_output_')[0]
                
                filed = open(j)
                filed_r = filed.read()
                filed.close()
                
                data = filed_r.split('\n')[1:-1]
                for k in data:
                    #if 'NA' not in k:           
                        tmp = k.split('\t')
                        ID = tmp[0]
                        pv = tmp[-1]
                        result_dic[data_set][ID][me] = pv   
    
    fileo = open('D:\\Job\\handover\\impli_expli\\result\\NA\\fdr\\'+sets+'_347_'+cutoff+'_fdrpv.txt','w')
    output = 'Data_set'+'\t'+'pathway'+'\t'+'Gscore'+'\t'+'NEA'+'\t'+'ROntoTools'+'\t'+'SPIA'+'\t'+'GSA'+'\t'+'GSEA'+'\t'+'PADOG'+'\t'+'ORA'+'\n'
    
    for i in result_dic:
        for j in result_dic[i]:
            output+= i+'\t'+j+'\t'+result_dic[i][j]['gscore']+'\t'+result_dic[i][j]['NEA']+'\t'+result_dic[i][j]['ronto']+'\t'+result_dic[i][j]['SPIA']+'\t'+result_dic[i][j]['GSA']+'\t'+result_dic[i][j]['GSEA']+'\t'+result_dic[i][j]['padog']+'\t'+result_dic[i][j]['ORA']+'\n'
    fileo.write(output)
    fileo.close()      
        
