##python3 correlation.py <input.path>
import numpy as np

np.seterr(divide='ignore', invalid='ignore')

def correlation(input_path,output_path):
    #read gene expression data
    fp=open(input_path+'/'+'GSE157103.txt','r',encoding='utf-8')
    lines = fp.readlines()
    fp.close()
    
    #read DEG list
    fp2=open(input_path+'/'+'GSE157103_DEG.txt','r',encoding='utf-8')
    DEGs = fp2.read().strip().split('\n')
    fp2.close()
    DEGs = set(DEGs)
    
    #extract DEG case(phenotype 1) samples' data
    case_sample = []
    phenotype = lines[0].strip().split('\t')
    for c in range(len(phenotype)):
        if phenotype[c]=='1':
            case_sample.append(c)
    
    Exp , geneID = [] , []
    for line in lines[1:]:
        temp=line.strip().split('\t')
        if temp[0] in DEGs:
            geneID.append(temp[0])
            for i in case_sample:
                Exp.append(temp[i])
    Exp = np.array(Exp,dtype=float).reshape(len(geneID),len(case_sample))
    PCC = np.round(np.corrcoef(Exp),6) #calculate pcc
    
    op=open(output_path+'/'+'GSE157103_pcc.txt','w',encoding='utf-8')
    op.write('gene1\tgene2\tcorrelation\n')
    for i in range(len(geneID)):
        for j in range(i,len(geneID)):
            op.write(geneID[i]+'\t'+geneID[j]+'\t'+str(PCC[i][j])+'\n')
    op.close()
    
if __name__ == '__main__' :
    input_path = 'sample_data_input/'
    output_path = 'sample_data_output/'
    correlation(input_path,output_path)
