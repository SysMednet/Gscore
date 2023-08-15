##python3 correlation.py <input.path>
import numpy as np
import argparse

np.seterr(divide='ignore', invalid='ignore')

def correlation(input_GEM, input_DEG_list, output_path):
    #read gene expression data
    fp=open(input_GEM,'r',encoding='utf-8')
    lines = fp.readlines()
    fp.close()
    
    #read DEG list
    fp2=open(input_DEG_list,'r',encoding='utf-8')
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
    
    op=open(output_path,'w',encoding='utf-8')
    op.write('gene1\tgene2\tcorrelation\n')
    for i in range(len(geneID)):
        for j in range(i,len(geneID)):
            op.write(geneID[i]+'\t'+geneID[j]+'\t'+str(PCC[i][j])+'\n')
    op.close()
    
if __name__ == '__main__' :
    parser = argparse.ArgumentParser(description="Manual")
    parser.add_argument("-e", type=str , default='./sample_data_input/example_GEM.txt' , help='The gene expression profile.')
    parser.add_argument("-g", type=str , default='./sample_data_input/example_allDEG.txt' , help='The DEG list file, the DEG IDs are seperated by "\n".')
    parser.add_argument("-o", type=str , default='./sample_data_output/example_pcc.txt' , help='The output file name for Pearson correlation network.')
    
    args = parser.parse_args()
    input_GEM, input_DEG_list, output_path = args.e, args.g, args.o
    
    correlation(input_GEM, input_DEG_list, output_path)
