# import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
# import seaborn as sns
# from pylab import *
from collections import Counter
import math
import sys
import os
import random
import seaborn as sns

def timeuse(func) :
    def wrapper(*args,**kwargs) :
        import datetime
        start = datetime.datetime.now()
        res = func(*args,**kwargs)
        finish = datetime.datetime.now()
        print("%s\tfinish time:%s"%(func.__name__+'\t'+'\t'.join(str(s) for s in args),str(finish - start)))
        return res
    return wrapper

@timeuse
def ks_pic(file_list,label_list,target_list,file_save,onlylabel):

    splitrange = 500
    
    filelen = len(file_list)
    ranklist , scorelist = {} , {}
    recall = {}
    
    for label in label_list :
        ranklist[label] = np.zeros(500+1)
        scorelist[label] = np.zeros(500)
        recall[label] = np.zeros(2,dtype=int)
    recall["random"]= np.zeros(2,dtype=int)
        
    ranscore , ranrank = np.zeros(500) , np.zeros(500+1)

    for file , l ,label,targetfile in zip(file_list,range(filelen),label_list,target_list) :

        #open target file
        target = set()
        with open(targetfile,mode='r') as r_line :
            for n_line in r_line :
                tem = n_line.strip('\n').split('\t')
                target.add(tem[0])
                
        gene , value = [] , []
        genescore = []
        rangenescore = []
        check = True
        with open(file,mode='r') as r_line :
            for n_line in r_line.readlines()[1:]:   #跳過第一行
                g ,qval,*val = n_line.strip('\n').split('\t')    #JunMao revised(qval)
                # print( g ,qval,val)
                value+=val
                gene.append(g)
                if check :
                    patlen = len(val)
                    check = False
        genelen = len(gene)
        value = np.array(value,dtype=float).reshape(genelen,patlen) 
        gene = np.array(gene,dtype=str)
        
        randommax = genelen+1
        
        for r in range(patlen) :
            val = value[:,r]
            loc = np.where(val>0)
            val = val[loc]
            vmin , vmax = np.min(val) , np.max(val)
            val = (val-vmin)/(vmax-vmin)
            val[np.isnan(val)] = 0
            # print(len(loc))
            recall[label][1] += len(set(target))
            recall["random"][1] +=len(set(target))   
            for g , v in zip(gene[loc],val) :
                if g in target :
                    genescore.append(v)
            
            #random
            ran = random.sample(range(1,randommax,1),genelen)
            ran = (np.array(ran)-1)/randommax
            for g , v in zip(gene,ran) :
                if g in target :
                    rangenescore.append(v)
 
        recall[label][0] += len(genescore)
        recall["random"][0] +=len(rangenescore)
       
                    
        rangec = 1/splitrange
        score, rank = np.histogram(genescore,bins=np.arange(0,1+rangec,rangec))

        ranklist[label] = rank
        scorelist[label] += score
        
        #random
        rantscore, rantrank = np.histogram(rangenescore,bins=np.arange(0,1+rangec,rangec))
        ranscore += rantscore
        ranrank = rantrank
        
    
    plt.figure(figsize=(14,10),dpi=300,linewidth = 3)
    sns.set(context='notebook', style='ticks', font_scale=1.3)
    file_w=open(f"{file_save}_auc_PCCchanged.txt","w")#
    file_w.write("method\tAUC\n")#auc_writte
    for label in onlylabel :
        xvalue = ranklist[label]
        yvalue = scorelist[label]
        csum = np.sum(yvalue)
        fscore = yvalue/csum 
        cumscore = np.cumsum(fscore)
        cumscore = np.r_[0,cumscore] - rank
        if 'fav' in label :
            sty = '-'
        elif 'adv' in label :
            sty = '-'
        else :
            sty = ':'
        
        if 'Gscore' in label and 'fav' in label :
            color_pick=sns.color_palette("tab10")[1]
        elif 'NEA' in label and 'fav' in label :
            color_pick=sns.color_palette("tab10")[0]
        elif 'Gscore' in label and 'adv' in label :
            color_pick=sns.color_palette("tab10")[1]
        elif 'NEA' in label and 'adv' in label :
            color_pick=sns.color_palette("tab10")[0]
        else :
            color_pick='orange'
        area = np.trapz(cumscore,xvalue)
        plt.plot(xvalue,cumscore,sty,color=color_pick,label=f'{label},AUC={area:.4}')
        file_w.write(label+"\t"+str(area)+"\n")#auc_writte
    
    #random
    rancsum = np.sum(ranscore)
    ranfscore = ranscore/rancsum 
    rancumscore = np.cumsum(ranfscore)
    rancumscore = np.r_[0,rancumscore] - ranrank
    area = np.trapz(rancumscore,ranrank)    
    plt.plot(ranrank,rancumscore,'-',color='black',label=f'random,AUC={area:.2}')
    file_w.write("random"+"\t"+str(area)+"\n")#auc_writte
    file_w.close()#auc_writte
    

    
        
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel("Rank", fontsize=40, labelpad = 15)
    plt.ylabel("Cumulative Probability-Rank", fontsize=30, labelpad = 20)
    # plt.legend(loc = "best", fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left',fontsize=20)
    plt.tight_layout()
    plt.savefig(file_save+".tiff")
    plt.close()

    with open(f'{file_save}_PCCchanged.txt',mode='w') as w_line :
        w_line.write('label\ttarget\tall_gene\n')
        for label in recall :
            w_line.write(label+'\t'+'\t'.join(str(s) for s in recall[label])+'\n')


#%%16個癌症整合
for dataset_type in ['TCGA']:
    #for path_cond in ['TR', 'OandTR', 'ORandTR']:
    for path_cond in ['OR']:
        for pw_num in ['347']:
            for prognostic_type in ['fariable', 'adverse']:
                file_list = []
                label_list = []
                onlylabel = []
                target_list = []
                
                if pw_num == '149':
                    f = '/data1/LanYun/GScore/dataset_and_pw/'+'10TCGA_target_pw.txt'
                if pw_num == '347':
                    f = '/data1/LanYun/GScore/dataset_and_pw/'+'13TCGA_target_pw.txt'
                    
                with open(f,mode='r') as rline :
                    for nline in rline :
                        tem = nline.strip('\n').split('\t')
                        typelabel = tem[0]
                        # Gscore_metaz_rank
                        file_list += ['/data1/LanYun/GScore/for_figs/fig5C/output/'+dataset_type+'_'+pw_num+'/'+'Fractional_RANK_Gscore/'+typelabel+'_'+path_cond+'_metaZ_FDRQ_delm0_rank_PCC0.5.txt']
                        #NEA_metaz_rank
                        file_list += ['/data1/LanYun/GScore/for_figs/fig5C/output/'+dataset_type+'_'+pw_num+'/'+'Fractional_RANK_NEA/'+typelabel+'_'+path_cond+'_metaZ_FDRQ_delm0_rank_PCC0.5.txt']
                        if prognostic_type=="fariable":
                            label_list += [ 'Gscore_fav']
                            label_list += [ 'NEA_fav']
                        elif prognostic_type=="adverse":
                            label_list += [ 'Gscore_adv']
                            label_list += [ 'NEA_adv']
                        
                        #prognostic_gene_set
                        target_list += ['/data1/yujia/handover/Gscore_paper/for_figure5/Figure5C_individual_CDF/survival_median/'+typelabel.lower()+f'_{prognostic_type}.txt']
                        target_list += ['/data1/yujia/handover/Gscore_paper/for_figure5/Figure5C_individual_CDF/survival_median/'+typelabel.lower()+f'_{prognostic_type}.txt']
                        print(typelabel)
                if prognostic_type=="fariable":
                    onlylabel += [ 'Gscore_fav']
                    onlylabel += [ 'NEA_fav']
                elif prognostic_type=="adverse":
                    onlylabel += [ 'Gscore_adv']
                    onlylabel += [ 'NEA_adv']
                file_save = f"/data1/LanYun/GScore/for_figs/fig5C/CDF_output/{path_cond}/{dataset_type}_{pw_num}_single_{path_cond}_CDF_delm0_PCC0.5_{prognostic_type}"
                ks_pic(file_list,label_list,target_list,file_save,onlylabel)

