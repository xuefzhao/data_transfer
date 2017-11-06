#!python
#script to extract inheritance rate / denovo rate from vcf file

import os
def vcf_readin(input_vcf):
        header=[]
        info=[]
        fin=open(input_vcf)
        for line in fin:
                pin=line.strip().split()
                if pin[0][:2]=='##':                    header.append(pin)
                else:
                        if pin[0][0]=='#':              info.append(pin)
                        else:
                                if pin[6]=='PASS':      info.append(pin)
        fin.close()
        return [header, info]

def denovo_extract(gt):
        out=[0,0,0]
        if gt[0]=='CN':
                if (int(gt[3])-int(gt[1]))*(int(gt[3])-int(gt[2]))>0:   out[0]+=1
                if (int(gt[6])-int(gt[4]))*(int(gt[6])-int(gt[5]))>0:   out[1]+=1
                if (int(gt[9])-int(gt[7]))*(int(gt[9])-int(gt[8]))>0:   out[2]+=1
        elif gt[0]=='GT':
                for i in gt[3].split('/'):
                        if not i in gt[1].split('/')+gt[2].split('/') and out[0]==0:    out[0]+=1
                for i in gt[6].split('/'):
                        if not i in gt[4].split('/')+gt[5].split('/') and out[1]==0:    out[1]+=1
                for i in gt[9].split('/'):
                        if not i in gt[7].split('/')+gt[8].split('/') and out[2]==0:    out[2]+=1
        return out

def denovo_calcu(info):
        out=[0,0,0,0,0,0,0]
        out2=[]
        for i in info:
                if i[0][0]=='#': continue
                if not i[0][0]=='#':    out[0]+=1
                tmp=denovo_calcu_from_raw_caller(i)
                tmp2=sv_check(i[8:])
                out[1]+=tmp2[0]
                out[2]+=tmp[0]
                out[3]+=tmp2[1]
                out[4]+=tmp[1]
                out[5]+=tmp2[2]
                out[6]+=tmp[2]
                out2.append(i[:8]+[tmp2[0],tmp[0],tmp2[1],tmp[1],tmp2[2],tmp[2]])
        return [out,out2]

def denovo_calcu_from_raw_caller(pin):
        info=pin_info_cha_extract(pin,'INFO_POS')
        if not info=='':
                sample=[i.split(':')[-2] for i in info.split(',')]
                gt=['GT']
                sample_name=['HG00512', 'HG00513', 'HG00514', 'HG00731','HG00732','HG00733','NA19239','NA19238','NA19240']
                for i in sample_name:
                        if not i in sample:     gt.append('0/0')
                        else:                   gt.append('0/1')
        else:
                gt=pin[8:]
        gt_reform=[]
        for i in gt:
                if i in ['.','./.']: gt_reform.append('0/0')
                else:   gt_reform.append(i)
        denovo_info=denovo_extract(gt_reform)                
        return denovo_info

def info_col_reform(pin):
        info_new=[pin[0],pin[1],pin_info_cha_extract(pin, 'END'), pin[6],pin_info_cha_extract(pin, 'SVLEN'),pin_info_cha_extract(pin, 'SVTYPE'),pin_info_cha_extract(pin, 'SOURCES')]
        return info_new

def info_matrix_reform(info):
        info_new=[]
        for i in info:
                info_new.append(info_col_reform(i)+i[8:])
        return info_new

def pin_info_cha_extract(pin, cha):
        out=''
        for i in pin[7].split(';'):
                if i.split('=')[0]==cha:
                        out=i.split('=')[1]
        return out

def sv_check(gt):
        out=[0,0,0]
        if gt[0]=='CN':
                if not int(gt[1])==int(gt[2])== int(gt[3]) ==2:  out[0]+=1
                if not int(gt[4])==int(gt[5])== int(gt[6]) ==2:  out[1]+=1
                if not int(gt[7])==int(gt[8])== int(gt[9]) ==2:  out[2]+=1
        elif gt[0]=='GT':
                if not gt[1]==gt[2]== gt[3] =='0/0': out[0]+=1
                if not gt[4]==gt[5]== gt[6] =='0/0': out[1]+=1
                if not gt[7]==gt[8]== gt[9] =='0/0': out[2]+=1
        return out

def write_denovo_matrics(inheritance_matrix,info,fileout):
        fo=open(fileout,'w')
        print('\t'.join(info[0][:8]+['chs','chs_denovo','pur','pur_denovo','yri','yri_denovo']), file=fo)
        for i in inheritance_matrix:
                if not i[0][0]=='#':
                        print('\t'.join([str(x) for x in i]), file=fo)
        fo.close()

def write_new_matrics(info_reform, fileout):
        fo=open(fileout,'w')
        print('\t'.join(['#CHR','POS','END','FILTER','SVLEN','SVTYPE','SOURCES','chs','chs_denovo','pur','pur_denovo',' yri','yri_denovo']), file=fo)
        for i in info_reform:
                print('\t'.join([str(j) for j in i]), file=fo)
        fo.close()

def main():
        input_vcf='HGSV_ILL_SV_Integration.sorted.vcf'
        [header,info]=vcf_readin(input_vcf)
        [denovo_info,inheritance_matrix]=denovo_calcu(info) #[#allSV, #CHSallSV, #CHSdenovo, #PURallSV, #PURdenovo, #YRIallSV, #YRIdenovo]
        write_denovo_matrics(inheritance_matrix,info,'HGSV_ILL_SV.inheritance.stat')
        info_reform=info_matrix_reform(inheritance_matrix)
        write_new_matrics(info_reform, 'HGSV_ILL_SV.inheritance.info')

main()

