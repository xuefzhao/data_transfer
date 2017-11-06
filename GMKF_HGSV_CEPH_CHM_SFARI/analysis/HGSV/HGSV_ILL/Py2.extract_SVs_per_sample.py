#!python
#script to extract HGSV calls for each individual
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

def pin_info_cha_extract(pin, cha):
        out=''
        for i in pin[7].split(';'):
                if i.split('=')[0]==cha:
                        out=i.split('=')[1]
        return out

def sample_all_name_extract(info):
        out=[]
        for i in info:
                if i[0][0]=='#': continue
                else:
                        out.append([i[0],i[1],pin_info_cha_extract(i,'END'),pin_info_cha_extract(i,'SVTYPE') ]+[sorted(unify_list(sample_name_extract(i)))]) 
        return out

def sample_name_extract(pin):
        sample=[]
        info_pin=pin_info_cha_extract(pin,'INFO_POS')        
        if not info_pin=='':
                sample=[i.split(':')[-2] for i in info_pin.split(',')]
        else:
                gt=pin[9:]
                sample_name= ['HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240']
                sample=[]
                for j in sample_name:
                        if not gt[sample_name.index(j)] in ['.','./.','0/0']:
                                sample.append(j)
        return sample

def unify_list(list):
        out=[]
        for i in list:
                if not i in out:        
                        out.append(i)
        return out

def write_individual_bed(info_per_samp):
        sample_name= ['HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240']
        for i in sample_name:
                fo=open(i+'.HGSV_ILL.bed','w')
                for j in info_per_samp:
                        if i in j[-1]:
                                print('\t'.join(j[:4]), file=fo)
                fo.close()

def main():
        input_vcf='HGSV_ILL_SV_Integration.sorted.vcf'
        [header,info]=vcf_readin(input_vcf)
        info_per_samp=sample_all_name_extract(info)
        write_individual_bed(info_per_samp)

import os
main()



