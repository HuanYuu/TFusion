#!/usr/bin/env python
#coding:utf-8

__author__ = 'happieryu@163.com'
__version__ = 'v1.2'
__date__ = '2023.02.03'
__comment__ = '''Filter and reform the result of LongGF. Filter parameters:
1. Fusion_pair dose not satisfy 5'-3' join rule.
2. Captured gene's breakpoint(20bp) average depth < 100.
3. Reads number supported fusion < 5.
4. Fusion rate(reads_support/captured_gene_breakpoint_20bp_depth) < 0.5%.
5. None of the fuion genes in captured gene list.
6. When any of the genes in the fusion_pairs has more than one strand information ,mark as "gene_multi_stand".
7. When any of the genes in the fusion_pairs dose not locate on chr1-22,X,Y, mark as "gene_not_chr_gene".
'''

import os
import argparse
from argparse import RawTextHelpFormatter
##software

samtools='samtools'

def safeopen(infile, mode):
    '''safe open file or gziped file'''
    ofile = os.path.abspath(infile)
    if not os.path.exists(ofile):
        exit('%s is not exists.' %ofile)
    elif ofile.endswith('.gz'):
        import gzip
        return gzip.open(ofile,mode)
    else:
        return open(ofile,mode)

def get_file_list(infile):
    '''get list form file'''
    filelist = []
    ofile = os.path.abspath(infile)
    if not os.path.exists(ofile):
        exit('%s is not exists.' %ofile)
    else:
        for line in safeopen(infile, 'r'):
            if line.strip() not in filelist:
                filelist.append(line.strip())
    return filelist

def getstrand(strand_file):
    '''get gene strand information'''
    strand_dic = {}
    open_sf = safeopen(strand_file, 'r')
    mulit_strand = []
    for line in open_sf:
        uline = line.strip().split('\t')
        if uline[0] not in strand_dic:
            strand_dic[uline[0]] = uline[1]
        else:
            ##gene in strand_dic, and strand not the same. (if there are repeat rows)
            if strand_dic[uline[0]] != uline[1]:
                mulit_strand.append(uline[0])
    return strand_dic, mulit_strand

def getdepth(bamfile, region, skippara):
    bamname = os.path.basename(bamfile) 
    if not os.path.exists('tmp_depth_info'):
        os.system('mkdir tmp_depth_info')
    if skippara == 'N':
        os.system('%s depth -aa -d 100000 -r %s %s > tmp_depth_info/%s.%s.depth.txt' %(samtools, region, bamfile, bamname, region))
    tempdepth = open('tmp_depth_info/'+bamname+'.'+region+'.depth.txt', 'r')
    depthlist = []
    for line in tempdepth:
        depthlist.append(int(line.strip().split('\t')[2]))
    avedepth = round(float(sum(depthlist))/len(depthlist),2)
    return avedepth

#def getstream(gene1, gene2, gene1_strand, gene2_strand, gene1_breakp, gene2_breakp, gene1_lrange, gene2_lrange, gene1_rrange, gene2_rrange):
def getstream(strand, breakpinfo, lrange, rrange):
    breakp = int(breakpinfo.split(':')[1])
    mean_lrange = sum(lrange)/len(lrange)
    mean_rrange = sum(rrange)/len(rrange)
    if strand == '+':
        if abs(mean_lrange-float(breakp)) < abs(mean_rrange-float(breakp)):
            'upstream breakpoint'
            return "downstream_gene"
        else:
            'downstream breakpoint'
            return "upstream_gene"
    if strand == '-':
        if abs(mean_lrange-float(breakp)) < abs(mean_rrange-float(breakp)):
            return "upstream_gene"
        else:
            return "downstream_gene"

'''
GF      ARL15:EPOR 2 2 supporting reads=2/2 54310546:1;54310549:1 11378593:2 chr5:54310547 0/2:8:2 chr19:11378593 0/2:8:502
        54310546(+chr5:54310430-54310546/805196bc-45df-4477-9727-2f21ed49f53a:442-555)1 11378593(-chr19:11378468-11378593/560-684)1
        54310549(-chr5:54310420-54310549/9b8715f8-ba26-4471-86ab-69960e9fbdde:188-311)1 11378593(+chr19:11378468-11378593/66-186)1
SumGF   ARL15:EPOR 2 chr5:54310547 chr19:11378593
'''
def get_fusion(log_file, gene_strand_dic, mulit_strand, skippara, bamfile='', capgenelist=''):
    '''
    fusion_info:
    0: gene_pairs
    1: support_reads
    2: 3'breakpoint gene(5'gene)
    3: 5'breakpoint gene(3'gene)
    4: 5' gene breakpoint
    5: 3' gene breakpoint
    6: 3' gene_break_depth
    7: fusion_rate
    8: filter_tag
    '''
    fusion_info = []
    open_log = safeopen(log_file, 'r')
    needed_info = False
    for line in open_log:
        if line.startswith('GF'):
            ## get information when meet 'GF' tag
            needed_info = True
            fusion_info_temp = []
            fgene_lrange =[]
            fgene_rrange = []
            sgene_lrange = []
            sgene_rrange = []
            upstreamgene = []
            upstreamgene_strand = []
            upstreambp = []
            downstreamgene = []
            downstreamgene_strand = []
            downstreambp = []
            uline = line.strip().split()
            gene_pairs = uline[1]
            support_reads_num = uline[2]
            fgene = gene_pairs.strip().split(':')[0]
            sgene = gene_pairs.strip().split(':')[1]
            fgene_breakpoint = uline[-4]
            sgene_breakpoint = uline[-2]
            continue
        if needed_info:
            ##get information after 'GF' tag, and do not deal with 'SumGF' line.
            conditionflag = False
            capgenenum = 0
            filtertag = []
            if not line.startswith('SumGF'):
                uline = line.strip().split()
                #try:
                fgene_lrange.append(int(uline[0].split(':')[1].split('-')[0]))
                fgene_rrange.append(int(uline[0].split(':')[1].split('/')[0].split('-')[1]))
                sgene_lrange.append(int(uline[1].split(':')[1].split('-')[0]))
                sgene_rrange.append(int(uline[1].split(':')[1].split('/')[0].split('-')[1]))
                #except IndexError:
                #    print(line)
            else:
                ## get information afer 'GF' tag, when meet 'SumGF' line, summarize all needed information.
                #print(gene_pairs)
                if fgene in mulit_strand:
                    filtertag.append(fgene+'_multi_stand')
                if sgene in mulit_strand:
                    filtertag.append(sgene+'_multi_stand')
                if fgene not in gene_strand_dic:
                    filtertag.append(fgene+'_not_chr_gene')
                if  sgene not in gene_strand_dic:
                    filtertag.append(sgene+'_not_chr_gene')
                if filtertag != []:
                    fusion_info_temp.extend([gene_pairs, fgene, sgene, fgene_breakpoint, sgene_breakpoint, support_reads_num, '/', '/', '/', ';'.join(filtertag)])
                else:
                    if getstream(gene_strand_dic[fgene], fgene_breakpoint, fgene_lrange, fgene_rrange) == 'upstream_gene':
                        upstreamgene.append(fgene)
                        upstreamgene_strand.append(fgene+'('+gene_strand_dic[fgene]+')')
                        upstreambp.append(fgene_breakpoint)
                    if getstream(gene_strand_dic[fgene], fgene_breakpoint, fgene_lrange, fgene_rrange) == 'downstream_gene':
                        downstreamgene.append(fgene)
                        downstreamgene_strand.append(fgene+'('+gene_strand_dic[fgene]+')')
                        downstreambp.append(fgene_breakpoint)
                    if getstream(gene_strand_dic[sgene], sgene_breakpoint, sgene_lrange, sgene_rrange) == 'upstream_gene':
                        upstreamgene.append(sgene)
                        upstreamgene_strand.append(sgene+'('+gene_strand_dic[sgene]+')')
                        upstreambp.append(sgene_breakpoint)
                    if getstream(gene_strand_dic[sgene], sgene_breakpoint, sgene_lrange, sgene_rrange) == 'downstream_gene':
                        downstreamgene.append(sgene)
                        downstreamgene_strand.append(sgene+'('+gene_strand_dic[sgene]+')')
                        downstreambp.append(sgene_breakpoint)
                    downavedp = '/'
                    upavedp = '/'
                    furate = '/'
                    tgenedepth = 100
                    if bamfile != '':
                        if len(downstreambp) == 1:   ##not multi_gene, get 3' depth information
                            if gene_strand_dic[downstreamgene[0]] == '+':
                                downregion = downstreambp[0].split(':')[0]+':'+str(int(downstreambp[0].split(':')[1])+1)+'-'+str(int(downstreambp[0].split(':')[1])+20)
                            if gene_strand_dic[downstreamgene[0]] == '-':
                                downregion = downstreambp[0].split(':')[0]+':'+str(int(downstreambp[0].split(':')[1])-20)+'-'+str(int(downstreambp[0].split(':')[1])-1)
                            downavedp = getdepth(bamfile, downregion, skippara)
                            if capgenelist != '':
                                if downstreamgene[0] in capgenelist:
                                    capgenenum += 1
                                    if downavedp < tgenedepth:
                                        conditionflag = True
                                        #print('downavedp filter, capture_gene_downstream_depth<100')
                        else:
                            conditionflag = True
                            #print('down_multi filter, downstream_multi_gene')
                        if len(upstreambp) == 1: ##not multi_gene, get 5' depth information
                            if gene_strand_dic[upstreamgene[0]] == '+':
                                upregion = upstreambp[0].split(':')[0]+':'+str(int(upstreambp[0].split(':')[1])-20)+'-'+str(int(upstreambp[0].split(':')[1])-1)
                            if gene_strand_dic[upstreamgene[0]] == '-':
                                upregion = upstreambp[0].split(':')[0]+':'+str(int(upstreambp[0].split(':')[1])+1)+'-'+str(int(upstreambp[0].split(':')[1])+20)
                            upavedp = getdepth(bamfile, upregion, skippara)
                            if capgenelist != '':
                                if upstreamgene[0] in capgenelist:
                                    capgenenum += 1
                                    if upavedp < tgenedepth:
                                        conditionflag = True
                                        #print('upavedp filter, capture_gene_upstream_depth<100')
                        else:
                            conditionflag = True
                            #print('up_multi filter, upstream_multi_gene')
                        if len(downstreambp) == 1 and len(upstreambp) == 1:
                            if downstreamgene[0] in capgenelist:
                                totaldepth = downavedp
                                supportdepth = upavedp
                            elif upstreamgene[0] in capgenelist:
                                totaldepth = upavedp
                                supportdepth = downavedp
                            else:
                                totaldepth = max(downavedp, upavedp)
                                supportdepth = min(downavedp, upavedp)
                            theorate = round(float(supportdepth)/totaldepth*100, 3)
                            furate = round(float(support_reads_num)/totaldepth*100, 3)
                            if furate < 0.5:  ## 0.5%, for any gene in capgenelist
                                conditionflag = True
                            if theorate < 1: ## 1%, theoretically support fusion rate (if all supportdepth reads are fusions, the detect limit is 1%)
                                conditionflag = True
                            if float(supportdepth) < float(support_reads_num)+1: # not confident, theoretically support reads number must larger than support_reads_num.
                                conditionflag = True
                    if int(support_reads_num) < 5:
                        conditionflag = True
                        #print('support num filter, support_reads_number<5')
                    if capgenenum == 0:
                        conditionflag = True
                        #print('cap_gene_num filter, no fusion genes in capture gene list.')
                    if  conditionflag:
                        filtertag.append('Filter')
                    fusion_info_temp = [','.join(upstreamgene)+'-'+','.join(downstreamgene), \
                                        ','.join(upstreamgene_strand), \
                                        ','.join(downstreamgene_strand), \
                                        ','.join(upstreambp), \
                                        ','.join(downstreambp), \
                                        str(support_reads_num), \
                                        str(upavedp), \
                                        str(downavedp), \
                                        str(furate), \
                                        ';'.join(filtertag)]
            if fusion_info_temp != []:
                fusion_info.append(fusion_info_temp)
    return fusion_info

def init():
    '''init parameter'''
    parser = argparse.ArgumentParser(description = 'Comment: %s\n\nAuthor: %s\nVersion: %s\nDate: %s' \
            %(__comment__,__author__,__version__,__date__), formatter_class = RawTextHelpFormatter)
    parser.add_argument('-i', '--inf', metavar = 'File', required = True, \
            help = 'The log_file generated by LongGF. Required.')
    parser.add_argument('-o', '--out',  metavar = 'File', required = False, \
            help = 'The out file. If not provide, add \"out\" to the input file\'s name.')
    parser.add_argument('-sf', metavar = 'File', required = True, \
            help = 'The "gene\tstand_information" file. Required')
    parser.add_argument('-bam', metavar = 'File', required = False, \
            help = 'Bam file of the sample. Influence FP mark if not provide.')
    parser.add_argument('-gene', metavar = 'File', required = False, \
            help = 'Capture or target gene list, one gene per line. Influence FP mark if not provide.')
    parser.add_argument('-skipc', choices=['Y', 'N'], default='N', required = False, \
            help = 'Skip samtools depth with Y when depth calculated before. Only use in test.')
    parser.add_argument('-v', '-V', help="Show version number and exit. Version = %s" \
            %__version__, action='version', version = __version__)
    argv = parser.parse_args()
    return argv

def main():
    argv = init()
    inputf = argv.inf
    outf = argv.out
    strandf = argv.sf
    skipset = argv.skipc
    bamf = ''
    genelist = ''
    if argv.bam:
        bamf = argv.bam
    if argv.gene:
        genelist = get_file_list(argv.gene)
        print("genelist:\t",genelist)

    if not outf:
        out = open(os.path.splitext(os.path.basename(os.path.abspath(inputf)))[0]+'.out.xls', 'w')
    else:
        out = open(outf, 'w')

    headinfo = [os.path.splitext(os.path.basename(os.path.abspath(inputf)))[0]+'_Gene_pairs', \
                'Upstream_gene', \
                'Downstream_gene', \
                'Upstream_gene_breakpoint', \
                'Downstream_gene_breakpoint', \
                'Support_reads_number', \
                'Upstream_gene_20bp_breakpoint_depth', \
                'Downstream_gene_20bp_breakpoint_depth', \
                'Fusion_rate(%)', \
                'Multiple_strand_filter']
    stranddic, multidic = getstrand(strandf)
    fusioncontent = get_fusion(inputf, stranddic, multidic, skipset, bamf, genelist)
    out.write('\t'.join(headinfo)+'\n')
    for item in fusioncontent:
        out.write('\t'.join(item)+'\n')
    out.close()
    
if __name__ == '__main__':
    main()