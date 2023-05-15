#!/usr/bin/env python
#coding:utf-8

#=================================================================================================================================
# Name                  : 
# Created On            : 2023-03-28
# Author                : yuhuan
# Version               : 1.0
# Last Modified By      : yuhuan
# Last Modified On      : 
#=================================================================================================================================
__author__ = 'happieryu@163.com'
__version__ = 'v1.0'
__date__ = '2023.03.28'
__comment__ ='''Reform the result of mosthdepth. 
Specify mosdepth -b with a bed with genename or marker name in 4th column 
and specify --thresholds 1,50,100,200,500.'''

import os
import argparse
from argparse import RawTextHelpFormatter

##software
samtools = 'samtools'

def safeopen(infile, mode):
    '''safe open file or gziped file'''
    ofile = os.path.abspath(infile)
    if not os.path.exists(ofile):
        exit('%s is not exists.' %ofile)
    elif ofile.endswith('.gz'):
        import gzip
        return gzip.open(ofile, mode)
    else:
        return open(ofile, mode)

def getnanostat(statfile):
    statdic = {}
    ostatfile = open(statfile, 'r')
    linenum = 0
    for line in ostatfile:
        linenum += 1
        if linenum < 16 and (':' in line):
            uline = line.strip().split(':')
            statdic[uline[0].strip()] = uline[1].strip()
    return statdic

def init():
    '''init parameter'''
    parser = argparse.ArgumentParser(description = 'Comment: %s\n\nAuthor: %s\nVersion: %s\nDate: %s' \
            %(__comment__,__author__,__version__,__date__), formatter_class = RawTextHelpFormatter)
    parser.add_argument('-d', '--dir', metavar = 'path', required = True, \
            help = 'The path of mosdepth result.')
    parser.add_argument('-p', '--prefix',  metavar = 'File', required = True, \
            help = 'The prefix of mosdepth.')
    parser.add_argument('-b', '--bed',  metavar = 'File', required = True, \
            help = 'The bedfile of mosdepth.')
    parser.add_argument('-rawstat', '--rawstat',  metavar = 'File', required = True, \
            help = 'The nanostat file of raw fq.')
    parser.add_argument('-qcstat', '--qcstat',  metavar = 'File', required = True, \
            help = 'The nanostat file of filter fq.')
    parser.add_argument('-bam', '--bam',  metavar = 'File', required = True, \
            help = 'The bam file of the sample.')
    parser.add_argument('-v', '-V', help="Show version number and exit. Version = %s" \
            %__version__, action='version', version = __version__)
    argv = parser.parse_args()
    return argv

def main():
    argv = init()
    pathdir = argv.dir
    prefix = argv.prefix
    bedfile = argv.bed
    rawstat = argv.rawstat
    qcstat = argv.qcstat
    bamfile = argv.bam
    reformout = open(prefix+'.mosdepth_reform.xls', 'w')
    reformdp = open(prefix+'.regions.bed.normlize.xls', 'w')
    outcontent = ''
    geneconten = ''
    sumfile = pathdir+'/'+prefix+'.mosdepth.summary.txt'
    regiondepth = pathdir+'/'+prefix+'.regions.bed.gz'
    regionthresholds = pathdir+'/'+prefix+'.thresholds.bed.gz'
    #if os.path.exists(regionthresholds):

    # nanostat item:
    #Mean read length:                1,015.7
    #Mean read quality:                  11.4
    #Median read length:                828.0
    #Median read quality:                11.4
    #Number of reads:               372,849.0
    #Read length N50:                 1,156.0
    #STDEV read length:                 949.3
    #Total bases:               378,686,237.0
    #Number, percentage and megabases of reads above quality cutoffs
    #>Q5:    372849 (100.0%) 378.7Mb
    #>Q7:    372849 (100.0%) 378.7Mb
    #>Q10:   245613 (65.9%) 237.2Mb
    #>Q12:   152165 (40.8%) 142.7Mb
    #>Q15:   30119 (8.1%) 21.2Mb

    # qc stat information
    rawstatdic = getnanostat(rawstat)
    qcstatdic = getnanostat(qcstat)

    outcontent += 'Sample'+'\t'+prefix+'\n'
    outcontent += '\t'.join(['Raw bases(Mb)', str(round(float(rawstatdic['Total bases'].replace(',',''))/1000000.0, 2))])+'\n'
    outcontent += '\t'.join(['Raw reads(M)', str(round(float(rawstatdic['Number of reads'].replace(',',''))/1000000.0, 2))])+'\n'
    outcontent += '\t'.join(['Raw mean read length', rawstatdic['Mean read length']])+'\n'
    outcontent += '\t'.join(['Raw read length N50', rawstatdic['Read length N50']])+'\n'
    outcontent += '\t'.join(['Raw mean read quality', rawstatdic['Mean read quality']])+'\n'
    outcontent += '\t'.join(['Raw read quality>Q7', rawstatdic['>Q7']])+'\n'
    outcontent += '\t'.join(['Raw read quality>Q10', rawstatdic['>Q10']])+'\n'
    #filter
    outcontent += '\t'.join(['QC bases(Mb)', str(round(float(qcstatdic['Total bases'].replace(',',''))/1000000, 2))])+'\n'
    outcontent += '\t'.join(['QC reads(M)', str(round(float(qcstatdic['Number of reads'].replace(',',''))/1000000, 2))])+'\n'
    outcontent += '\t'.join(['QC mean read length', qcstatdic['Mean read length']])+'\n'
    outcontent += '\t'.join(['QC read length N50', qcstatdic['Read length N50']])+'\n'
    outcontent += '\t'.join(['QC mean read quality', qcstatdic['Mean read quality']])+'\n'
    outcontent += '\t'.join(['QC read quality>Q7', qcstatdic['>Q7']])+'\n'
    outcontent += '\t'.join(['QC read quality>Q10', qcstatdic['>Q10']])+'\n'
    outcontent += 'QC pass rate of bases(%)' + '\t' + str(round(float(qcstatdic['Total bases'].replace(',',''))/float(rawstatdic['Total bases'].replace(',',''))*100, 2))+'\n'
    outcontent += 'QC pass rate of reads(%)' + '\t' + str(round(float(qcstatdic['Number of reads'].replace(',',''))/float(rawstatdic['Number of reads'].replace(',',''))*100, 2))+'\n'

    total_reads = float(qcstatdic['Number of reads'].replace(',',''))
    umapped_reads = int(os.popen('%s view -c -f 0x0004 %s' %(samtools, bamfile)).read().strip())
    map_rate = round((total_reads-umapped_reads)/total_reads*100, 2)
    outcontent += 'Mapping rate of reads(%)' + '\t' + str(map_rate)+ '\n'

    genelist = []
    tmbed = open(bedfile, 'r')
    for line in tmbed:
        uline = line.strip().split('\t')
        genename = uline[3]
        if genename not in genelist:
            genelist.append(genename)

    osumfile = safeopen(sumfile, 'r')
    for line in osumfile:
        uline = line.strip().split('\t')
        if uline[0] == 'total':
            totalbase = int(uline[2])
        if uline[0] == 'total_region':
            regionlen = int(uline[1])
            regionbase = int(uline[2])
            ave_regiondepth = float(uline[3])
            base_capture = round(float(regionbase/totalbase)*100,2)
    osumfile.close()
    
    outcontent += 'Target region length (bp)'+'\t'+str(regionlen)+'\n'
    outcontent += 'Capture rate of bases (%)'+'\t'+str(base_capture)+'\n'
    outcontent += 'Average depth on target (x)'+'\t'+str(ave_regiondepth)+'\n'

    genedpdic = {}
    oregiondepth = safeopen(regiondepth, 'r')
    for line in oregiondepth:
        uline = line.decode().strip().split('\t')
        genename = uline[3]
        perlinelen = int(uline[2]) - int(uline[1])
        perlinedp = float(uline[4])
        normdp = round(perlinedp/ave_regiondepth, 4)
        reformdp.write(line.decode().strip()+'\t'+str(normdp)+'\n')
        if genename not in genedpdic:
            genedpdic.setdefault(genename,{})['len'] = 0
            genedpdic.setdefault(genename,{})['linecount'] = 0
            genedpdic.setdefault(genename,{})['linedp'] = 0
            genedpdic.setdefault(genename,{})['1x'] = 0
            genedpdic.setdefault(genename,{})['50x'] = 0
            genedpdic.setdefault(genename,{})['100x'] = 0
            genedpdic.setdefault(genename,{})['200x'] = 0
            genedpdic.setdefault(genename,{})['500x'] = 0
        genedpdic[genename]['len'] += perlinelen
        genedpdic[genename]['linecount'] += 1
        genedpdic[genename]['linedp'] += perlinedp

    oregionthresholds = safeopen(regionthresholds, 'r')
    for line in oregionthresholds:
        uline = line.decode().strip().split('\t')
        if line.decode().startswith('#'):continue
        genename = uline[3]
        perlinecov = int(uline[4])
        cov50x = int(uline[5])
        cov100x = int(uline[6])
        cov200x = int(uline[7])
        cov500x = int(uline[8])
        genedpdic[genename]['1x'] += perlinecov
        genedpdic[genename]['50x'] += cov50x
        genedpdic[genename]['100x'] += cov100x
        genedpdic[genename]['200x'] += cov200x
        genedpdic[genename]['500x'] += cov500x
    
    totalcovlen = 0
    totalcov50len = 0
    totalcov100len = 0
    totalcov200len = 0
    totalcov500len = 0
    for geneinfo in genelist:
        totalcovlen += float(genedpdic[geneinfo]['1x'])
        totalcov50len += float(genedpdic[geneinfo]['50x'])
        totalcov100len += float(genedpdic[geneinfo]['100x'])
        totalcov200len += float(genedpdic[geneinfo]['200x'])
        totalcov500len += float(genedpdic[geneinfo]['500x'])
        generegion = genedpdic[geneinfo]['len']
        geneavdp = round(float(genedpdic[geneinfo]['linedp']/genedpdic[geneinfo]['linecount']), 2)
        genecov = round(float(genedpdic[geneinfo]['1x']/genedpdic[geneinfo]['len'])*100, 2)
        genecov50 = round(float(genedpdic[geneinfo]['50x']/genedpdic[geneinfo]['len'])*100, 2)
        genecov100 = round(float(genedpdic[geneinfo]['100x']/genedpdic[geneinfo]['len'])*100, 2)
        genecov200 = round(float(genedpdic[geneinfo]['200x']/genedpdic[geneinfo]['len'])*100, 2)
        genecov500 = round(float(genedpdic[geneinfo]['500x']/genedpdic[geneinfo]['len'])*100, 2)
        geneconten += geneinfo+' length (bp)'+'\t'+str(generegion)+'\n'
        geneconten += geneinfo+' average depth (x)'+'\t'+str(geneavdp)+'\n'
        geneconten += geneinfo+' covreage ≥1x (%)'+'\t'+str(genecov)+'\n'
        geneconten += geneinfo+' covreage ≥50X (%)'+'\t'+str(genecov50)+'\n'
        geneconten += geneinfo+' covreage ≥100X (%)'+'\t'+str(genecov100)+'\n'
        geneconten += geneinfo+' covreage ≥200X (%)'+'\t'+str(genecov200)+'\n'
        geneconten += geneinfo+' covreage ≥500X (%)'+'\t'+str(genecov500)+'\n'
    totalcov = round(float(totalcovlen/regionlen)*100, 2)
    totalcov50 = round(float(totalcov50len/regionlen)*100, 2)
    totalcov100 = round(float(totalcov100len/regionlen)*100, 2)
    totalcov200 = round(float(totalcov200len/regionlen)*100, 2)
    totalcov500 = round(float(totalcov500len/regionlen)*100, 2)
    outcontent += 'Coverage ≥ 1x (%)' + '\t' + str(totalcov)+'\n'
    outcontent += 'Coverage ≥ 50x (%)' + '\t' + str(totalcov50)+'\n'
    outcontent += 'Coverage ≥ 100x (%)' + '\t' + str(totalcov100)+'\n'
    outcontent += 'Coverage ≥ 200x (%)' + '\t' + str(totalcov200)+'\n'
    outcontent += 'Coverage ≥ 500x (%)' + '\t' + str(totalcov500)+'\n'

    reformout.write(outcontent)
    reformout.write(geneconten)

if __name__ == '__main__':
    main()