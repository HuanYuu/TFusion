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
__comment__ = '''QC, mapping, fusion, stats'''

import os
import sys
import argparse
from argparse import RawTextHelpFormatter

##scriptpath
scriptpath = os.path.dirname(os.path.abspath(sys.argv[0]))

###Software
##********************** User can modify all these softwares' and databse files' path *************************
##python
python3 = 'python'
#remove adapter and debarcode
#user can modify porechop-runner.py path based on their install
porechop = scriptpath+'/Porechop*/porechop-runner.py'
##Filt
nanofilt = 'NanoFilt'
##Mappng
minimap2 = 'minimap2'
samtools = 'samtools'
sambamba = 'sambamba'
##statistic and plot
nanostat = 'NanoStat'
#nanoplot = 'NanoPlot'
mosdepth = 'mosdepth'
statreform = scriptpath+'/src/mosdepth_bedstat_reform.py'
##fusion 
longgf = 'LongGF'
reformlgf = scriptpath+'/src/reform_LongGF_Result.py'
###databse
humangref = scriptpath+'/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
humangreffai = scriptpath+'/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
humangtf = scriptpath+'/reference/Homo_sapiens.GRCh38.109.gtf'
genestrand = scriptpath+'/reference/Homo_sapiens.GRCh38.109.gene.strand.txt'


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
    
def cpath(pathdir, namestr):
    pathid = pathdir+'/'+namestr
    if not os.path.exists(pathid):
        os.system('mkdir -p %s' %pathid)
    return pathid

def tdic(file):
    filedic = {}
    tmpfile = safeopen(file, 'r')
    for line in tmpfile:
        uline = line.strip().split('\t')
        filedic.setdefault(uline)

def saminfo(samlist, outdir):
    '''
    samdic[sampleid][sampleid||libid||barcodeid] = [rawfq1, rawfq2]
    multidic[poolibid][fq]: [[sampleid, barcodeid], [sampleid, barcodeid]]
    '''
    osamlist = safeopen(samlist, 'r')
    newsamlist = open(outdir+'/sample.list.reconfig', 'w')
    samdic = {}
    samplelist = []
    pooldic = {}
    barcodelist = []
    for line in osamlist:
        if line.startswith('#'):continue
        uline = line.strip().split()
        sampleid = uline[0]
        libid = uline[1]
        barcodeid = uline[2]
        fqdata = uline[3]
        if not os.path.isabs(fqdata):
            fqdata = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(samlist)), fqdata))
        prefix = '_'.join([sampleid, libid, barcodeid])
        prefixtmp = '||'.join([sampleid, libid, barcodeid])
        if sampleid not in samplelist:
            samplelist.append(sampleid)
        ##check if there is pooling barcode id
        if barcodeid == '' or barcodeid=='-' or barcodeid=='/':
            barcodeid = 'none'
            newline = '\t'.join([sampleid, libid, barcodeid, fqdata])+'\n'
            prefix = '_'.join([sampleid, libid, barcodeid])
            bcfqdata = fqdata
        else:
            bcfqdata = outdir+'/Rawdata/debarcode/'+libid+'/'+prefix+'.fastq'
            newline = '\t'.join([sampleid, libid, barcodeid, bcfqdata])+'\n'
        ##write debarcode fq path
        newsamlist.write(newline)
        ##barcode list for judge if need do debarcode, construct barcode dic info.
        if barcodeid not in barcodelist:
            barcodelist.append(barcodeid)
        ##for pool samples, one fq can links to many sample, pooldic save rawfq.
        if libid not in pooldic:
            pooldic.setdefault(libid, {})[fqdata] = [[sampleid, barcodeid]]
        if fqdata not in pooldic[libid]:
            pooldic[libid][fqdata] = [[sampleid, barcodeid]]
        if [sampleid, barcodeid] not in pooldic[libid][fqdata]:
            pooldic[libid][fqdata].append([sampleid, barcodeid])
        #construct sample dic, one prefixtmp may links to many fastq(many cell).
        #samdic save bcfqdata
        if sampleid not in samdic:
            samdic.setdefault(sampleid,{})[prefixtmp] = [bcfqdata]
        if prefixtmp not in samdic[sampleid]:
            samdic[sampleid][prefixtmp] = [bcfqdata]
        if bcfqdata not in samdic[sampleid][prefixtmp]:
            samdic[sampleid][prefixtmp].append(bcfqdata)
    return samdic, samplelist, pooldic, barcodelist

def getnanostat(statfile):
    statdic = {}
    statitem = []
    ostatfile = open(statfile, 'r')
    linenum = 0
    for line in ostatfile:
        linenum += 1
        if linenum < 16 and (':' in line):
            uline = line.strip().split(':')
            statitem.append(uline[0].strip())
            statdic[uline[0].strip()] = uline[1].strip()
    return statitem, statdic

def init():
    '''init parameter'''
    parser = argparse.ArgumentParser(description = 'Comment: %s\n\nAuthor: %s\nVersion: %s\nDate: %s' \
            %(__comment__,__author__,__version__,__date__), formatter_class = RawTextHelpFormatter)
    parser.add_argument('-s', '--sample_list', metavar = 'file', required = True, \
            help = 'The sample list strore : sampleid, libraryid, barcodeid, fq_data_path.')
    parser.add_argument('-o', '--outdir', metavar = 'str', required = True, \
            help = 'The path to save result.')
    parser.add_argument('-b', '--bed', metavar = 'file', required = True, \
            help = 'The bed file to statistic, with 4th coulum tag, eg gene/region_mark.')
    parser.add_argument('-t', '--thread', metavar = 'int', required = False, default=8, \
            help = 'The thread to run mapping.')
    parser.add_argument('-rna', '--rna', default=False, action = 'store_true', \
            help = 'RNA sequence data input.')
    parser.add_argument('-lgf_set', '--longgf_set', metavar = 'str', required = False, default='90,20,90,0,0,5', \
            help = 'min-overlap-len, bin_size, min-map-len, pseudogene, Secondary_alignment, min_sup_read for LongGF, \
            default=90,20,90,0,0,5.')
    parser.add_argument('-gene', '--target_gene', metavar = 'file', required = True, \
            help = 'The gene list to detect fusions, one gene per line.')
    parser.add_argument('-v', '-V', help="Show version number and exit. Version = %s" \
            %__version__, action='version', version = __version__)
    argv = parser.parse_args()
    return argv

def main():
    argv = init()
    slist = os.path.abspath(argv.sample_list)
    outdir = os.path.abspath(argv.outdir)
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' %outdir)
    bedfile = os.path.abspath(argv.bed)
    genefile = os.path.abspath(argv.target_gene)
    threadnum = argv.thread
    lgf_set = ' '.join(argv.longgf_set.strip().split(','))

    rawpath = cpath(outdir, 'Rawdata')
    qcpath = cpath(outdir, 'QC')
    mappath = cpath(outdir, 'Mapping')
    fusionpath = cpath(outdir, 'Fusion')
    reportpath = cpath(outdir, 'Report')
    shellpath = cpath(outdir, 'Shell')
    shell_list = []

    '''multidic[poolibid][fq]: [[sampleid, barcodeid], [sampleid, barcodeid]]'''
    samdic, samlist, pooldic, barcodeli = saminfo(slist, outdir)
    #print(samdic)
    #print(pooldic)
    if len(barcodeli) > 1:
        barcodedir = cpath(rawpath, 'debarcode')
        write_debarcode_shf = shellpath+'/'+'step0_debarcode.sh'
        write_debarcode_sh = open(write_debarcode_shf, 'w')
        for poolid in pooldic:
            merge_bcfq_cmd = ''
            mergebcfqtmp = poolid+'_merge_pool.fastq'
            mergebcfq = poolid+'_merge_pool.fastq.gz'
            if len(pooldic[poolid]) > 1:
                for poolfqtmp in pooldic[poolid]:
                    if poolfqtmp.endswith('.gz'):
                        merge_bcfq_cmd += '\nzcat %s >> %s' %(poolfqtmp, mergebcfqtmp)
                    else:
                        merge_bcfq_cmd += '\ncat %s >> %s' %(poolfqtmp, mergebcfqtmp)
            else:
                if list(pooldic[poolid].keys())[0].endswith('.gz'):
                    merge_bcfq_cmd += '\nln -sf %s %s\n' %(list(pooldic[poolid].keys())[0], mergebcfq)
                    merge_bcfq_cmd += '\ngzip -dc %s > %s \n' %(mergebcfq, mergebcfqtmp)
                else:
                    merge_bcfq_cmd += '\nln -sf %s %s\n' %(list(pooldic[poolid].keys())[0], mergebcfqtmp)
            for fqdata in pooldic[poolid]:
                if len(pooldic[poolid][fqdata]) > 1: ##with multi_barcode, need do de_barcode
                    ##merge multi_fq to one
                    debarcode_path = cpath(barcodedir, poolid)
                    barcodeshf = debarcode_path+'/'+poolid+'_debarcode.sh'
                    barcodesh = open(barcodeshf, 'w')
                    debarcode_cmd = 'set -e \necho %s de_barcode start: `date`\n' %poolid
                    debarcode_cmd += 'cd %s\n' %debarcode_path
                    debarcode_cmd += merge_bcfq_cmd
                    debarcode_cmd += '\n\n'+python3+' '+porechop+ \
                                    ' \\\n    -i %s' %mergebcfqtmp+ \
                                    ' \\\n    --format fastq'+ \
                                    ' \\\n    -b %s\n' %debarcode_path
                    #debarcode_cmd += '\n\n'+python3+' '+porechop+' \\\n'.join(['-i %s ', '--format fastq ', '-b %s '] %(mergebcfqtmp, debarcode_path))
                    for item in pooldic[poolid][fqdata]:##move barcode id to sampleid
                        prefixid = '_'.join([item[0], poolid, item[1]]) #sample, libid, barcodeid
                        debarcode_cmd += '\nmv %s %s \n' %(item[1]+'.fastq', prefixid+'.fastq')
                    debarcode_cmd += '\necho %s debarcode end: `date` \n' %poolid
                    barcodesh.write(debarcode_cmd)
        ##****** poolid_debarcode.sh ******
        shell_list.append(barcodeshf)
        write_debarcode_sh.write('set -e \n##debarcode\n')
        for bashell in shell_list:
            write_debarcode_sh.write('\nsh %s \n' %bashell)
        write_debarcode_sh.close()
    
    '''samdic[sampleid][sampleid||libid||barcodeid] = [rawfq1, rawfq2]'''

    for samid in samlist:
        samshell_list = []
        # merge fq, samid_mergefq.sh
        #reprefix = '_'.join(samdic[samid].split('||'))
        rawsampath = cpath(rawpath, samid)
        mergefqtmp = samid+'.merge.fastq'
        mergefq = samid+'.merge.fastq.gz'
        mergeshf = rawsampath+'/'+samid+'_mergefq.sh'
        mergesh = open(mergeshf, 'w')
        mergecmd = 'set -e \necho Raw_merge_fastq start: `date` \n'
        mergecmd += 'cd %s\n' %rawsampath
        if len(samdic[samid]) > 1: ##more than one sam_lib_bc id
            for pretmp in samdic[samid]:
                for fqtmp in samdic[samid][pretmp]:
                    if fqtmp.endswith('.gz'):
                        mergecmd += 'zcat %s >> %s\n' %(fqtmp, mergefqtmp)
                    else:
                        mergecmd += 'cat %s >> %s\n' %(fqtmp, mergefqtmp)
            mergecmd += 'gzip -f %s \n' %(mergefqtmp)
        else: ## only one sam_lib_bc id
            for pretmp in samdic[samid]:
                if len(samdic[samid][pretmp]) == 1: ##only one fq for the lib
                    if samdic[samid][pretmp][0].endswith('.gz'):
                        mergecmd += 'cp %s %s\n' %(samdic[samid][pretmp][0], mergefq)
                    else:
                        mergecmd += 'cp %s %s\n' %(samdic[samid][pretmp][0], mergefqtmp)
                        mergecmd += 'gzip -f %s \n' %(mergefqtmp)
                else: ## more than one fq for same lib
                    for fqtmp in samdic[samid][pretmp]:
                        if fqtmp.endswith('.gz'):
                            mergecmd += 'zcat %s >> %s\n' %(fqtmp, mergefqtmp)
                        else:
                            mergecmd += 'cat %s >> %s\n' %(fqtmp, mergefqtmp)
                    mergecmd += 'gzip -f %s \n' %(mergefqtmp)
        mergecmd += 'md5sum %s > %s.md5.txt \n' %(mergefq, mergefq)
        mergecmd += 'echo Raw_merge_fastq end: `date` \n'
        mergesh.write(mergecmd)
        mergesh.close()
        #******write mergeshf ******
        samshell_list.append(mergeshf)
        # raw stat
        #file : rawpath/samid/samid.rawfq.nanostat.xls
        rawstatshf = rawsampath+'/'+samid+'_raw.stat.sh'
        rawstatsh = open(rawstatshf, 'w')
        rawstatcmd = 'set -e \necho Raw_statistic start: `date` \n'
        rawstatcmd += 'cd %s\n\n' %rawsampath
        rawstatcmd += nanostat+' -t %s' %threadnum+ \
                    ' \\\n    --fastq %s' %mergefq+ \
                    ' \\\n    -o %s' %rawsampath+ \
                    ' \\\n    -n %s\n' %(samid+'.rawfq.nanostat.xls')
        rawstatcmd += 'echo Raw_statistic fq end: `date` \n'
        rawstatfile = rawsampath+'/'+samid+'.rawfq.nanostat.xls'
        rawstatsh.write(rawstatcmd)
        rawstatsh.close()
        #****** write rawstatshf ******
        samshell_list.append(rawstatshf+' &')
        # porochop deadapter, samid_qc.sh
        filterpath = cpath(qcpath, samid)
        absmergefq = rawsampath+'/'+mergefq
        deadapterfq = samid+'.merge.deadapter.fastq'
        filterfq = samid+'.filter.fastq.gz'
        filtershf = filterpath+'/'+samid+'_qc.sh'
        filtersh = open(filterpath+'/'+samid+'_qc.sh', 'w')
        filtercmd = 'set -e \necho QC_remove_adapter start: `date` \n'
        filtercmd += 'cd %s\n\n' %filterpath
        filtercmd += python3+' '+porechop + \
                    ' -t %s' %threadnum+ \
                    ' \\\n    -i %s' %absmergefq+ \
                    ' \\\n    -o %s' %deadapterfq+ \
                    ' \\\n    --format fastq\n'
        filtercmd += 'echo QC_remove_adapter end: `date` \n\n'
        # and nanofilt filter
        filtercmd += 'echo QC_filter start: `date` \n'
        filtercmd += nanofilt+' -q 7 -l 100 \\\n    %s | \ngzip > %s \n' %(deadapterfq, filterfq)
        filtercmd += '\nmd5sum %s > %s.md5.txt \n' %(filterfq, filterfq)
        filtercmd += '\nrm -f %s \n' %(deadapterfq)
        filtercmd += '\necho QC_filter end: `date` \n'
        filtersh.write(filtercmd)
        filtersh.close()
        #****** write filtershf ******
        samshell_list.append(filtershf)
        # nanoplot, samid_nanoplot.sh
        # statfile: qcpath/samid/samid_nanoplot.NanoStats.txt
        plotshf = filterpath+'/'+samid+'_qc_nanostat.sh'
        plotsh = open(filterpath+'/'+samid+'_qc_nanostat.sh', 'w')
        plotcmd = 'echo QC_statistic start: `date` \n'
        plotcmd += 'cd %s\n' %filterpath
        #plotcmd += nanoplot+' -t %s' %threadnum+ \
        #        ' \\\n    --fastq %s' %filterfq+ \
        #        ' \\\n    --plots dot'+ \
        #        ' \\\n    --maxlength 30000'+ \
        #        ' \\\n    --prefix %s' %(samid+'_nanoplot.')+ \
        #        ' \\\n    --info_in_report \n'
        plotcmd += nanostat+' -t %s' %threadnum+ \
                    ' \\\n    --fastq %s' %filterfq+ \
                    ' \\\n    -o %s' %filterpath+ \
                    ' \\\n    -n %s\n' %(samid+'.qc.nanostat.xls')
        plotcmd += 'echo QC_statistic end: `date` \n'
        plotsh.write(plotcmd)
        plotsh.close()
        #****** write plotshf ******
        samshell_list.append(plotshf+' &')
        #qcstatfile = filterpath+'/'+samid+'_nanoplot.NanoStats.txt'
        qcstatfile = filterpath+'/'+samid+'.qc.nanostat.xls'
        # mapping
        mapsampath = cpath(mappath, samid)
        absfilterfq = filterpath+'/'+filterfq
        mapshf = mapsampath+'/'+samid+'_mapping.sh'
        mapsh = open(mapsampath+'/'+samid+'_mapping.sh', 'w')
        mapcmd = 'set -e \necho Mapping start: `date` \n\n'
        mapcmd += 'cd %s\n\n' %mapsampath
        samrawbam = samid+'.raw.bam'
        sambam = samid+'.sorted.bam'
        if argv.rna:
            print('Data type: RNA seq data.')
            mapcmd += minimap2 + \
                ' -t %s' %threadnum+ \
                ' -ax splice'+ \
                ' \\\n    %s' %humangref+ \
                ' \\\n    %s |\n%s view -@ %s -hbS' %(absfilterfq, samtools, threadnum)+ \
                ' \\\n    -t %s' %humangreffai+ \
                ' \\\n    -o %s\n' %samrawbam
        else:
            print('Data type: gDNA seq data.')
            mapcmd += minimap2 + \
                ' -t %s -ax map-ont' %threadnum+ \
                ' \\\n    %s' %humangref+ \
                ' \\\n    %s |\n%s view -@ %s -hbS' %(absfilterfq, samtools, threadnum)+ \
                ' \\\n    -t %s' %humangreffai+ \
                ' \\\n    -o %s\n' %samrawbam
        mapcmd += '\n'+sambamba +' sort'+ \
                ' \\\n    -t %s' %threadnum+ \
                ' \\\n    -m 6G'+ \
                ' \\\n    --tmpdir %s_sambamba_sort.tmp' %samid+ \
                ' \\\n    -o %s' %sambam+ \
                ' \\\n    %s\n' %samrawbam
        mapcmd += '\n'+sambamba +' index -t %s %s' %(threadnum, sambam)
        mapcmd += '\n\nrm -f %s' %samrawbam
        mapcmd += '\n\nmd5sum %s > %s.md5.txt' %(sambam, sambam)
        mapcmd += '\n\necho Mapping end: `date` \n'
        mapsh.write(mapcmd)
        mapsh.close()
        #****** write mapshf ******
        samshell_list.append(mapshf)
        # LongGF fusion
        fusampath = cpath(fusionpath, samid)
        lgf_shf = fusampath+'/'+samid+'_fusion_lgf.sh'
        lgf_sh = open(lgf_shf, 'w')
        lgf_log = fusampath+'/'+samid+'.LongGF.log'
        lgf_result = fusampath+'/'+samid+'.LongGF.reform.xls'
        sambam_id = samid+'.readname.sorted.bam'
        lgf_cmd = 'set -e \necho LongGF fusion start: `date` \n\n'
        lgf_cmd += 'cd %s\n' %mapsampath
        lgf_cmd += samtools + ' sort -n -O BAM %s -o %s\n\n' %(sambam, sambam_id)
        lgf_cmd += 'cd %s\n' %fusampath
        lgf_cmd += longgf + \
                ' \\\n    %s' %(mapsampath+'/'+sambam_id)+ \
                ' \\\n    %s' %humangtf+ \
                ' \\\n    %s > %s\n' %(lgf_set, lgf_log)
        lgf_cmd += '\n'+python3+' '+reformlgf+ \
                ' \\\n -i %s' %lgf_log+ \
                ' \\\n -bam %s' %(mapsampath+'/'+sambam)+ \
                ' \\\n -gene %s' %genefile+ \
                ' \\\n -sf %s' %genestrand+ \
                ' \\\n -o %s\n' %lgf_result
        lgf_cmd += '\n\ncd %s' %reportpath
        lgf_cmd += '\nln -sf ../Fusion/%s/%s.LongGF.reform.xls .\n' %(samid, samid)
        lgf_cmd += '\nrm -f %s' %mapsampath+'/'+sambam_id
        lgf_cmd += '\n\necho LongGF fusion end: `date` \n'
        lgf_sh.write(lgf_cmd)
        #****** write lgf_fusion_shf ******
        samshell_list.append(lgf_shf +' &')
        # mapping stat
        depthshf = mapsampath+'/'+samid+'_mapstat.sh'
        depthsh = open(mapsampath+'/'+samid+'_mapstat.sh', 'w')
        depthcmd = 'set -e \necho Mapping_statistic start: `date` \n'
        depthcmd += 'cd %s\n\n' %mapsampath
        depthcmd += mosdepth+' -t %s' %threadnum+ \
                ' -T 1,50,100,200,500'+ \
                ' \\\n    -b %s' %bedfile+ \
                ' \\\n    %s' %samid+ \
                ' \\\n    %s\n' %sambam
        depthcmd += '\necho Mapping_statistic end: `date` \n'
        depthsh.write(depthcmd)
        depthsh.close()
        #****** write depthshf ******
        samshell_list.append(depthshf +' &')
        # reform stat
        statshf = mapsampath+'/'+samid+'_reformstat.sh'
        statsh = open(mapsampath+'/'+samid+'_reformstat.sh', 'w')
        #statreformf = mapsampath+'/'+samid+'.mosdepth_reform.xls'
        statcmd = 'set -e\necho Statistic_result_reform start: `date` \n'
        statcmd += '\ncd %s\n' %mapsampath
        statcmd += python3+' '+statreform+ \
                ' \\\n    -d %s' %mapsampath+ \
                ' \\\n    -b %s' %bedfile+ \
                ' \\\n    -rawstat %s' %rawstatfile+ \
                ' \\\n    -qcstat %s' %qcstatfile+ \
                ' \\\n    -bam %s' %sambam+ \
                ' \\\n    -p %s\n' %samid
        statcmd += '\ncd %s' %reportpath
        statcmd += '\nln -sf ../Mapping/%s/%s.mosdepth_reform.xls .\n' %(samid, samid)
        statcmd += 'echo Statistic_result_reform end: `date` \n'
        #reformout_file: prefix+'.mosdepth_reform.xls'
        #reformdp_file: prefix+'.regions.bed.normlize.xls'
        statsh.write(statcmd)
        statsh.close()
         #****** write statshf ******
        samshell_list.append(statshf)
        tmpshc = 0
        write_all_shell = open(shellpath+'/'+'run_'+samid+'_analysis.sh', 'w')
        write_all_shell_cmd = 'set -e\necho run all analysis of %s start: `date` \n' %samid
        for shellinfo in samshell_list:
            tmpshc += 1
            write_all_shell_cmd += '\n## step%s' %str(tmpshc)
            if tmpshc == len(samshell_list):
                write_all_shell_cmd += '\nwait'
            write_all_shell_cmd += '\nsh %s \n' %shellinfo
        write_all_shell_cmd += '\necho run all analysis of %s end: `date` \n' %samid
        write_all_shell.write(write_all_shell_cmd)
        write_all_shell.close()


if __name__ == '__main__':
    main()