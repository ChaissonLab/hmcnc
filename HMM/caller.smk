import os
import tempfile
import subprocess
import os.path
import json

# Snakemake and working directories
#SD = os.path.dirname(workflow.snakefile)
#print(SD)

RD = config["scr"]


if "r" in config:
    REP=config['r']
else:
    REP=RD+"/annotation/repeatMask.merged.bed"

if "index" in config:
    ASM=config['index']
else:
    ASM=RD+"/annotation/hg38.fa.fai"



GEN=RD+"/annotation/gencode.gene.bed"

#config("hmm_caller.json")


#configfile: config["configfile"]
if config["outdir"]==".":
    outdir="hmm"
else:
    outdir=config["outdir"]+"/hmm"

fai= open(ASM)

contigs = [l.split()[0].strip().replace("|","_") for l in fai]

#fofn

bamt =  config["bam"]

bamm = bamt.split("/")[-1]
prefix_bam = os.path.splitext(bamm)[0]


ep=config["epsi"]


cov=config["coverage"]
#mask=config["repeatMask"]

sub=config["subread"]
mq=config["mq"]
ml=config["minL"]

rule all:
    input:
        bed=expand("{outdir}/cov.bed",outdir=outdir),
        covbed=expand("{outdir}/cov.no_subread.bed",outdir=outdir),
        windows=expand("{outdir}/ref_windows.bed",outdir=outdir),

        bins=expand("{outdir}/coverage.bins.bed.gz",outdir=outdir),
        avg=expand("{outdir}/mean_cov.txt",outdir=outdir),
        cov=expand("{outdir}/{ctg}.viterout.txt", ctg=contigs,outdir=outdir),
        cn=expand("{outdir}/copy_number.{ctg}.bed", ctg=contigs,outdir=outdir),
        allCN=expand("{outdir}/copy_number.tsv",outdir=outdir),
        call=expand("{outdir}/DUPcalls.copy_number.tsv",outdir=outdir),
        delcall=expand("{outdir}/DELcalls.copy_number.tsv",outdir=outdir),

        maskedcall=expand("{outdir}/DUPcalls.masked_CN.tsv",outdir=outdir),
    #    inter="{outdir}/DUPcalls.rep_int.bed",

        compCall=expand("{outdir}/DUPcalls.composite.bed",outdir=outdir),
        maskedCompcall=expand("{outdir}/DUPcalls.masked_CN.composite.tsv",outdir=outdir),
    #    interC="{outdir}/DUPcalls.rep_int.composite.bed",

        plot=expand("{outdir}/{bm}.noclip.{ep}.pdf",bm=prefix_bam,outdir=outdir,ep=ep),
        sumcall=expand("{outdir}/CallSummary.{ep}.tsv",outdir=outdir,ep=ep),
        Vsumcall=expand("{outdir}/CallSummary.verbose.{ep}.tsv",outdir=outdir,ep=ep),
        geneCount=expand("{outdir}/DUP.gene_count.bed",outdir=outdir),

rule MakeCovBed:
    input:
        bam=config["bam"],
    output:
        bed=protected("{outdir}/cov.bed"),
    params:
        rd=RD,
        out={outdir},
        mapq={mq},
    shell:"""
mkdir -p {params.out}
samtools view -q {params.mapq} -F 2304 -@ 3 {input.bam} | {params.rd}/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""

if sub==0:
    rule FilterSubreads0:
        input:
            bed="{outdir}/cov.bed",
        output:
            covbed="{outdir}/cov.no_subread.bed",
        params:
            rd=RD,
        shell:"""
    	cp {input.bed} {output.covbed}
	"""
else:
    rule FilterSubreads1:
        input:
            bed="{outdir}/cov.bed",
        output:
            covbed="{outdir}/cov.no_subread.bed",
        params:
            rd=RD,
            mapq={mq},
        shell:"""

           cat {input.bed} | python {params.rd}/filter.subread.py --mapq {params.mapq} |sort -k1,1 -k2,2n > {output.covbed}
    """

rule MakeIntersect:
    input:
        covbed="{outdir}/cov.no_subread.bed",
    output:
        windows="{outdir}/ref_windows.bed",
        bins=protected("{outdir}/coverage.bins.bed.gz"),
    params:
        rd=RD,
        asm=ASM,#config['index'],
    shell:"""
bedtools makewindows -b <( awk 'BEGIN {{OFS="\\t"}} {{print $1,0,$2}} ' {params.asm} ) -w 100 | intersectBed -v -a stdin -b {params.rd}/annotation/lc.bed {params.rd}/annotation/hg38.telomere.extended.bed {params.rd}/annotation/gap.bed {params.rd}/annotation/centromere.bed |bedtools sort > {output.windows}

intersectBed -sorted -c -a {output.windows} -b {input.covbed}| bgzip -c > {output.bins}
tabix -C {output.bins}
"""
### annotation?



if cov !="No":
    rule GetMeanCoverage0:
        input:
            bins="{outdir}/coverage.bins.bed.gz",
        output:
            avg="{outdir}/mean_cov.txt",
        params:
            rd=RD,
            cove={cov},
        shell:"""
    echo '{params.cove}' > {output.avg}
    """
else:
    rule GetMeanCoverage:
        input:
            bins="{outdir}/coverage.bins.bed.gz",
        output:
            avg="{outdir}/mean_cov.txt",
        params:
            rd=RD,
        shell:"""
    zcat {input.bins} | awk 'BEGIN{{OFS="\\t";c=0;sum=0;}} sum=sum+$4;c=c+1;END{{print sum/c;}}' | tail -1> {output.avg}
    """

rule RunVitter:
    input:
        avg="{outdir}/mean_cov.txt",
        bins="{outdir}/coverage.bins.bed.gz",
    output:
        cov=protected("{outdir}/{contig}.viterout.txt"),
    params:
        rd=RD,
        contig_prefix="{contig}",
        eps={ep},
    shell:"""
echo {params.contig_prefix}
mean=$(cat {input.avg})
touch {output.cov}
tabix {input.bins} {wildcards.contig} | cut -f 4 | \
  {params.rd}/viterbi  /dev/stdin $mean {outdir}/{params.contig_prefix} 100 {params.eps} 0

"""
## viterbi params ?



rule orderVitter:
    input:
        CopyNumber="{outdir}/{contig}.viterout.txt",
        bins="{outdir}/coverage.bins.bed.gz",
    output:
        cn="{outdir}/copy_number.{contig}.bed",
    shell:"""

paste <(tabix  {input.bins} {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule combineVitter:
    input:
        allCopyNumberBED=expand("{outdir}/copy_number.{contig}.bed", contig=contigs,outdir=outdir),
    output:
        allCN="{outdir}/copy_number.tsv",
    run:
        import sys
        sys.stderr.write("Writing " + str(output.allCN) + "\n")
        cn=open(output.allCN,'w')
        for fn in input.allCopyNumberBED:
            sys.stderr.write("writing " + str(fn) + "\n")
            f=open(fn)
            l=f.readlines()
            cn.write("".join(l))
        cn.close()


rule call:
    input:
        allCN="{outdir}/copy_number.tsv",
        #aln=config["aln"],
        avg="{outdir}/mean_cov.txt",
    output:
        int="{outdir}/total_calls.bed",
        call="{outdir}/DUPcalls.copy_number.tsv",
        delcall="{outdir}/DELcalls.copy_number.tsv",
    params:
        rd=RD,
        genome_prefix=prefix_bam,
    shell:"""
cat {input.allCN} | awk '
BEGIN{{OFS="\t";s=0;e=0;cns=0;o=3;}}
{{
    if(NR==1)
    {{
        s=$1"\t"$2; e=$3; cns=$5;
        chr = $1;
    }}
    else if ($5==cns && chr==$1)
    {{
        e=$3;
        o=0;
    }}
    else
    {{
        print s,e,cns;
        s=$1"\t"$2;
        e=$3;
        cns=$5;
        chr=$1;
        o=1;
    }}
}}
END{{
if(o==0){{print s,e,cns;}}
}} ' > {output.int}

awk '$4<2' {output.int} > {output.delcall}
awk '$4>2' {output.int} > {output.call}

"""
mask="Yes"
#print(mask)
if mask != "No":
    rule repeatMask:
        input:
            call="{outdir}/DUPcalls.copy_number.tsv",
        output:
            maskedcall="{outdir}/DUPcalls.masked_CN.tsv",
            inter="{outdir}/DUPcalls.rep_int.bed",
        params:
            rep=REP,#config['repeatMask'],
            mnl=config["minL"],
            rd=RD,
        shell:"""

        intersectBed -wa -wb -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |sort -k1,1 -k2,2n | python {params.rd}/repeatMask.py | groupBy -g 1,2,3,10 -c 9| awk 'BEGIN{{OFS="\t"}} $6=$5/$4' > {output.inter}

    intersectBed -v -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |awk 'BEGIN{{OFS="\t"}}{{print$1,$2,$3,$3-$2,0,0}}'>>{output.inter}

    intersectBed -wa -wb -a <(awk '$6<0.8' {output.inter} ) -b <( awk '$4>2' {input.call} ) |awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$10,$6;}}'|awk '$3-$2>{params.mnl}' > {output.maskedcall}
    """

else:
    rule repeatMask0:
        input:
            call="{outdir}/DUPcalls.copy_number.tsv",
        output:
            maskedcall="{outdir}/DUPcalls.masked_CN.tsv",
        params:
            rd=RD,
        shell:"""
        cp {input.call} {output.maskedcall}
        """


rule compositeCall:
    input:
        call="{outdir}/DUPcalls.copy_number.tsv",
    output:
        compCall="{outdir}/DUPcalls.composite.bed",
    params:
        rd=RD,
        asm=ASM,#config["index"],
        mnl=config["minL"],
    shell:"""
mergeBed -i  {input.call} |sort -k1,1 -k2,2n |intersectBed -wb -b stdin -a {input.call} |groupBy -g 5,6,7 -c 2,3,4 -o collapse,collapse,collapse|awk '$3-$2 >{params.mnl}' >{output.compCall}
"""


mask="Yes"
if mask != "No":
#    repp=config["repeatMask"] #"{RD}/annotation/repeatMask.merged.bed",
    rule repeatMask2:
        input:
            call="{outdir}/DUPcalls.composite.bed",
        output:
            maskedCompcall="{outdir}/DUPcalls.masked_CN.composite.tsv",
            interC="{outdir}/DUPcalls.rep_int.composite.bed",
        params:
            rep=REP,#config["repeatMask"],
            rd=RD,
        shell:"""

        intersectBed -wa -wb -a <( cat {input.call} |cut -f 1-3) -b {params.rep} |sort -k1,1 -k2,2n | python {params.rd}/repeatMask.py | groupBy -g 1,2,3,10 -c 9| awk 'BEGIN{{OFS="\t"}} $6=$5/$4' > {output.interC}

    intersectBed -v -a <( cat {input.call} |cut -f 1-3) -b {params.rep} |awk 'BEGIN{{OFS="\t"}}{{print$1,$2,$3,$3-$2,0,0}}'>>{output.interC}

    intersectBed -wa -wb -a <(awk '$6<0.8' {output.interC} ) -b <( cat {input.call} ) |awk 'BEGIN{{OFS="\t"}} {{print $7,$8,$9,$10,$12,$6;}}' > {output.maskedCompcall}
    """

else:
    rule repeatMask3:
        input:
            call="{outdir}/DUPcalls.copy_number.tsv",
        output:
            maskedcall="{outdir}/DUPcalls.masked_CN.tsv",
        params:
            rd=RD,
        shell:"""
        cp {input.call} {output.maskedcall}
        """



rule PlotBins:
    input:
        allCN="{outdir}/copy_number.tsv",
        #aln=config["aln"],
        avg="{outdir}/mean_cov.txt",
    output:
        plot="{outdir}/{prefix_bam}.noclip.{ep}.pdf",
    params:
        rd=RD,
        genome_prefix="{prefix_bam}",
    shell:"""
touch {output.plot}
#plot every 50000 points ~ 5MB
$sum/../anaconda3/bin/Rscript {params.rd}/plot.HMM.noclip.R {input.allCN} {params.genome_prefix} 50000 {input.avg}
touch {output.plot}
"""




rule callSummary:
    input:
        call="{outdir}/DUPcalls.masked_CN.tsv",
    output:
        sumcall="{outdir}/CallSummary.{ep}.tsv",
        Vsumcall="{outdir}/CallSummary.verbose.{ep}.tsv",
    shell:"""
sort -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 4 -c 4,6 -o count,mean>{output.sumcall}
sort -k1,1 -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 1,4 -c 4,6 -o count,mean> {output.Vsumcall}
    """


rule GeneCount:
    input:
        call="{outdir}/DUPcalls.masked_CN.composite.tsv",
    output:
        geneCount="{outdir}/DUP.gene_count.bed",
    params:
        gen=GEN,
    shell:"""
        intersectBed -F 1 -wb -wa -a {input} -b {params.gen} | groupBy -g 1,2,3,6 -c 5,10,10 -o collapse,collapse,count > {output}
        """
