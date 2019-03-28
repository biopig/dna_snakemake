#这是一个材料的基因组重测序的数据分析流程，使用的是双端测序
#首先需要先下载参考基因组，并构建索引

# by：詹为民
# email：630950832@qq.com
# data：2019.3.24
#需要的软件有：fastqc, trimmomatic, bowtie2, samtools, sambamba, qualimap, deeptools, GATK, VEP 

Vep="/home/zz/biosoft/vep/ensembl-vep-release-92/vep"
Gatk="~/biosoft/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"
Trimmomatic="~/biosoft/Trimmomatic-0.38/trimmomatic-0.38.jar"
Java="java -jar -Xss1G -Xmx5G -XX:ParallelGCThreads=3"

Index="/disks/reference/index/Zea_mays.AGPv4_index"
Reference="/disks/reference/index/Zea_mays.AGPv4.dna.toplevel.fa"
Reference_dir="/disks/reference/index/"
Known_vcf="/disks/reference/annotation/zea_mays.vcf"
Known_vcf_dir="/disks/reference/annotation/"
Cache_dir="/disks/reference/vep/"
Truseq="/home/zz/biosoft/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa"

samples=["text"]

ALL_fastqc=expand("1_fastqc/{sample}_1_fastqc.zip",sample=samples)

ALL_trimmomatic=[]
for Sample in samples:
    ALL_trimmomatic.append("2_trimm/{}_1_paired.fq.gz".format(Sample))
    ALL_trimmomatic.append("2_trimm/{}_1_unpaired.fq.gz".format(Sample))
    ALL_trimmomatic.append("2_trimm/{}_2_paired.fq.gz".format(Sample))
    ALL_trimmomatic.append("2_trimm/{}_2_unpaired.fq.gz".format(Sample))

ALL_bowtie2=expand("3_mapping/{sample}.bam",sample=samples)
ALL_sort=expand("4_sambamba/{sample}_sort.bam",sample=samples)
ALL_markdup=expand("4_sambamba/{sample}_rmdup.bam",sample=samples)
ALL_qualimap_summary=expand("5_summary/{sample}_rmdup_stats/qualimapReport.html",sample=samples)
ALL_bam_coverage=expand("5_summary/{sample}_coverage.pdf",sample=samples)
ALL_transform=expand("5_summary/{sample}_bamcoverage.bw",sample=samples)
ALL_addhead=expand("6_GATK/{sample}_addhead.bam",sample=samples)
ALL_recalibration=expand("6_GATK/{sample}_baserecalibrator.list",sample=samples)
ALL_applyBQSR=expand("6_GATK/{sample}_applybqsr.bam",sample=samples)
ALL_haplotypecaller=expand("6_GATK/{sample}_haplotypecaller.vcf",sample=samples)
ALL_vep=expand("7_VEP/{sample}_vep_annotate.vcf",sample=samples)

TARGET=[]
TARGET.extend(ALL_fastqc)
TARGET.extend(ALL_trimmomatic)
TARGET.extend(ALL_bowtie2)
TARGET.extend(ALL_sort)
TARGET.extend(ALL_markdup)
TARGET.extend(ALL_qualimap_summary)
TARGET.extend(ALL_bam_coverage)
TARGET.extend(ALL_transform)
TARGET.extend(ALL_addhead)
TARGET.extend(ALL_recalibration)
TARGET.extend(ALL_applyBQSR)
TARGET.extend(ALL_haplotypecaller)
TARGET.extend(ALL_vep)

rule all:
    input:
        TARGET

rule fastqc:
    input:
        "raw_data/{sample}_1.fq.gz",
        "raw_data/{sample}_2.fq.gz"
    output:
        "1_fastqc/{sample}_1_fastqc.zip"
    threads:
        4
    shell:
        "fastqc -o 1_fastqc -t {threads} -f fastq {input}"
 
rule trimmomatic:
    input:
        "1_fastqc/{sample}_1_fastqc.zip",
        "raw_data/{sample}_1.fq.gz",
        "raw_data/{sample}_2.fq.gz"
    output:
        "2_trimm/{sample}_1_paired.fq.gz",
        "2_trimm/{sample}_1_unpaired.fq.gz",
        "2_trimm/{sample}_2_paired.fq.gz",
        "2_trimm/{sample}_2_unpaired.fq.gz"
    shell:
        "{Java} {Trimmomatic} PE {input[1]} {input[2]} {output} ILLUMINACLIP:{Truseq}:2:30:10:8:True SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50"

rule bowtie2:
    input:
        "2_trimm/{sample}_1_paired.fq.gz",
        "2_trimm/{sample}_1_unpaired.fq.gz",
        "2_trimm/{sample}_2_paired.fq.gz",
        "2_trimm/{sample}_2_unpaired.fq.gz"
    output:
        "3_mapping/{sample}.bam"
    threads:
        4
    shell:
        "bowtie2 -p {threads} -x {Index} -1 {input[0]} -2 {input[2]} -U {input[1]} -U {input[3]} |samtools view -S -O bam -q 30 -F 4 -o {output} -"
#在转换的时候就把质量值低的和没有比对上的去掉

#sambamba在运行完后会自动构建索引的
rule sambamba_sort:
    input:
        "3_mapping/{sample}.bam"
    output:
        "4_sambamba/{sample}_sort.bam"
    threads:
        4
    shell:
        "sambamba sort -t {threads} -o {output} {input}"

rule markdup:
    input:
        "4_sambamba/{sample}_sort.bam"
    output:
        "4_sambamba/{sample}_rmdup.bam"
    threads:
        4
    shell:
        "sambamba markdup -r -t {threads} {input} {output}"
        
rule qualimap_summary:
    input:
        "4_sambamba/{sample}_rmdup.bam"
    output:
        "5_summary/{sample}_rmdup_stats/qualimapReport.html"
    params:
        directory("5_summary/{sample}_rmdup_stats")
    shell:
        "qualimap bamqc -bam {input} -outdir {params} -c -nw 400 -hm 3"

rule bam_coverage:
    input:
        "4_sambamba/{sample}_rmdup.bam",
        "5_summary/{sample}_rmdup_stats/qualimapReport.html"
    output:
        "5_summary/{sample}_coverage.pdf"
    threads:
        4
    shell:
        "plotCoverage -b {input[0]} -o {output} --plotFileFormat pdf -p {threads}"

rule deeptools_transform:
    input:
        "4_sambamba/{sample}_rmdup.bam",
        "5_summary/{sample}_coverage.pdf"
    output:
        "5_summary/{sample}_bamcoverage.bw"
    threads:
        4
    shell:
        "bamCoverage -b {input[0]} -o {output} -of bigwig -p {threads}" 

rule GATK_addgroups:
    input:
        "4_sambamba/{sample}_rmdup.bam"
    output:
        "6_GATK/{sample}_addhead.bam"
    params:
        LB="{sample}ID",
        PU="{sample}PU",
        SM="{sample}"
    shell:
        "{Java} {Gatk} AddOrReplaceReadGroups -I {input} -O {output} -LB {params.LB} -PL illumina -PU {params.PU} -SM {params.SM}"

rule GATK_recalibration:
    input:
        "6_GATK/{sample}_addhead.bam"
    output:
        "6_GATK/{sample}_baserecalibrator.list"
    shell:
        "{Java} {Gatk} BaseRecalibrator -I {input} --known-sites {Known_vcf} -O {output} -R {Reference} "

rule GATK_applyBQSR:
    input:
        "6_GATK/{sample}_baserecalibrator.list",
        "6_GATK/{sample}_addhead.bam"
    output:
        "6_GATK/{sample}_applybqsr.bam"
    shell:
        "{Java} {Gatk} ApplyBQSR -bqsr {input[0]} -I {input[1]} -O {output}"

rule GATK_HaplotypeCaller:
    input:
        "6_GATK/{sample}_applybqsr.bam"
    output:
        "6_GATK/{sample}_haplotypecaller.vcf"
    shell:
        "{Java} {Gatk} HaplotypeCaller --emit-ref-confidence GVCF -I {input} -O {output} -R {Reference}"

rule VEP_annotation:
    input:
        "6_GATK/{sample}_haplotypecaller.vcf"
    output:
        "7_VEP/{sample}_vep_annotate.vcf"
    shell:
        "{Vep} --fasta {Reference} --offline --species zea_mays --cache_version 39 --dir_cache {Cache_dir} -i {input} -o {output}"
