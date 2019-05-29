import os
import re
from glob import glob

configfile: "config.yaml"


print('CONFIG\n{}'.format('\n'.join(['{}: {}'.format(k, v) for k, v in config.items()])))


def fastq_files():
    return glob(os.path.join(config["fastq_dir"], "*.f*q"))


def fastq_names():
    return [os.path.splitext(os.path.basename(fastq_file))[0] for fastq_file in fastq_files()]


def fastq_common_names_paired():
    basenames = {os.path.splitext(fastq_file)[0] for fastq_file in fastq_files()}
    paired_candidates = [basename[:-2] for basename in basenames if basename[-2:] == "_1"]
    return [os.path.basename(common_name) for common_name in paired_candidates if common_name + "_2" in basenames]


def fastq_names_single():
    common_names = fastq_common_names_paired()
    return [fastq_name for fastq_name in fastq_names()
            if fastq_name[-2:] not in ["_1", "_2"] or fastq_name[:-2] not in common_names]


def fastq_aligned_names():
    return fastq_common_names_paired() + fastq_names_single()


def cell_names_dict():
    cell_names = {re.sub('^GSM[0-9]+_', '', fastq_aligned_name) for fastq_aligned_name in fastq_aligned_names()}
    return {cell_name: [fastq_aligned_name for fastq_aligned_name in fastq_aligned_names()
                        if re.match('^GSM[0-9]+_{cell_name}$'.format(cell_name=cell_name),
                                    fastq_aligned_name)]
            for cell_name in cell_names}


workdir: config["work_dir"]

rule fastqc:
    input: os.path.join(config["fastq_dir"], "{sample}.fastq")
    output:
          html="qc/fastqc/{sample}_fastqc.html",
          zip="qc/fastqc/{sample}_fastqc.zip"
    log: "logs/fastqc/{sample}.log"
    wrapper: "0.31.1/bio/fastqc"

rule multiqc_fastq:
    input: expand("qc/fastqc/{sample}_fastqc.zip", sample=fastq_names())
    output: "multiqc/fastqc/multiqc.html"
    log: "multiqc/fastqc/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"

rule trim_adapters:
    input:
         r1=os.path.join(config["fastq_dir"], "{sample}_1.fastq"),
         r2=os.path.join(config["fastq_dir"], "{sample}_2.fastq")
    output:
          r1="trimmed/{sample}_1.fastq.gz",
          r2="trimmed/{sample}_2.fastq.gz",
          r1_unpaired="trimmed/{sample}_1_unpaired.fastq.gz",
          r2_unpaired="trimmed/{sample}_2_unpaired.fastq.gz"
    params:
          trimmer=["TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:25"],
          extra="-phred33 -threads 1"
    log: "logs/trimmomatic/{sample}.log"
    wrapper: "0.31.1/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input: "trimmed/{sample}.fastq.gz"
    output:
          html="qc/fastqc_trimmed/{sample}_fastqc.html",
          zip="qc/fastqc_trimmed/{sample}_fastqc.zip"
    log: "logs/fastqc_trimmed/{sample}.log"
    wrapper: "0.31.1/bio/fastqc"

rule multiqc_fastq_trimmed:
    input: expand("qc/fastqc_trimmed/{sample}_fastqc.zip", sample=fastq_names())
    output: "multiqc/fastqc_trimmed/multiqc.html"
    log: "multiqc/fastqc_trimmed/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"

rule bowtie2:
    input:
         sample=["trimmed/{sample}_1.fastq.gz", "trimmed/{sample}_2.fastq.gz"]
    output: "mapped/{sample}.bam"
    log: "logs/bowtie2/{sample}.log"
    params:
          index=config["indexes"],
          extra="-X 2000 --dovetail"
    threads: 8
    wrapper: "0.31.1/bio/bowtie2/align"

rule multiqc_bowtie2:
    input: expand("logs/bowtie2/{sample}.log", sample=fastq_aligned_names())
    output: "multiqc/bowtie2/multiqc.html"
    log: "multiqc/bowtie2/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"

rule samtools_filter:
    input: "mapped/{sample}.bam"
    output: "filtered/{sample}.bam"
    params: "-bhF 4 -f2 -q30"
    wrapper: "0.31.1/bio/samtools/view"

rule download_blacklist:
    output: "blacklist/blacklist.bed"
    shell:
         'wget -O - "{config[blacklist]}" | gunzip -c > {output}'

rule remove_blacklist:
    input:
         bam="filtered/{sample}.bam",
         blacklist=rules.download_blacklist.output
    output: "blacklisted/{sample}.bam"
    shell: 'intersectBed -v -abam {input.bam} -b {input.blacklist} > {output}'

rule sort_bams:
    input: "blacklisted/{sample}.bam"
    output: "sorted/{sample}.bam"
    threads: 4
    wrapper: '0.31.1/bio/samtools/sort'

rule remove_duplicates:
    input: "sorted/{sample}.bam"
    output:
          bam="deduplicated/{sample}.bam",
          metrics="qc/picard/{sample}.txt"
    log: "deduplicated/{sample}.log"
    params: "REMOVE_DUPLICATES=True"
    wrapper: "0.31.1/bio/picard/markduplicates"

rule sort_deduplicated_bams:
    input: "deduplicated/{sample}.bam"
    output: "deduplicated_sorted/{sample}.bam"
    threads: 4
    wrapper: '0.31.1/bio/samtools/sort'

rule index_bams:
    input: "{anywhere}/{sample}.bam"
    output: "{anywhere}/{sample, [^/]*}.bam.bai"
    wrapper: "0.31.1/bio/samtools/index"

rule bam2bw:
    input:
        bam="{anywhere}/{filename}.bam",
        bai="{anywhere}/{filename}.bam.bai"
    output: "{anywhere}/bw/{filename, [^/]*}.bw"
    conda: "envs/deeptools.environment.yaml"
    threads: 4
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output}'

rule clean_bams:
    input:
         bam="deduplicated_sorted/{sample}.bam",
         bai="deduplicated_sorted/{sample}.bam.bai"
    output: "cleaned/{sample}.bam"
    params:
         chrs=' '.join(config['chromosomes'])
    conda: "envs/samtools.env.yaml"
    shell:
         'samtools view -b {input.bam} {params.chrs} > {output}'

rule bam_stats:
    input: "cleaned/{sample}.bam"
    output: "qc/samtools_stats/{sample}_samtools_stats.txt"
    wrapper: "0.31.1/bio/samtools/stats"

rule bam_stats_multiqc:
    input: expand("qc/samtools_stats/{sample}_samtools_stats.txt", sample=fastq_aligned_names())
    output: "multiqc/samtools_stats/multiqc.html"
    log: "multiqc/samtools_stats/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"

rule merge_cells:
    input: lambda wildcards: expand("cleaned/{sample}.bam", sample=cell_names_dict()[wildcards.cell_name])
    output: "cleaned_cells/{cell_name}.bam"
    threads: 4
    wrapper: "0.31.1/bio/samtools/merge"

rule merge_all:
    input: expand("cleaned/{sample}.bam", sample=fastq_aligned_names())
    output: "cleaned_all/pooled.bam"
    threads: 8
    wrapper: "0.31.1/bio/samtools/merge"

rule sort_cleaned_bams_cells:
    input: "cleaned_cells/{cell_name}.bam"
    output: "cleaned_cells_sorted/{cell_name}.bam"
    threads: 4
    wrapper: '0.31.1/bio/samtools/sort'

rule sort_cleaned_bams_all:
    input: "cleaned_all/pooled.bam"
    output: "cleaned_all_sorted/pooled.bam"
    threads: 4
    wrapper: '0.31.1/bio/samtools/sort'

rule call_peaks_macs2_cell:
    input: "cleaned_cells/{cell_name}.bam"
    output:
          "cleaned_cell_peaks/macs2/{cell_name}_{macs2_suffix}_peaks.narrowPeak",
          "cleaned_cell_peaks/macs2/{cell_name}_{macs2_suffix}_summits.bed"
    params:
          outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
          macs2_stats=config.get('macs2_stats', '-p 0.0001')
    conda: "envs/macs2.env.yaml"
    shell:
         'macs2 callpeak -t {input} --outdir {params.outdir} -n {wildcards.cell_name}_{wildcards.macs2_suffix} ' \
         '{params.macs2_stats} -g hs -f BAMPE --nomodel --nolambda -B --keep-dup all --call-summits'

rule call_peaks_macs2_all:
    input: "cleaned_all/pooled.bam"
    output:
          "cleaned_all_peaks/macs2/pooled_{macs2_suffix}_peaks.narrowPeak",
          "cleaned_all_peaks/macs2/pooled_{macs2_suffix}_summits.bed"
    params:
          outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
          macs2_stats=config.get('macs2_stats', '-p 0.0001')
    conda: "envs/macs2.env.yaml"
    shell:
         'macs2 callpeak -t {input} --outdir {params.outdir} -n pooled_{wildcards.macs2_suffix} '
         '{params.macs2_stats} -g hs -f BAMPE --nomodel --nolambda -B --keep-dup all --call-summits'

rule download_span:
    output: "bin/span-0.11.0.jar"
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.11.0.4882.jar'

rule download_chrom_sizes:
    output: "{config[genome]}.chrom.sizes"
    shell: 'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes'

rule call_peaks_span_cell:
    input:
         bam="cleaned_cells/{cell_name}.bam",
         span=rules.download_span.output,
         chrom_sizes=rules.download_chrom_sizes.output
    output: "cleaned_cell_peaks/span/{cell_name}_{bin, [0-9]+}.span"
    params:
          outdir=lambda wildcards, output: os.path.dirname(str(output)),
          xmx=lambda wildcards: str(800 // int(wildcards.bin))
    threads: 8
    shell:
         'java -Xmx{params.xmx}G -jar {input.span} analyze -t {input.bam} --workdir {params.outdir} --fragment 0 '
         '--bin {wildcards.bin} --cs {input.chrom_sizes} --threads {threads} --model {output} --keep-dup true --debug; '

rule call_peaks_span_all:
    input:
         bam="cleaned_all/pooled.bam",
         span=rules.download_span.output,
         chrom_sizes=rules.download_chrom_sizes.output
    output: "cleaned_all_peaks/span/pooled_{bin, [0-9]+}.span"
    params:
          outdir=lambda wildcards, output: os.path.dirname(str(output)),
          xmx=lambda wildcards: str(800 // int(wildcards.bin))
    threads: 8
    shell:
         'java -Xmx{params.xmx}G -jar {input.span} analyze -t {input.bam} --workdir {params.outdir} --fragment 0 '
         '--bin {wildcards.bin} --cs {input.chrom_sizes} --threads {threads} --model {output} --keep-dup --debug; '

rule all:
    input:
         expand("cleaned_cell_peaks/macs2/{cell_name}_{macs2_suffix}_peaks.narrowPeak",
                cell_name=cell_names_dict().keys(), macs2_suffix=config.get('macs2_suffix', 'p0.0001')),

         expand("cleaned_all_peaks/macs2/pooled_{macs2_suffix}_peaks.narrowPeak",
                macs2_suffix=config.get('macs2_suffix', 'p0.0001')),

         expand("cleaned_cell_peaks/span/{cell_name}_100.span", cell_name=cell_names_dict().keys()),
         "cleaned_all_peaks/span/pooled_100.span",

         expand("cleaned_cells_sorted/bw/{cell_name}.bw", cell_name=cell_names_dict().keys()),
         "cleaned_all_sorted/bw/pooled.bw"
