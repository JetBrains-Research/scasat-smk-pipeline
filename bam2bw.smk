rule bam2bdg:
    input:
        bam="{anywhere}/{filename}.bam",
        chrom_sizes="hg19.chrom.sizes"
    output: temp("{anywhere}/bdg/{filename, [^/]*}.bdg")
    conda: "bio.environment.yml"
    shell: 'bedtools genomecov -ibam {input.bam} -bg -g {input.chrom_sizes} > {output}'


rule bdg2clip:
    input:
        bdg="{anywhere}/bdg/{filename}.bdg",
        chrom_sizes="hg19.chrom.sizes"
    output: temp("{anywhere}/clip/{filename, [^/]*}.clip")
    conda: "bio.environment.yml"
    shell: 'bedtools slop -i {input.bdg} -g {input.chrom_sizes} -b 0 | '
           'bedClip stdin {input.chrom_sizes} {output}'


rule sort_clip:
    input: "{anywhere}/clip/{filename}.clip"
    output: temp("{anywhere}/clip/{filename, [^/]*}.clip.sort")
    shadow: "shallow"
    shell: 'LC_COLLATE=C sort -k1,1 -k2,2n -T . {input} > {output}'


rule clip2bw:
    input:
        clip="{anywhere}/clip/{filename}.clip.sort",
        chrom_sizes="hg19.chrom.sizes"
    output: "{anywhere}/bw/{filename, [^/]*}.bw"
    conda: "bio.environment.yml"
    shell: 'bedGraphToBigWig {input.clip} {input.chrom_sizes} {output}'
