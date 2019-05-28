rule bam2bw:
    input:
        bam="{anywhere}/{filename}.bam",
        bai="{anywhere}/{filename}.bam.bai"
    output: "{anywhere}/bw/{filename, [^/]*}.bw"
    conda: "envs/deeptools.environment.yaml"
    shell: 'bamCoverage -b {input.bam} -o {output}'