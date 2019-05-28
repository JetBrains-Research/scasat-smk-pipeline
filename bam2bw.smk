rule bam2bw:
    input: "{anywhere}/{filename}.bam"
    output: "{anywhere}/bw/{filename, [^/]*}.bw"
    conda: "envs/deeptools.environment.yaml"
    shell: 'bamCoverage -b {input} -o {output}'