# sc-atacseq-smk-pipeline
Single cell ATAC-Seq Snakemake pipeline.

Adapted from [Scasat](https://github.com/ManchesterBioinference/Scasat/blob/master/Scasat_Pre-process_rep1_BAMPE.ipynb).

We assume that the genome build is `hg19`. The pipeline can be easily adapted to
other builds and to other organisms.

Data preprocessing
--------------------------------------
Example: [GSE74310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74310).\
Download `SRR.txt` and `SraRunTable.txt` files from NCBI SRA Run Selector [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSEXXXXXX&go=go](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE74310&go=go)

Download fastq files
--------------------
```bash
cat SRR.txt | while read -r ID; do echo $ID; fastq-dump --split-files $ID; done
```

Renaming files
--------------
Here we assume that each GSM file corresponds to a single SRR file.
```bash
for FILE in $(find . -name "*.fastq"); do 
  SRR=$(echo ${FILE} | sed 's#./##g' | sed 's#_1.fastq##g' | sed 's#_2.fastq##g'); 
  NAME=$(cat SraRunTable.txt | grep $SRR | awk -v FS='\t' '{ printf("%s_%s", $9, $10)}' |\
    sed 's#/#_#g' | sed 's#-#_#g' | sed 's#,##g' | sed 's# #_#g' ); 
  SUFFIX=$(echo $FILE | sed 's#.*SRR.*_##g'); 
  echo "${FILE} -> ${NAME}_${SUFFIX}"; 
  mv ${FILE} ${NAME}_${SUFFIX}; 
done
```

Launching Pipeline
------------------
The only tool required to launch the pipeline is `conda`. You'll also need `bowtie2`
indexes for the appropriate genome build.
* If `conda` is not installed,
follow the instructions at
[Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* Navigate to repository directory.
```bash
$ cd sc-atacseq-smk-pipeline
```
Create a Conda environment for `snakemake`:
```bash
$ conda env create --file envs/snakemake.yaml --name snakemake
```
Activate the newly created environment:
```bash
$ source activate snakemake
```
Run the pipeline:
```bash
$ snakemake all [--cores <cores>] --use-conda --config work_dir=<work_dir> fastq_dir=<fastq_dir> indexes=<bowtie2_indexes>

```
Here `<work_dir>` is the directory where the output and temporary files will be created,
`<fastq_dir>` is the directory that houses the raw FASTQ files that have been downloaded
and renamed in the previous steps, and `<bowtie2_indexes>` is the path to `bowtie2` index
files together with the genome build.
For example, if the indexes reside at `/home/user/index`
and are named `hg19.1.bt2` etc., then you should substitute `<bowtie2_indexes>` with the
exact string `/home/user/index/hg19`.

The pipeline can take a significant time to complete. To speed it up, you can provide
more computational power using `<cores>` parameter. This will allow the pipeline
to launch independent jobs in parallel, and to run some jobs in multithreaded mode.

Pipeline Highlights
-------------------

- trimming adapters (`trimmomatic`)
- aligning (`bowtie2`)
- filtering mappings by quality
- deduplication (`picard`)
- filtering by [ENCODE mapability blacklist](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability)
- pooling by cell types and the entire dataset
- calling peaks (`macs2`)
- generating quality reports (`multiqc`)

Pipeline Output
---------------
Notable output files for further downstream analysis:
- quality control:
  - `multiqc/fastqc/multiqc.html` -- FASTQC summary for the original dataset
  (see `qc/fastqc` for individual reports)
  - `multiqc/fastqc_trimmed/multiqc.html` -- FASTQC summary after trimming adapters
  (see `qc/fastqc_trimmed` for individual reports)
  - `multiqc/bowtie2/multiqc.html` -- alignment summary
  (see `logs/bowtie2` for individual reports)
  - `multiqc/samtools_stats/multiqc.html` -- BAM summary after deduplication, cleaning
  and filtering (see `qc/samtools_stats` for individual reports)
- aligned reads:
  - `cleaned` -- coordinate-sorted BAM files for individual single cells
  - `cleaned_cells_sorted` -- coordinate-sorted BAM files for pooled cell type data
  - `cleaned_all_sorted` -- coordinate-sorted BAM file for total pooled data  
- visualization:
  - `cleaned_cells_sorted/bw` -- bigWig files for pooled cell type data
  - `cleaned_all_sorted/bw` -- bigWig file for total pooled data
- peaks:
  - `cleaned_cell_peaks/macs2` -- MACS2 peaks for pooled cell type data 
  - `cleaned_all_peaks/macs2` -- MACS2 peaks for total pooled data