# sc-atacseq-smk-pipeline
Single cell ATAC-Seq Snakemake pipeline

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
