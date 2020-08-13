#!/bin/bash

# Get all files to be processed
selected_files=$(ls ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/ | grep "DLG2\|GFP" -)
selected_files=$(ls ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/ | grep "DLG2\|GFP" - | grep "_1\." - | sed 's/_1\..*//g' -)

# Process one by one
for f in $selected_files; 
do 
	echo 'Processing' $f;
  file1=$f'_1.clean.fq.gz'
  file2=$f'_2.clean.fq.gz'
  
  # Check whether processed already and skip if it does
  if [ -f ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts' ];
  then
      echo "~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts' exists already, skipping ...";
      continue;
  fi;
  
  # Align to grch38 using hisat: +/-1h
	echo '  Aligning ...';
  (~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/$file1 -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/$file2 | samtools view -bS - > ~/Projects/pub/DLG2/raw/$f'.bam') >& ~/Projects/pub/DLG2/raw/hisat2_log/$f'_log.txt'

  # Quantify to gencode 29 using htseq-count: +/-2h30
	echo '  Quantifying ...';
  # (samtools view ~/Projects/pub/DLG2/raw/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_name - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/htseq_log/$f'_log.txt'
  (samtools view ~/Projects/pub/DLG2/raw/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_id - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/htseq_log/$f'_log.txt'
  
  # Get bam stats: includes insert sizes, read lengths, ...
  samtools stats ~/Projects/pub/DLG2/raw/$f'.bam' > ~/Projects/pub/DLG2/raw/bam_stat/$f'_stat_log.txt'
  
  # Rm bam file
	echo '  Finishing ...';
  rm ~/Projects/pub/DLG2/raw/$f'.bam'
  echo ''

done
