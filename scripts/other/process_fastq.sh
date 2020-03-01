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

# Check in IGV whether differences between RA and DLG2 on DLG2?
# selected_files=$(ls ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/ | grep "DLG2_5Y1\|GFP_5Y5" - | grep "_1\." - | sed 's/_1\..*//g' -)
selected_files=$(ls ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/ | grep "GFP_5Y1\|DLG2_5Y5" - | grep "_1\." - | sed 's/_1\..*//g' -)

# Process one by one
for f in $selected_files; 
do 
	echo 'Processing' $f;
  file1=$f'_1.clean.fq.gz'
  file2=$f'_2.clean.fq.gz'
  
  # # Check whether processed already and skip if it does
  # if [ -f ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts' ];
  # then
  #     echo "~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts' exists already, skipping ...";
  #     continue;
  # fi;
  #
  # Align to grch38 using hisat: +/-1h
	echo '  Aligning ...';
  (~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/$file1 -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/$file2 | samtools view -bS - > ~/Projects/pub/DLG2/raw/$f'.bam')

  # Sort and index
  samtools sort ~/Projects/pub/DLG2/raw/$f'.bam'  -o ~/Projects/pub/DLG2/raw/$f'.bam'
  samtools index ~/Projects/pub/DLG2/raw/$f'.bam'

	# Limit to DLG2 region
	samtools view -H raw/$f'.bam'  > raw/$f'_chr11.sam' # Output header
	samtools view raw/$f'.bam' "11:82345532-86838367" >> raw/$f'_chr11.sam' # extract region
	samtools view -bS raw/$f'_chr11.sam' > raw/$f'_chr11.bam' # Convert to bam
	samtools sort raw/$f'_chr11.bam' > raw/$f'_chr11_sorted.bam' # Sort
	samtools index raw/$f'_chr11_sorted.bam' # Index


#   # Quantify to gencode 29 using htseq-count: +/-2h30
# 	echo '  Quantifying ...';
#   # (samtools view ~/Projects/pub/DLG2/raw/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_name - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/htseq_log/$f'_log.txt'
#   (samtools view ~/Projects/pub/DLG2/raw/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_id - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/htseq_log/$f'_log.txt'
#   
#   # Get bam stats: includes insert sizes, read lengths, ...
#   samtools stats ~/Projects/pub/DLG2/raw/$f'.bam' > ~/Projects/pub/DLG2/raw/bam_stat/$f'_stat_log.txt'
#   
#   # Rm bam file
# 	echo '  Finishing ...';
#   rm ~/Projects/pub/DLG2/raw/$f'.bam'
#   echo ''

done

# Manually moves files to bam_DLG2/, only keep region specific ones



# 
# 
# # 
# # echo $1: starting
# # 
# # 
# # # Bowtie2: used before
# # start=`date +%s`
# # bowtie2 -p 8 -x ~/downloads/data/genomes/hg19_ind -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz -S ~/Projects/pub/DLG2/raw/test_file.sam >& ~/Projects/pub/DLG2/raw/test_file_log.txt
# # end=`date +%s`
# # runtime=$((end-start))
# # 
# # #1h10' but only 75% ???
# # 
# # # Retry with indexes from igenome
# # start=`date +%s`
# # bowtie2 -p 8 -x ~/downloads/data/genomes/igenomes/hg19 -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz -S ~/Projects/pub/DLG2/raw/test_file.sam >& ~/Projects/pub/DLG2/raw/test_file_log2.txt
# # end=`date +%s`
# # runtime=$((end-start))
# # 
# # #Again the same!
# # 
# # # TRY hisat: used by ari now
# # hisat -p 4 --no-mixed --no-discordant --rna-strandness RF --no-unal --known-splicesite-infile splicesites.txt -x ~/erik_home/erik/genomes/hg19_ucsc_minus_hap -1 180403_NB501037_0252_AHLTNGBGX5/Fastq/KE210218-$1_R1_001.fastq.gz,180404_NB501037_0253_AHLTVTBGX5/Fastq/KE210218-$1_R1_001.fastq.gz -2 180403_NB501037_0252_AHLTNGBGX5/Fastq/KE210218-$1_R2_001.fastq.gz,180404_NB501037_0253_AHLTVTBGX5/Fastq/KE210218-$1_R2_001.fastq.gz | samtools view -bS - > bam/$1.bam) >& logs/$1.hisat_log
# # # They recommend Hisat2? Try with hg38?
# # # See http://ccb.jhu.edu/software/hisat2/index.shtml
# 
# # HISAT 2
# # Dowload to hisat2 to tools/
# # Download index hg38 from igenome
# start=`date +%s`
# # ~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz -S ~/Projects/pub/DLG2/raw/test_file.sam >& ~/Projects/pub/DLG2/raw/test_file_hisat2_log.txt
# # 92.37%! and only  20 mins!!!
# # ~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz -S ~/Projects/pub/DLG2/raw/test_file.sam >& ~/Projects/pub/DLG2/raw/test_file_hisat2_log.txt
# # 90.5%
# (~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz | samtools view -bS - > ~/Projects/pub/DLG2/raw/test_file.bam) >& ~/Projects/pub/DLG2/raw/test_file_hisat2_bam_log.txt
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# # 1h5'
# 
# 
# # Sort and index
# start=`date +%s`
# samtools sort ~/Projects/pub/DLG2/raw/test_file.bam  -o ~/Projects/pub/DLG2/raw/test_file_sorted.bam 
# samtools index ~/Projects/pub/DLG2/raw/test_file_sorted.bam 
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# # +/- 20'
# # Visualize in IGV
# 
# # Quantify using HTSeq
# # Download latest gencode (gencode 29)
# # From Ari: (samtools view bam/$1.bam | htseq-count -m intersection-strict -s reverse - wgEncodeGencodeMerged_reformatted_LINC00942_fixed.gtf > results/$1.gene_counts) >& /dev/stdout | tail -n200 > logs/$1.htseq_log
# # More efficiet from sam?
# start=`date +%s`
# (samtools view ~/Projects/pub/DLG2/raw/test_file.bam | htseq-count -m intersection-strict -s reverse - ~/downloads/data/gencode29/gencode.v29.annotation.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts) >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/test_htseq_log
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# # 45' (2864 sec)
# # No counts in DESEQ2???
# # rename chr in gtf
# # sed 's/chr//g' gencode.v29.annotation.gtf > gencode.v29.annotation_noChr.gtf
# start=`date +%s`
# (samtools view ~/Projects/pub/DLG2/raw/test_file.bam | htseq-count -m intersection-strict -s reverse - ~/downloads/data/gencode29/gencode.v29.annotation_noChr.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts) >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/test_htseq_log
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# # 5733: +/-1h30
# # MT genome still a problem?
# # samtools view raw/test_file.sam | cut -f3 | sort -u --> MT, X, Y, ...
# # cat ~/downloads/data/gencode29/gencode.v29.annotation_noChr.gtf | cut -f1 | sort -u --> M
# # rename chrM in gtf
# # sed 's/M/MT/g' gencode.v29.annotation_noChr.gtf > gencode.v29.annotation_noChr_MT.gtf
# start=`date +%s`
# (samtools view ~/Projects/pub/DLG2/raw/test_file.bam | htseq-count -m intersection-strict -s reverse - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts) >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/test_htseq_log
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# 
# 
# # Normalize using DSeq2, also for diff expr!!!
# 
# # Faster without bam files?
# start=`date +%s`
# (~/tools/hisat2/hisat2-2.1.0/hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x ~/downloads/data/genomes/igenomes/grch38/genome -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz | htseq-count -m intersection-strict -s reverse - ~/downloads/data/gencode29/gencode.v29.annotation_noChr.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts2) >& ~/Projects/pub/DLG2/raw/test_file_log.txt
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# # 7883: 2h10'
# 
# # Have to be sorted!!!! ????
# # Add gene_name: --additional-attr gene_name --> Only in later versions of ht-seq count, just set as primary attribute for the moment
# start=`date +%s`
# # (samtools view ~/Projects/pub/DLG2/raw/test_file_sorted.bam | htseq-count -m intersection-strict -s reverse -i gene_id --additional-attr=gene_name - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts) >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/test_htseq_log
# (samtools view ~/Projects/pub/DLG2/raw/test_file.bam | htseq-count -m intersection-strict -s reverse -i gene_name - ~/downloads/data/gencode29/gencode.v29.annotation_noChr_MT.gtf > ~/Projects/pub/DLG2/raw/test.gene_counts) >& /dev/stdout | tail -n200 > ~/Projects/pub/DLG2/raw/test_htseq_log
# end=`date +%s`
# runtime=$((end-start))
# echo $runtime
# 
# # # TRY tophat2: used by ari before
# # start=`date +%s`
# # tophat2 -p 8 -r 100 -G /home/arghavan/bengts_project/wgEncodeGencodeMerged_reformatted.gtf -o ~/Projects/pub/DLG2/raw/test_tophat /home/arghavan/erik_home/erik/genomes/hg19_ucsc_minus_hap -1 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_1.clean.fq.gz -2 ~/Projects/pub/DLG2/raw/fastq_NOBACKUP/Ctrl_10_2.clean.fq.gz 
# # end=`date +%s`
# # runtime=$((end-start))
# 
# 
# 
# 
# 
# 
# echo $1: aligned
# 
# samtools view -b -q1 ~/bengts_project/tophat_out_$1/accepted_hits.bam > ~/bengts_project/bam/$1.bam
# echo $1: filtered
# 
# # Quantify GENCODE genes
# samtools sort -n ~/bengts_project/bam/$1.bam ~/bengts_project/bam/$1.sorted
# echo $1: sorted
# 
# samtools view ~/bengts_project/bam/$1.sorted.bam > ~/bengts_project/$1.sam
# echo $1: converted to SAM
# 
# htseq-count -m intersection-strict -s no ~/bengts_project/$1.sam ~/bengts_project/wgEncodeGencodeMerged_reformatted.gtf > ~/bengts_project/results/$1.gene_counts >& /dev/stdout | tail -n100 > ~/bengts_project/logs/$1.htseq_log
# echo $1: genes quantified
# 
# #sort and index the bam file for IGV:
# samtools sort ~/bengts_project/bam/$1.sorted.bam ~/bengts_project/bam/$1.sort
# samtools index ~/bengts_project/bam/$1.sort.bam
# 
# rm ~/bengts_project/bam/$1.bam
# rm ~/bengts_project/bam/$1.sorted.bam 
# rm ~/bengts_project/$1.sam
# 
# echo $1: processing complete
# 
# # ALTERNATIVE
# #!/bin/bash -l
# #all the gtf, bed files, etc. were copied form rudy/erik/TCGA_rnaseq/expr_realign/
# 
# echo $1: starting
# 
# # note: no-mixed, no-disc, no-unal makes very little difference for htseq in intersection-strict mode - but makes sense to do this anyway, such that coverageBed counts will be equally strict. --rna-strandness RF is used because the data is strand-specific!
# (hisat -p 4 --no-mixed --no-discordant --rna-strandness RF --no-unal --known-splicesite-infile splicesites.txt -x ~/erik_home/erik/genomes/hg19_ucsc_minus_hap -1 180403_NB501037_0252_AHLTNGBGX5/Fastq/KE210218-$1_R1_001.fastq.gz,180404_NB501037_0253_AHLTVTBGX5/Fastq/KE210218-$1_R1_001.fastq.gz -2 180403_NB501037_0252_AHLTNGBGX5/Fastq/KE210218-$1_R2_001.fastq.gz,180404_NB501037_0253_AHLTVTBGX5/Fastq/KE210218-$1_R2_001.fastq.gz | samtools view -bS - > bam/$1.bam) >& logs/$1.hisat_log
# echo $1: aligned and converted to bam
# 
# # Quantify tiles
# #(samtools view -b -q10 bam/$1.bam | ~/erik_home/erik/bin/BEDTools-Version-2.16.2/bin/coverageBed -split -counts -abam stdin -b hg18_hg19_1kb_tiles.bed | cut -f4,5 > results/$1.tile_counts) >& logs/$1.bedtools_log
# echo $1: tiles quantified
# 
# # Quantify exons                                                                                                   
# #(samtools view -b -q10 bam/$1.bam | ~/erik_home/erik/bin/BEDTools-Version-2.16.2/bin/coverageBed -split -counts -abam stdin -b gencode_exons.bed | cut -f4,5 > results/$1.exon_counts) >& /dev/null
# echo $1: exons quantified
# 
# # Quantify genes - q10 threshold is default with latest htseq-version (0.6.0 used here) (-s <yes/no/reverse> , whether the data is from a strand-specific assay (default: yes))
# (samtools view bam/$1.bam | htseq-count -m intersection-strict -s reverse - wgEncodeGencodeMerged_reformatted_LINC00942_fixed.gtf > results/$1.gene_counts) >& /dev/stdout | tail -n200 > logs/$1.htseq_log
# echo $1: genes quantified
# 
# #sort by coordinates for indexing
# samtools sort bam/$1.bam -o bam/$1.sorted.bam
# echo $1: sorted
# 
# samtools index bam/$1.sorted.bam
# echo $1: indexed
# 
# echo $1: processing complete
