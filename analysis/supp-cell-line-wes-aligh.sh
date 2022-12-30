#!/bin/bash
# conda activate biosoft
# cd /data3/wsx_data/raw_cell_line/fq
# bash supp-cell-line-wes-aligh.sh /data3/wsx_data/raw_cell_line/fq
cd $1

cores=20

mkdir bam
sn=$(ls *.fastq.gz | tr ' ' '\n' | sed 's/_[12].fastq.gz//' | sort | uniq)
INDEX=/data1/database/human/hg38/bwa_index/hg38_p7

for id in ${sn}; do

  if [ ! -f bam/${id}.bam ]
  then
    if [ ! -f bam/${id}.sam ]
    then
      fq1=${id}_1.fastq.gz
      fq2=${id}_2.fastq.gz
      echo "Start aligning for ${id}"
      bwa mem -M -t $cores -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" ${INDEX} ${fq1} ${fq2} \
          > bam/${id}.sam 2>bam/${id}_bwa.log
    else
        echo "BWA align for ${id} is done before, directly go to sam > bam step"
    fi
  
    if [ $? -eq 0 ]
    then
        if [ ! -f bam/${id}.bam ]
        then
            samtools sort -@ $cores bam/${id}.sam -o bam/${id}.bam 2>bam/${id}_bam.log
            if [ $? -eq 0 ]
            then
                echo "Removing sam files"
                rm bam/${id}.sam
            else
                echo "Failed when using samtools sort, please check"
                exit 1
            fi
        fi
        echo "Done for ${id}." `date`
    else
        echo "Failed for ${id} in bwa." `date`
    fi
  
  fi

done
