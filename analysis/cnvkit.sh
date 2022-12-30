# https://cnvkit.readthedocs.io/en/stable/quickstart.html#install-cnvkit
dir=tumor_only
mkdir /data3/wsx_data/cnvkit_result/$dir
cnvkit.py batch /data3/wsx_data/raw_cell_line/own_data/bam/*.bam \
  -p 10 \
  -n -t /home/wangsx/data/S07604514_hs_hg38/S07604514_Regions.bed -f /data1/database/human/hg38/genome.fa \
  --access /home/wangsx/data/access.hg38.bed --annotate /home/wangsx/data/refFlat.txt \
  -d /data3/wsx_data/cnvkit_result/$dir

# * Use this result
dir=tumor_vs_sampleB
mkdir /data3/wsx_data/cnvkit_result/$dir
cnvkit.py batch /data3/wsx_data/raw_cell_line/own_data/bam/*.bam \
  -p 20 \
  -n /data3/wsx_data/raw_cell_line/fq/bam/SRR13084972.bam /data3/wsx_data/raw_cell_line/fq/bam/SRR13084973.bam \
  -t /home/wangsx/data/S07604514_hs_hg38/S07604514_Regions.bed -f /data1/database/human/hg38/genome.fa \
  --access /home/wangsx/data/access.hg38.bed --annotate /home/wangsx/data/refFlat.txt \
  -d /data3/wsx_data/cnvkit_result/$dir
  
# cnvkit.py call pc3-1.cns --purity 1 -o pc3-1.call.cns2 -y -g male
# cnvkit.py call SNU16-1.cns --purity 1 -o SNU16-1.call.cns2 -y -g female

dir=tumor_vs_HCCBL
mkdir /data3/wsx_data/cnvkit_result/$dir
cnvkit.py batch /data3/wsx_data/raw_cell_line/own_data/bam/*.bam \
  -p 20 \
  -n /data3/wsx_data/raw_cell_line/fq/bam/SRR925779.bam \
  -t /home/wangsx/data/S07604514_hs_hg38/S07604514_Regions.bed -f /data1/database/human/hg38/genome.fa \
  --access /home/wangsx/data/access.hg38.bed --annotate /home/wangsx/data/refFlat.txt \
  -d /data3/wsx_data/cnvkit_result/$dir
