cnvkit.py call -g female snu16.cnvkit.segment.cns --purity 1 -o snu16.cnvkit.segment.cns2
cnvkit.py call -g male pc3.cnvkit.segment.cns --purity 1 -o pc3.cnvkit.segment.cns2
# https://www.mskcc.org/research-advantage/support/technology/tangible-material/caki-2-human-renal-cell-line
cnvkit.py call -g male CAKI2.cnvkit.segment.cns --purity 1 -o CAKI2.cnvkit.segment.cns2
# https://en.wikipedia.org/wiki/K562_cells
cnvkit.py call -g female K562.cnvkit.segment.cns --purity 1 -o K562.cnvkit.segment.cns2
# https://en.wikipedia.org/wiki/T-47D
cnvkit.py call -g female T47D.cnvkit.segment.cns --purity 1 -o T47D.cnvkit.segment.cns2

# GBM39 was derived from the primary GBM tumor of a 51-year-old male prior to radiation, chemotherapy (Temozolomide), craniotomy, and Carmustine (BCNU) wafers (Gliadel)
# https://www.mdpi.com/2077-0383/8/9/1269/htm
cnvkit.py call -g male GBM39.cnvkit.segment.cns --purity 1 -o GBM39.cnvkit.segment.cns2

# https://en.wikipedia.org/wiki/NCI-60
# https://www.cellosaurus.org/CVCL_0459
cnvkit.py call -g male H460.cnvkit.segment.cns --purity 1 -o H460.cnvkit.segment.cns2

# HCC827
# https://www.cellosaurus.org/CVCL_2063
cnvkit.py call -g female HCC827_Erlotinib_removed.cnvkit.segment.cns --purity 1 -o HCC827_Erlotinib_removed.cnvkit.segment.cns2
cnvkit.py call -g female HCC827_naive.cnvkit.segment.cns --purity 1 -o HCC827_naive.cnvkit.segment.cns2
# Keep naive as HCC827 represent

# HK301
# https://academic.oup.com/neuro-oncology/article/18/10/1367/2222625#supplementary-data
cnvkit.py call -g male HK301.cnvkit.segment.cns --purity 1 -o HK301.cnvkit.segment.cns2
