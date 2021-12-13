### Data for SMM Clustering

Author: Shankara Anand (sanand@broadinstitute.org)

**Reference Dir** (`./ref`)
* `./ref/gmts/`
  * GMTs from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/) - need to be downloaded
  * `./ref/gmts/mm_sigs_custom.gmt.txt`: custom GMT created from multiple myeloma signatures (see `../../Fig2/1_process_rna.ipynb`)
  * `./ref/gmts/staudt_2020.gmt.txt`: custom GMT created from SignatureDB file from [Staudt](https://lymphochip.nih.gov/signaturedb/)
* `./ref/de_highlights/`: genes to highlight in Fig2 in clusters C3 (FMD) & C4 (CND)
* `./ref/sonneveld_*_sigs.txt`
  * Pulled from [here](https://ashpublications.org/blood/article/116/14/2543/27550/Gene-expression-profiling-for-molecular)
  * (see `../../Fig2/1_process_rna.ipynb`)
*  `./ref/zhan_*_sigs.txt`
   * Pulled from [here](https://ashpublications.org/blood/article/108/6/2020/22665)
   * (see `../../Fig2/1_process_rna.ipynb`)
*  `./ref/SignatureDB_030920.txt`
   * Pulled from [here](https://lymphochip.nih.gov/signaturedb/)

**Raw Data Dir** (`./raw`)
* Binarized file final study list for DNA features
* Annotations/groupings may be found `../funcs/features.py`

**RNA Data Dir** (`./rna`)
* Output files from RNA-SeQC (in progress, please contact sanand@broadinstitute.org)
* Processing details may be found here: https://github.com/broadinstitute/gtex-pipeline
