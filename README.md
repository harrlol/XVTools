Project Ex Vivo internal application for running various bioinformatic utilities on Microsoft compute. Set up your amulet project and workspace as guided [here](https://deep-acapella-f1a.notion.site/Guide-to-your-sandbox-201d503097798026848be26e307b619e) prior to using this tool. 

## Installation and Usage
```
git clone https://github.com/harrlol/XVTools.git && cd XVTools && pip install .
conda activate amlt10
az login
```
Download data from Terra to Microsoft sandbox
```
xvtools download terra \
  --secret sa.json \
  --gcs-dir gs://<bucket>/PRISM_pilot/analysis/ \
  --dest ./data/patient123 \
  --pattern '*filtered*matrix*' \
  --exclude '.bam,.bai,.fastq,.fastq.gz,.fq,.fq.gz'
```
Submit infercnv job to azure
```
xvtools submit infercnv-aml \
  --data ./data/patient123 \
  --out ./out/patient123 \
  --ref-group-names-str 'T_Cell, Macrophage, B_Cell' \
  --n-parallel 4 --n-threads 2 --sku 8C15
```

## Notes
- Data folder should either be 1) a folder of h5ad's, or 2) a folder of sample folders each containing the necessary infercnv files for each sample. Tree structure shown below.
```
1) ├── data/
   │    ├── sample1.h5ad
   │    ├── sample2.h5ad
   │    ├── sample3.h5ad
   │    └── ...

2) ├── data/
   │    ├── sample1/
   │         ├── singleCell.counts.matrix
   │         ├── cellAnnotations.txt
   │         └── gene_ordering_file.txt
   │    ├── sample2/
   │         ├── singleCell.counts.matrix
   │         ├── cellAnnotations.txt
   │         └── gene_ordering_file.txt
   │    ├── sample3/
   │         ├── ...
   │    └── ...
```

- As shown in 2), please strictly follow the naming convention as outlined [here](https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV#infercnv-2-step-execution-overview) on the InferCNV wiki.
- I recommend running this in a [tmux](https://github.com/tmux/tmux/wiki) session for infercnv output to be automatically fetched to local. Please also monitor Azure UI for debugging.
  
## Considerations on Parallelization

- Each sample is allowed a default of 4 threads. When samples are abundant, set N_THREADS to a lower number so that N_THREADS * N_PARALLEL < 48.
- Set N_PARALLEL close to the number of samples in your run. A high N_PARALLEL significantly boosts efficiency; on the contrary high N_THREADS with low N_PARALLEL behaves similar to single core processing.
- For each job, try to maximize thread usage as much as possible.
