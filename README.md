Ex Vivo internal API application for running InferCNV on Microsoft compute. Using this requires MSR credentials and Azure to be already setup.

## Installation and Usage
- Data folder should either be 1) a folder of h5ad's, or 2) a folder of sample folders each containing the necessary infercnv files for each sample. Tree structure shown below.
```
- 1) ├── data
     │    ├── sample1.h5ad
     │    ├── sample2.h5ad
     │    ├── sample3.h5ad
     │    └── ...
- 2) ├── data
     │    ├── sample1
     │         ├── singleCell.counts.matrix
     │         ├── cellAnnotations.txt
     │         ├── gene_ordering_file.txt
     │    ├── sample2
     │         ├── singleCell.counts.matrix
     │         ├── cellAnnotations.txt
     │         ├── gene_ordering_file.txt
     │    └── ...
```

- As shown in 2), please strictly follow the naming convention as outlined [here](https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV#infercnv-2-step-execution-overview) on the InferCNV wiki.
- I recommend running this in a [tmux](https://github.com/tmux/tmux/wiki) session for infercnv output to be automatically fetched to local. Please also monitor Azure UI for debugging.
  
```
# make sure you are setup with authentication
conda activate amlt10
az login

# install and run
git clone https://github.com/harrlol/InfercnvAzureAPI.git
cd InfercnvAzureAPI
bash ./docker-shell.sh \
  -I PATH_TO_DATA_FOLDER \
  -O PATH_TO_OUTPUT_FOLDER \
  # below is optional
  -N JOB_NAME \
  -P N_PARALLEL \
  -T N_THREADS \
  -R "NORMAL1 NORMAL2 NORMAL3" \
  -M "malignant"
```

## Considerations on Parallelization

- Each sample is allowed a default of 4 threads. When samples are abundant, set N_THREADS to a lower number so that N_THREADS * N_PARALLEL < 48.
- Set N_PARALLEL close to the number of samples in your run. A high N_PARALLEL significantly boosts efficiency; on the contrary high N_THREADS with low N_PARALLEL behaves similar to single core processing.
- For each job, try to maximize thread usage as much as possible.
