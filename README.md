Ex Vivo internal API application for running InferCNV on Microsoft compute

## Installation and Usage

I recommend running this in a tmux session for infercnv output to be automatically fetched to local. Please also monitor Azure UI for debugging.
```
git clone https://github.com/harrlol/InfercnvAzureAPI.git
cd InfercnvAzureAPI
sh ./docker-shell.sh \
  -I PATH_TO_DATA_FOLDER \
  -O PATH_TO_OUTPUT_FOLDER
```

## Considerations on Parallelization

- Each sample is allowed a default of 4 threads. When samples are abundant, set N_THREADS to a lower number so that N_THREADS * N_PARALLEL < 48.
- Set N_PARALLEL close to the number of samples in your run. A high N_PARALLEL significantly boosts efficiency; on the contrary high N_THREADS with low N_PARALLEL behaves similar to single core processing.
- For each job, try to maximize thread usage as much as possible.
