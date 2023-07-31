### Validation with timeloop

This folder reproduces the experiment in Fig.7 a/b

- Folder Description
    - `data.pkl`: RTL simluation result of different mappings of GEMM operation on systolic architecture.
    - `arch/`, `map/`, `prob/`: the architecture/mapping/workload description.
    - `sample_output`: the output figures shown in the paper.
- Run Script
```sh
python script.py
```