### Validation with accelerator

This folder produces the result of Fig.7 c/d. 

- Folder Description
    - `data/`: RTL simulation result of a systolic based hardware.
    - `prob/`,`map/`,`arch`: the workload, mapping, architecture descriptions.
    - `sample_outputs`: the Fig.7 c/d in the paper.

- Run Scripts
```sh
python ./validation.py
```

- Run time: 5 seconds on 112 cores CPU. 