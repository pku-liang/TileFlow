# Frontend Syntax of TileFlow 

TileFlow uses yaml for input configuration. There are 3 required fields for an  application: `architecture` field for archtecture specification, `problem` field for problem specification, and `mapping` field problem-architecture mapping specification. Normally, we implement three fields in 3 separate `yaml` files. A typical config file is like: 

```
problem:
    ...
architecture:
    ...
mapping: 
    ...
check (optional): 
    ...
tileflow-mapper (optional):
    ...
macro (optional): 
    ...
output (optional): 
verbose (optional):  
```

Please see `tests/cases` for examples.

## Arch Scope 

The `architecture` field uses the syntax of `Timeloop`, see [this](https://timeloop.csail.mit.edu/timeloop/input-formats/design/architecture) for description. 

## Prob Scope 

The `problem` field extends `Timeloop`'s syntax to support multi-op. The file is organized like:
```
problem:
    io:
        ...
    dimensions:
        ...
    instance:
        ...
    ops:    
    - TIMELOOP-OP1
    - TIMELOOP-OP2
    ... 
```

- `io`: the input and output of the function.
- `dimensions`: all the dimensions appeared in describing the tensors. 
- `instance`: the specification for parameters, see [this](https://timeloop.csail.mit.edu/timeloop/input-formats/problem#problem-shape) for description. 
- `ops`: a list of tensor operations, follow the same syntax with [shape](https://timeloop.csail.mit.edu/timeloop/input-formats/problem#problem-shape) in timeloop without the instance field. An extra field of each op is the `ins` and `out` field to specify the IO of each operation.

## Mapping Scope 

The `mapping` field decribed the mapping in a tree. A Node in a tree is like:

```yaml
node-type: TILE|Scope|Op
# optional attributes
type: [Sharing|Temporal|Spatial|Pipeline|temporal|spatial]
factors:
permutation:
target:
split: 

subtree:   
- CHILD1
- CHILD2
...
```

There are three kinds of nodes:

- Scope Node: to specify the boundary of memory hierarchy; The only attribute of a scope node is its sub-types: Sharing/Temporal/Spatial/Pipeline.  

- Tile Node: to specify the temporal/spatial mapping of loops. The key attributes include factors, permutations, target, split, etc.. See [this](https://timeloop.csail.mit.edu/timeloop/input-formats/mapping) for illustration.
    - Key knobs:
        - multicast [true|false]: used for spatial tile to specify whether the higher memory level's bandwidth can perform multicast. 

- Op Node: to specify the arithmetic operations; Attributes:
    - name: the name of operation;

- To enable `tileflow-mapper`, user can simply replace the number in the specification for tile factors with arbitrary string. 

## Check Scope:
- To cutomize different kinds of checking;
- Attributes: 
    - `mem`(bool): whether or not to enable memory capticy check;
    - `loopcount`(bool): whether or not to enable the loopcount check (whether the multiplication of tile factors equal the shape);
    - `spatial`(bool): whether the spatial core usage is exceeded.

## Mapper Scope:
- Specify the configuration for mapper
- Attributes: 
    - `alg`[random, mtcs]: the searching algorithm for mapper;
    - `timeout`[INT]: the searching timeout in seconds;
    - `topk`(unsigned): record topK candidates.

## Macro Scope: 
    - key(string): value(int) pairs of macros. The macros can be used for instanciation of `factor` scope of `tile` nodes, and the instanciation of `instance scope` of `problem scope`. 

## Others

- `macro` attribute: list some constant values that can be used as tile factors/tensor shapes
- `verbose` attribute: specify the verbose level;
- `output`: the prefix for output; including 1. `$(output).csv` for cycle/energy/profiling results; 2. `$(output).mapping.csv` for searched best dataflow; 3. `$(output).tuning.csv` for mapper tunning log; 

