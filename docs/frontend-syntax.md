# Frontend Syntax of TileFlow 

TileFlow uses yaml for input configuration. There are 3 required fields for an  application: `architecture` field for archtecture specification, `problem` field for problem specification, and `mapping` field problem-architecture mapping specification. Normally, we implement three fields in 3 separate `yaml` files. Please see `tests/cases` for examples.   

## Arch File 

The `architecture` field uses the syntax of `Timeloop`, see [this](https://timeloop.csail.mit.edu/timeloop/input-formats/design/architecture) for description. 

## Prob File 

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

## Mapping File 

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

- Op Node: to specify the arithmetic operations; Attributes:
    - name: the name of operation;
    - binding: iteration-name-in-loop -> dimension-name-in-op 




