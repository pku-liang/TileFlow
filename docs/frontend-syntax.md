# Frontend Syntax of TileFlow 

TileFlow uses yaml for input configuration. There are 3 required fields for an  application: `architecture` field for archtecture specification, `problem` field for problem specification, and `mapping` field problem-architecture mapping specification. Normally, we implement three fields in 3 separate `yaml` files. Please see `tests/cases` for examples.   

## Arch File 

The `architecture` field uses the syntax of `Timeloop`. The memory hierarchy is described in a recursive style from the outermost level (Main Memory) to the innermost level (ALU). Every level has following subfields:
- name (required) -> string: a unique tag for the level. If there are multiple instances, uses NAME[0..N] in the name field. For example, PE[0..15] is used to describe a memory level with 16 PE cores.
- class (optional) -> string: [DRAM|SRAM|regfile]
- local: the local memory level;
- subtree: a list of sub memory/ALU units. 
- attributes: describe memory-related features. Useful ones (details in `3rdparty/timeloop/src/model/buffer.cpp::model::ParseSpecs`):
    - width: the number of bits in memory; (width = depth * block-size * word-bits)
    - depth: the number of blocks in a memory;
    - block-size: the number of word in a block;
    - word-bits: the word size in bits;
    - read_bandwidth: word/cycle?
    - write_bandwidth: word/cycle?

## Prob File 

The `problem` field describes the workload as a function. There are some required sub-filedsL:
- `io`: the input and output of the function.
- `dimensions`: all the dimensions appear in describing the tensors. 
- `instance`: the problem size for each dimension. 
- `ops`: a list of operations.

For `ops` field, the user need to specify the unique `name` for every operation, the `dimensions` used for every operation, and the access patterns `data-spaces` for every tensors. For the data-space, the user specify the tensor name in `name` and the access pattern in `projection`. For `projection`, the syntax is Product of Sum of Product. For example, T[A+B*2][C] is described as 

```
projection: 
    - [[A], [B,2]]
    - [C]
```

 A read-write tensor is marked manually by a flag `read-write`. Also, for each `op` the user needs to specify the `ins` field for input tensor and `out` field for a single output tensor. This information is used for def-use analysis. 

## Mapping File 

The `mapping` field decribed the mapping in a hierarchical tree. The children of every node in the tree is a list under the subfield `subtree`. There are three kinds of nodes in the tree (specified in `node-type`): 

- Scope Node: to specify the boundary of memory hierarchy; The only attribute of a scope node is its sub-types: Sharing/Temporal/Spatial/Pipeline (specified in `type`).  

- Tile Node: to specify the temporal/spatial mapping of loops. Attributes:
    - type: spatial/temporal;
    - permutation: A permutation of iterations; from inner to outer. 
    - target: A memory level name; 
    - factors: dimension name -> int. 

- Op Node: to specify the arithmetic operations; Attributes:
    - name: the name of operation;
    - binding: iteration-name-in-loop -> dimension-name-in-op 




