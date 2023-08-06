# Tutorial on GEMM 
In `TileFlow`, users specifies the architecture, mapping, and problem description in `TileFlow` 's [frontend syntax](../../docs/frontend-syntax.md). 
This folder uses the general matrix multiply (GEMM) to demonstrate `TileFlow`'s workflow. 

In this example, we described a spatial [accelerator](arch/arch.yaml) with three memory levels. Further, we described the computation in `prob/prob.yaml`, and the dataflow (mapping) in `map/map.yaml`. To instantiate the mapping problem, we wrote a macro file (marco.yaml) to instantiate the shape of the problem.

To run `TileFlow`, please ensure the binary is in your system's path, and simply append all configuration files as paramers (order-blind):

```sh
tileflow arch/arch.yaml prob/prob.yaml map/map.yaml macro.yaml 
```

> Tips: All configuration file can be combined as one file.

In less than a second, TileFlow will output the currently optimal dataflow found, along with the profiling metrics, including latency, energy, data movement volume, etc. Example output is shown in `gemm.csv`, `gemm.mapping.txt`.

Next, we will illustrate how a dataflow is described in TileFlow, i.e. the mapping description. In TileFlow, we design dataflow by mapping each operator on the software side to each memory level of the hardware. For the matrix multiply example, the operator's mapping is represented in a chain of tile nodes, where every tile node describe the mapping of a memory level. For example, we map the computation of the MainMemory Level using a temporal tile:

```
    node-type: Tile 
    type: temporal 
    factors: M = MO N = NO K= KO
    permutation: KMN 
    target: MainMemory 
```, 

And map the PE arrays using a spatial Tile:

```
    node-type: Tile 
    type: spatial  
    factors: M=16 K=16
    permutation: MK
    split: 1
    target: Cache
    multicast: true
```
> `multicast` metric indicates the hardware is able to perform multicast. 

You are free to change the tiling factors and permutations of the tile node. Or, you can replace concrete tile sizes with unspecified macros. These macros can be automatically decided by TileFlow.