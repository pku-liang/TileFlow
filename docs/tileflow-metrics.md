# Profile Metrics

- Basic:
`Cycle`: the total cycle
`Energy`: total energy consumption

- Data Movement profiling:
    - Key: `Level::Metric[::Tensor]`. MemScope is the name of the user specified storage/compute unit. 
    
        - Metric:
            - For compute level:
                - `Flop`: the total flops;
                - `Energy`: the compute energy;
            - For storage level:
                - `Read`: the total read from peer/child level;
                - `Update`: the write back from the peer/child level
                - `Fill`: the total fill from the parent level;
                - `Read|Update|Fill::t`: the metric for tensor t;
                - `SpatialUtil`: the max utilization of PE units.
                - `CapUtil`: the max utilized capcity X utilized PE / (total capcity X #PE)
                - `Energy`: the memory access energy;
                - `SlowDown` >= 1: the slowdown of this level compared to the child level. > 1 slowdown indicates this level is bottleneck compared to the child levels.
    - Value: `double` value.
- Legacy: Per tile utilization
    - Key: `ConstraintType::Level`. `ConstraintType` = [MEM|Spatial]
    - Value: `double` value indicating utilization.