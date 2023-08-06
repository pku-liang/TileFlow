# Fusion Dataflow design

TileFlow is a framework focusing on modeling fusion dataflow. In this tutorial, we would delve into the design process of fusion dataflow on the workload of self-attention. 

An example of the fusion dataflow design is available in map/map.yaml. You might be scaerd by its length. But don't worry, we will guide you through it! 

It is always a good option to start with running the script.

```sh 
tileflow arch/arch.yaml prob/prob.yaml map/map.yaml macro/macro.yaml
```

You can observe TileFlow tuning the dataflow continuously. 

Next we introduce how we implement the fusion dataflow in TileFlow. In TileFlow, fusion is performed by combining the iteration space of multiple operators. To do that, we introduce scope node in the mapping tree. For this example, an `Sequential` scope is introduced into the mapping descrition to indicates sequential execution of two operators:
```sh 
  ... 
  - node-type: Scope 
    type: Sequential 

    subtree: 
        ... 
```
The scope node is inserted below the tile for Mainmemory mapping and above tiles for the Cache's mapping, indicating the tiles are fused at the cache level. 

