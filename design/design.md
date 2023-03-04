We want to simulate general loopnest (mostly the imperfect loopnest) in tensor program on user-specified hardware. 

### The mapping IR design

The IR is represented in a tree. This tree is composed of three kinds of nodes:
- Scope Node: to specify the boundary of memory hierarchy;
- Tile Node: to specify the permutation;
- Op Node: to specify the arithmetic operations; 

#### Op Node 
Op nodes are leaves in the tree. They are used to specify the computation: 1. what tensors are accessed, 2. how the tensor is indexed. 

#### Tile Node 
Tile nodes are used to specify how loop is structured for a certain memory level. There are two types of Tile node, T-Tile (temporal), S-Tile (Spatial). T-tile tells that the computation is sequential, S-Tile tells that the computation is bond to spatial hardware (PE array, CPU cores, SMs, etc.). The key information carried by the a tile are: 
1. memory level: what level is this tile mapped to; 
2. loop transforamtion specifications: permutation/tiling factor; 
3. type: T/S-Tile; if is S-Tile, user may need to specify the binding information. 

#### Scope Node 
Scope Node is used to specify how `Tile Nodes` in a memory level are scheduled/executed in a memory level. 

There are 2 types of Scope Nodes, the P-Scope (Parallel Scope) and S-Scope (Sequential Scope). Children of P-Scope are either S-Scope or Tile node, who are executed in parallel and share the resource. Children of S-Scope are either P-Scope or Tile node, who are executed sequentially and can occupy full resource when executing. We require P-Scope to have more than 1 child to distinguish P-Scope from S-Scope for the degradation case of 1-child. 

#### An example


```
for m in [0, M):    
    for l in [0, L):
        fo k in [0, K):
            C[m,l] += A[m,k] * B[k,l]

for m in [0, M):
    for l in [0, L):
        C[m, l] = exp(C[m,l])

for m in [0, M):
    for l in [0, L):
        for n in [0, N):
            E[m,n] += C[m,l] * D[l,n]
```

```
for mo in [0, M, Tm):    
    for lo in [0, L, Tl):
        # S-Scope 
        {
            for mi in [0, Tm):
                for li in [0, Tl):
                    for k in [0, K, Tk):
                        for ki in [0, Tk):
                            C[m,l] += A[m,k] * B[k,l]
            for mi in [0, Tm): 
                for li in [0, Tn):
                    C[m, l] = exp(C[m,l])
            for mi in [0, Tm):
                for li in [0, Tl):
                    for ni in [0, Tn):
                        E[m,n] += C[m,l] * D[l,n]
        }
```

```
# S-Scope: MainMemory, Temporal   
{
    for mo in [0, M, Tm):    
        for lo in [Tl, L, Tl):
            # P-Scope: MainMemory, SyncAfterComputation = True  
            {
                # S-Scope: MainMemory
                {
                    for mm in [0, Tm, Ttm): # bind to spatial-X
                        for lm in [0, Tl, Ttl):  # bind to spatial-Y 
                            # S-Scope: PE 
                            {
                                for mi in [0, Ttm):
                                    for li in [0, Ttl):
                                        # Op node
                                        {
                                            l = lo + lm + li-Tl, m = mo + mm + mi
                                            C[m, l] = exp(C[m,l])
                                        }
                            }
                    for mm in [0, Tm, Ttm): # bind to spatial-X
                        for lm in [0, Tl, Ttl):  # bind to spatial-Y 
                            # S-Scope: PE 
                            {
                                for mi in [0, Ttm):
                                    for li in [0, Ttl):
                                        for ni in [0, Tn):
                                            # op node 
                                            {
                                                l = lo + li-Tl, m = mo + mi, n = no + ni 
                                                E[m,n] += C[m,l] * D[l,n]
                                            }
                            }
                }

                for mm in [0, Tm, Ttm): # bind to spatial-X
                    for lm in [0, Tl, Ttl):  # bind to spatial-Y
                        # S-Scope: PE temporal 
                        {
                            for mi in [0, Ttm):
                                for li in [0, Ttl):
                                    for k in [0, K):
                                        # Op-Node 
                                        {
                                            l = lo + li, m = mo + mi, k = ko + ki 
                                            C[m,l] += A[m,k] * B[k,l]
                                        }
                        }
            }
}
```