# MCTS
- monte caro tree search;
- problem description: 
    - A constraint optimization problem:
    $$
    Min_{variables} obj \\
    s.t.\ Constraints(variables) = 1 
    $$
    - variables:
        - tiling factors; 
        - scope type (not realized yet); 
        - permutation (not realized); 
    - constraints:
        - loop count constraint: example: $\Pi_{i} t^j_i \leq LoopCount_j$
        - memory constraint: example: $\Sigma_j\Pi_i t^k_{ij} \leq MemSize_k$
        - spatial constaint: example: $\Sigma_i Max_j...\Sigma_k<t^x_{ijk}, t^y_{ijk}> \leq <fanoutX, fanoutY>$ 
    - objective:
        - energy/latency;
        - this is calculated as a black box by simulation;
- Algorithm: for the tile size only
    - How the algorithm works can be seen in [this](https://hci.iwr.uni-heidelberg.de/system/files/private/downloads/297868474/report_robert-klassert.pdf)
    - Encoding in TileFlow:
        - State: 
            - A `symbol table` recording if a variable is fixed. If it is fixed, record the fixed value; else record the candidate values.
            - The candidate values is derived by the loop count constraint.
        - Action: 
            - Choose the next variable to be fixed;
                - Use heuristic to choose the variable with minimum feasible candidate values.
            - Decide the variable's value:  
                - Use the MCTS's UCB method to decide the value.
        - State Transition: Fix the variable with given value; update the symbol table using all constraints:
            - For loop count constraints, use it to give concrete candidate values easily. 
            - For other two type of constraints, use them as 0/1 pruning condition.
        - Terminate condition: all variables are fixed or no feasible solutions;
        - reward: energy/cycle from TileFlow's simulation

