dataflows = [
    'naive',
    'isos', 
    'fused-layer',
    'tileflow',
]

# dataflows = ['flash-attention'] #, 'flash-attention', 'flat', 'pipeline']
'''
CC1 64 112 112 192 128
CC2 32 147 147 64 80
CC3 64 56 56 128 64
CC4 128 28 28 256 128
CC5 16 227 227 64 16
'''
shapes = [
    # B H W C L K R S U V
    ('CC1', [1,112,112,64,192,128,3,3,3,3]),
    ('CC2', [1,144,144,32,64,80,3,3,1,1]),
    ('CC3', [1,56,56,64,128,64,3,3,1,1]),
    ('CC4', [1,28,28,128,256,128,3,3,1,1]),
    ('CC5', [1,224,224,16,64,16,3,3,1,1])
]

architectures = ['cloud', 'edge']

import os 
import multiprocessing 
            
def worker(prefix, workload, shape, dataflow, objective, architecture):
    tag = f'{workload}-{dataflow}-{objective}-{architecture}'
    file_path = os.path.join(prefix, 'macro', f'{tag}.yaml') 
    output_path = os.path.join(prefix, 'log', f'{tag}.txt')
    if not os.path.isfile(file_path):
        B,H,W,C,L,K,R,S,U,V = shape 
        with open(file_path, 'w') as f:
            f.write(f'''
macro:
  B: {B}
  H: {H} 
  W: {W} 
  C: {C} 
  L: {L}
  K: {K}
  R: {R}
  S: {S}
  U: {U}
  V: {V}      
output: {prefix}/out/{tag}
tileflow-mapper:
  timeout: 600
  objective: {objective}
verbose: 1
check: 
    loopcount: False
''') 
            
    cmd = f"source ./run.sh {dataflow} {file_path} {output_path} {architecture}"
    if not os.path.isfile(f'{prefix}/out/{tag}.csv'):
        os.system(cmd)
    else: 
        print(cmd+ ' already executed')

def main(
    prefix,
    objectives,
    dataflows, 
    shapes, 
    architectures
    ):
    
    os.system(f'mkdir -p result')
    os.system(f'mkdir -p result/{prefix}')
    os.system(f'mkdir -p result/{prefix}/pics')
    os.system(f'mkdir -p result/{prefix}/out')
    os.system(f'mkdir -p result/{prefix}/log')
    os.system(f'mkdir -p result/{prefix}/macro')
    
    procs = []
    
    for workload, shape in shapes: 
        for dataflow in dataflows: 
            for objective in objectives: 
                for architecture in architectures: 
                    p = multiprocessing.Process(target = worker,
                        args = [
                            f'result/{prefix}',
                            workload, 
                            shape, 
                            dataflow, 
                            objective, 
                            architecture])
                    procs.append(p)
                    p.start()

    for p in procs:
        p.join()
        
    for objective in objectives:
        visualize(f'result/{prefix}', dataflows, shapes, objective, architectures)

import pandas as pd 

def visualize(prefix, dataflows, shapes, objective, architectures, styles = {'workload': 'bar', 'architecture': 'line'}):
    dfs = []
    keys = []
    for dataflow in dataflows:
        for workload, _ in shapes:
            for architecture in architectures:
                tag = f'{workload}-{dataflow}-{objective}-{architecture}'
                filename = f'{prefix}/out/{tag}.csv' 
                if not os.path.isfile(filename):
                    print (filename + ' does not exists')
                    continue 
                keys.append((workload, dataflow, architecture))
                dfs.append(pd.read_csv(filename))

    df = pd.concat(dfs, keys = keys, names = ['workload', 'dataflow', 'architecture'])
    
    
    os.system(f'mkdir -p {prefix}/pics/{objective}')        
    df.to_csv(f'{prefix}/all_data.csv')
    
    def worker(group, key, other, metric, kind):
        ax = group.plot(kind = kind)
        os.system(f'mkdir -p {prefix}/pics/{objective}/{other}')
        os.system(f'mkdir -p {prefix}/pics/{objective}/{other}/{key}')
        ax.get_figure().savefig(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.png', bbox_inches='tight')
        ax.get_figure().savefig(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.pdf', bbox_inches='tight')
        group.to_csv(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.csv')
    
    workloads = [x[0] for x in shapes]
    
    procs = []
    for x, other, order in [('architecture', 'workload', workloads), ('workload', 'architecture', architectures)]:
        for (key, metric), group in df.groupby([x, 'metric']):
            try: 
                group = group.drop('metric', axis = 1).droplevel(x)
                group = group.reset_index().set_index([other, 'dataflow'])['value']
                group = group.unstack()
                group = group.reindex(order)
                group = group[dataflows]
                proc = multiprocessing.Process(target=worker, args = [group, key, other, metric, styles[other]])
                proc.start()
                procs.append(proc)
            except: 
                print (f'error with {key} {metric}')
    for proc in procs: proc.join()

os.system('echo "#!/bin/bash" > error.sh')

if __name__ == "__main__":
    configs = [
        ["Standard", ['cycle', 'energy'], dataflows, shapes, ['cloud']],
    ]
    procs = []
    for config in configs:
        p = multiprocessing.Process(target=main, args = config)
        procs.append(p)
        p.start()
    for p in procs: p.join()