dataflows = [
    'no-fuse',
    'chimera', 
    'pipeline', 
    'flash-attention', 
    'flat', 
    'tileflow' 
]

# dataflows = ['flash-attention'] #, 'flash-attention', 'flat', 'pipeline']

shapes = [
    # B H M N A L 
    ('BERT-small', [1,8,512,64,64,512]),
    ('BERT-base', [1,12,512,64,64,512]),
    ('BERT-large', [1,16,512,64,64,512]),
    ('ViT_14-Base', [1,12,256,64,64,256]),
    ('ViT_14-Large', [1,16,256,64,64,256]),
    ('ViT_14-High', [1,16,256,80,80,256]),
    ('ViT_16-Base', [1,12,192,64,64,192]),
    ('ViT_16-Large', [1,16,192,64,64,192]),
    ('ViT_16-High', [1,16,192,80,80,192]),
]

seqlen_shapes = [
    ('SeqLen-128', [1,8,128,64,64,128]),
    ('SeqLen-256', [1,8,256,64,64,256]),
    ('SeqLen-512', [1,8,512,64,64,512]),
    ('SeqLen-1024', [1,8,1024,64,64,1024]),
    ('SeqLen-2048', [1,8,2048,64,64,2048]),
]

# architectures = ['cloud-lowBW', 'cloud-midBW', 'cloud-highBW', 'cloud-largeBW']

import os 
import multiprocessing 
            
def worker(prefix, workload, shape, dataflow, objective, architecture):
    tag = f'{workload}-{dataflow}-{objective}-{architecture}'
    file_path = os.path.join(prefix, 'macro', f'{tag}.yaml') 
    output_path = os.path.join(prefix, 'log', f'{tag}.txt')
    if not os.path.isfile(file_path):
        B,H,M,N,A,L = shape 
        with open(file_path, 'w') as f:
            f.write(f'''
macro:
  B: {B}
  H: {H}
  M: {M}
  N: {N}
  A: {A}
  L: {L}        
output: {prefix}/out/{tag}
tileflow-mapper:
  timeout: 1200
  objective: {objective}
verbose: 1
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
    architectures,
    styles
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
        visualize(f'result/{prefix}', dataflows, shapes, objective, architectures, styles)

import pandas as pd 
import numpy as np

def gen_breakdown(df, prefix):
    
    os.system(f'mkdir -p {prefix}/breakdown')
    
    interested_metrics = [f'{level}::{metric}' for level in ['L0', 'L1', 'L2'] for metric in ['Read', 'Write', 'Update', 'Fill', 'SlowDown']]
    
    for (workload, arch), g in df.groupby(['workload', 'architecture']):
        g = g.drop(g[g.metric.map(lambda x: x not in interested_metrics)].index).reset_index()
        g = g.set_index(['metric', 'dataflow'])['value'].unstack(level=0)
        g['L2-L1'] = g['L2::Read'] + g['L2::Write']
        g['L1-RegFile'] = g['L1::Read'] + g['L1::Write']
        g['RegFile-ALU'] = g['L0::Read'] + g['L0::Write']
        g1 = g[['L2-L1', 'L1-RegFile', 'RegFile-ALU']]
        ax = g1.plot.bar(stacked=True)
        g2 = g[['L0::SlowDown', 'L1::SlowDown', 'L2::SlowDown']]
        g2['bottleNeck'] = pd.DataFrame(np.argmax(g2.to_numpy(), axis = -1), index = g2.index)
        g2['bottleNeck'].plot(ax = ax.twinx())
        ax.get_figure().savefig(f'{prefix}/breakdown/{workload}-{arch}.png', bbox_inches= 'tight')
        g1.to_csv(f'{prefix}/breakdown/{workload}-{arch}-breakdown.csv')
        g2.to_csv(f'{prefix}/breakdown/{workload}-{arch}-bottleneck.csv')
        

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
    
    df.to_csv(f'{prefix}/all.csv')
    
    gen_breakdown(df, prefix)        
    
    os.system(f'mkdir -p {prefix}/pics/{objective}')        
    
    def worker(group, key, other, metric, kind):
        ax = group.plot(kind = kind)
        os.system(f'mkdir -p {prefix}/pics/{objective}/{other}')
        os.system(f'mkdir -p {prefix}/pics/{objective}/{other}/{key}')
        ax.get_figure().savefig(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.png', bbox_inches='tight')
        ax.get_figure().savefig(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.pdf', bbox_inches='tight')
        group.to_csv(f'{prefix}/pics/{objective}/{other}/{key}/{metric}.csv')
        # group.to_csv(f'tmp/{name}.csv')
    
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

os.system('echo "" > error.sh')

if __name__ == "__main__":
    configs = [
        # ["Standard", ['cycle', 'energy'], dataflows, shapes, ['edge', 'cloud'], {'workload': 'bar', 'architecture': 'line'}],
        ["Sensitive-PE", ['cycle'], dataflows, shapes, ['cloud-4X4', 'cloud-12X14', 'cloud-16X16', 'cloud-32X32']
            , {'workload': 'bar', 'architecture': 'line'}],
        # ["SeqLen", ['cycle'], dataflows, seqlen_shapes, ['cloud-midBW'], {'workload': 'line', 'architecture': 'bar'}],
        # ["BandWidth", ['cycle'], dataflows, shapes, architectures, {'workload': 'line', 'architecture': 'line'}]
    ]
    procs = []
    for config in configs:
        p = multiprocessing.Process(target=main, args = config)
        procs.append(p)
        p.start()
    for p in procs: p.join()