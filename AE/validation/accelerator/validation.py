import pickle as pkl
from functools import reduce
import multiprocessing
import time
import pandas as pd 
import tempfile
import os

# n_proc = multiprocessing.cpu_count()
n_proc = os.cpu_count()
# data = [
#     ((256, 256, 256, 128, 128, 128), 78283),
#     ((256, 384, 256, 128, 128, 128), 111271),
#     ((256, 128, 96, 64, 32, 48), 38206),
#     ((1024, 512, 256, 64, 128, 64), 793662),
#     ((64, 64, 1024, 32, 32, 256), 67070)
# ]

io_data = pd.read_csv('./data/io_data.csv')
datafile = './data/data.pkl'

with open(datafile, 'rb') as f:
    data = pkl.load(f)

def run_gemm(M,N,K,MO,NO,KO,SX,SY):
    config_file = 'config/config-gemm-' + '-'.join([str(_) for _ in [M,N,K,MO,NO,KO,SX,SY]]) + '.yaml'
    output_file = tempfile.mktemp(suffix='.csv')
    with open(config_file, 'w') as f:
        f.write(f'''
output: {output_file}
verbose: 0
macro: 
  M: {M}
  N: {N}
  K: {K} 
  
  KO: {KO}
  MO: {MO}
  NO: {NO}

  SX: {SX}
  SY: {SY}  
''')
    cmd = f'tileflow arch/arch.yaml map/map-gemm.yaml prob/prob-gemm.yaml {config_file}'
    print (cmd)
    os.system(cmd)
    cycle = float(pd.read_csv(output_file).set_index('metric').loc['Cycle'].values)
    return cycle

def run(data, id):
    M1,N1,K1,micro_M1,micro_N1,micro_K1,M2,N2,K2,micro_M2,micro_N2,micro_K2 = data
    assert(N1 == micro_N1)
    assert(M1 == M2)
    assert(M1 % micro_M1 == 0)
    assert(K1 % micro_K1 == 0) 
    assert(N2 % micro_N2 == 0)
    assert(micro_M2 == micro_M1)
    assert(micro_N2 == N2)
    assert(K2 == N1)
    M = M1 
    N = N2 
    K = K1 
    L = N1
    MO = M1 // micro_M1
    MM = micro_M1 // 16 
    KO = K1 // micro_K1 
    KM = micro_K1 // 16 
    LO = K2 // micro_K2 
    LM = micro_K2 // 16 
    NO = N2 // micro_N2 
    
    SX = SY = 16
    
    filename = 'config-' + '-'.join([str(_) for _ in data]) + '.yaml'
    with open(f"config/{filename}", 'w') as f:
        f.write(
f'''
output: /tmp/_tmp-{id}.csv
verbose: 0
macro: 
  M: {M}
  N: {N}
  K: {K} 
  L: {L}
  
  KO: {KO}
  KM: {KM}
  LO: {LO}
  LM: {LM}
  MO: {MO}
  MM: {MM}
  NO: {NO}

  SX: {SX}
  SY: {SY}   
'''
        )
    first_gemm_cycle = run_gemm(M, L, K, MO, LO, KO, SX, SY) # for graph-based
    second_gemm_cycle = run_gemm(M, N, L, MO, NO, LO, SX, SY) # for graph-based
    os.system(f'tileflow arch/arch.yaml map/map.yaml prob/prob.yaml config/{filename}')
    tileflow = pd.read_csv(f'/tmp/_tmp-{id}.csv')
    df = pd.DataFrame(data = [M, N, K, L, micro_M1, micro_K1, micro_K2] + list(tileflow.value) + [first_gemm_cycle, second_gemm_cycle],
                    index = ['M', 'N', 'K', 'L', 'micro_M1', 'micro_K1', 'micro_K2'] + list(tileflow.metric) + ['first_gemm_cycle', 'second_gemm_cycle'],
                    columns=['value'])
    df = df.T 
    
    return df
    
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

def analyze(df, metrics = 'Cycle', target = 'real'):
    model = LinearRegression()
    xs = df[metrics].values 
    if type(metrics) == str:
        metric = metrics  
        metrics = [metrics] 
        xs = xs.reshape(-1,1)
    else: 
        metric = ','.join(metrics)
    model.fit(xs, df[target])
    if len(metrics) > 1: df[metric] = model.predict(xs)
    print (f'relation: {metric} v.s. {target}')
    print ('\tscore:', model.score(xs, df[target]))
    print ('\tk: ', model.coef_)
    print ('\tb: ', model.intercept_)
    ave_err = np.mean(np.abs((df[metric] - df[target]) / df[target]))
    print ('\taverage error: ', ave_err)
    # df['tmp'] = df[metric]
    # ax = df.plot(x=metric, y = 'tmp', label= 'y=x', c = 'b')
    # df.plot.scatter(ax=ax, x=metric, y = target)
    # plt.legend()
    # ax.get_figure().savefig(f'result-linear-{metric}-{target}.png', bbox_inches = 'tight')
    
    df['ratio'] = df[metric] / df[target]
    df['dataflow'] = np.arange(len(df.index))
    
    _, ax = plt.subplots(1,1,figsize=(6,4.5))
    df.plot.scatter(ax = ax,x = 'dataflow', 
                         y = 'ratio', 
                        #  x_label = 'dataflows', 
                        #  y_label = 'TileFlow v.s. Real', 
                         title = metric, 
                         ylim=(0,1.3), 
                         xlim=(0,len(df.index)),
                         marker='^',
                         color = 'white',
                         edgecolors = 'blue',
                         zorder = 2,
                         grid = True,
                         s = 50
                        )
    if metric == 'Cycle':
        df['graph-ratio'] = df['graph-based'] / df[target] 
        df.plot.scatter(ax = ax,x = 'dataflow', 
                         y = 'graph-ratio', 
                        #  x_label = 'dataflows', 
                        #  y_label = 'TileFlow v.s. Real', 
                         title = metric, 
                         ylim=(0,2), 
                         xlim=(0,len(df.index)),
                         marker='o',
                         color = 'white',
                         edgecolors = 'orange',
                         zorder = 2,
                         grid = True,
                         s = 50
                        )
    ax.hlines(y=1, xmin=0, xmax = len(df.index), color='red', linestyle='--', zorder = 1, linewidth=3)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    plt.legend(fontsize="x-large")
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.xaxis.offsetText.set_fontsize(20)
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)
    ax.get_figure().savefig(f'result-absolute-{metric}-{target}.pdf', bbox_inches = 'tight')
    ax.get_figure().savefig(f'result-absolute-{metric}-{target}.png', bbox_inches = 'tight')
    
def do_work(procnum, return_dict):
    rets = []
    stride = (len(data) + n_proc - 1) // n_proc
    start = stride*procnum
    end = min(len(data), start + stride)
    for i in range(start, end):
        shape, cycle = data[i]
        io = io_data.iloc[i:i+1]
        ret = run(shape, procnum) 
        ret['cycle'] = cycle
        for key in io:
            ret[key] = io[key].iloc[0]
        rets.append(ret)
    
    return_dict[procnum] = rets


def main():
    os.system('mkdir -p config')
    if not os.path.isfile('out.csv'):
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        jobs = []
        for i in range(n_proc):
            p = multiprocessing.Process(target = do_work, args = (i, return_dict))
            jobs.append(p)
            p.start()
        for proc in jobs: proc.join()
        
        ret = pd.concat(reduce(lambda x, y: x+y, return_dict.values(), []))
         # Compute the accelerator's enenrgy;
        ret['flops'] = ret['N'] * ret['K'] * ret['L'] + ret['M'] * ret['N'] * ret['L']
        reg_energy = 249.643
        cache_energy = 13.8392
        memory_energy = 200
        energy_table = {
            'mac::Flops':1,
            'RegFile::Fill': reg_energy,
            'RegFile::Read': reg_energy,
            'RegFile::Update': reg_energy,
            'Cache::Read': cache_energy,
            'Cache::Update': cache_energy,
            'Cache::Fill': cache_energy,
            'MainMemory::Read': memory_energy,
            'MainMemory::Fill': memory_energy,
            'MainMemory::Update': memory_energy,
        }
        ret['estimated_energy'] = 0.0
        for k in energy_table:
            ret['estimated_energy'] += energy_table[k] * ret[k]
        
        # Compute the accelerator's enenrgy;
        energy_estimation_table = {
            'mem_to_buf': memory_energy / 16 + cache_energy / 16,
            'buf_to_mem': memory_energy / 16 + cache_energy / 16,
            'buf_to_reg': (cache_energy / 16 + reg_energy / (256 * 16)),
            'reg_to_buf': (cache_energy / 16 + reg_energy / (256 * 16)),
            'flops': 1 + reg_energy / 256 * 4
        }
        ret['energy'] = 0.0
        for k in energy_estimation_table:
            ret['energy'] += energy_estimation_table[k] * ret[k]
        ret['energy'] *= 1.21
        
        ret.to_csv('out.csv')
    else: 
        ret = pd.read_csv('out.csv')
        
    # Compute the estimation of graph-based approach
    ret['graph-based'] = ret['first_gemm_cycle'] + ret['second_gemm_cycle']
    ret['mid_tensor'] = ret['M'] * ret['L']
    ret['graph-based'] -= ret['mid_tensor'] / 4
    print ('n_data:', len(ret.index))
    analyze(ret, 'Cycle', 'cycle')
    analyze(ret, 'Energy', 'energy')
    return ret 


if __name__ == '__main__':
    # run_gemm(512, 512, 64, 4, 4, 2, 16, 16)
    with plt.style.context('seaborn-white'):
        main()