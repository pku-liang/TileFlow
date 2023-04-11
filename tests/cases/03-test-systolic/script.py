import pickle as pkl
from functools import reduce
import multiprocessing
import time

# n_proc = multiprocessing.cpu_count()
n_proc = 112
# data = [
#     ((256, 256, 256, 128, 128, 128), 78283),
#     ((256, 384, 256, 128, 128, 128), 111271),
#     ((256, 128, 96, 64, 32, 48), 38206),
#     ((1024, 512, 256, 64, 128, 64), 793662),
#     ((64, 64, 1024, 32, 32, 256), 67070)
# ]

with open("data.pkl", 'rb') as f:
    data = pkl.load(f)

import os 
import pandas as pd 

def run(data, id):
    M,N,K,MO,MM,NO,NI,KO,KM = data
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
    MO: {MO}
    MM: {MM}
    NO: {NO}
    NI: {NI}
    KO: {KO}
    KM: {KM}        
'''
        )
    os.system(f'tileflow arch/arch.yaml map/map.yaml prob/prob.yaml config/{filename}')
    tileflow = pd.read_csv(f'/tmp/_tmp-{id}.csv').set_index('metric').T
    # df = pd.DataFrame(data = [M, N, K, MO, NO, KO, tileflow['Cycle'].iloc[0], tileflow['Energy'].iloc[0]],
    #                   index = ['M', 'N', 'K', 'MO', 'NO', 'KO','tileflow-cycle', 'tileflow-ener'],
    #                   columns = ['value'])
    os.system(f'timeloop-model arch/arch.yaml map/map-timeloop.yaml prob/prob-timeloop.yaml config/{filename}')
    timeloop = pd.read_csv(f'/tmp/_tmp-{id}.csv').set_index('metric').T
    df = pd.DataFrame(data = [M, N, K, MO, NO, KO, timeloop['Cycle'].iloc[0], tileflow['Cycle'].iloc[0],timeloop['Energy'].iloc[0], tileflow['Energy'].iloc[0]],
                    index = ['M', 'N', 'K', 'MO', 'NO', 'KO','timeloop-cycle', 'tileflow-cycle', 'timeloop-ener', 'tileflow-ener'],
                    columns=['value'])
    df = df.T 
    return df
    
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

def analyze(df, metric = 'tileflow-cycle', target = 'real'):
    model = LinearRegression()
    model.fit(df[metric].values.reshape(-1,1), df[target])
    print (f'relation: {metric} v.s. {target}')
    print ('\tscore:', model.score(df[metric].values.reshape(-1,1), df[target]))
    print ('\tk: ', model.coef_)
    print ('\tb: ', model.intercept_)
    ave_err = np.mean(np.abs((df[metric] - df[target]) / df[target]))
    print ('\taverage error: ', ave_err)
    df['pred'] = model.predict(df[metric].values.reshape(-1,1))
    ax = df.plot.scatter(x=metric, y = target, label=target)
    df.plot.line(ax = ax, x=metric, y = 'pred', label='pred', c = 'black')
    plt.legend()
    ax.get_figure().savefig(f'result-{metric}-{target}.png')
    
def do_work(procnum, return_dict):
    rets = []
    stride = (len(data) + n_proc - 1) // n_proc
    start = stride*procnum
    end = min(len(data), start + stride)
    for (M,N,K,micro_M,micro_N,micro_K), cycle in data[start:end]:
        MO = M // micro_M 
        MM = micro_M // 16 
        NO = N // micro_N
        NI = micro_N
        KO = K // micro_K 
        KM = micro_K // 16 
        
        ret = run((M,N,K,MO,MM,NO,NI,KO,KM), procnum) 
        ret['real'] = cycle   
        rets.append(ret)
    
    return_dict[procnum] = rets


def main():
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for i in range(n_proc):
        p = multiprocessing.Process(target = do_work, args = (i, return_dict))
        jobs.append(p)
        p.start()
    for proc in jobs: proc.join()
    
    print (return_dict.values())
    ret = pd.concat(reduce(lambda x, y: x+y, return_dict.values(), []))
    print(ret)
    # ret.to_csv('out.csv')
    analyze(ret, 'tileflow-cycle', 'timeloop-cycle')
    analyze(ret, 'tileflow-cycle', 'real')
    analyze(ret, 'timeloop-cycle', 'real')
    analyze(ret, 'timeloop-ener', 'tileflow-ener')
    return ret 


if __name__ == '__main__':
    main()