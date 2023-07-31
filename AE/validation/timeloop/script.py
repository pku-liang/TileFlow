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

def analyze(df, metric = 'tileflow-cycle', target = 'real', ylabel = None, title = "", xy = ["x", "y"]):
    if ylabel is None: 
        ylabel = metric + ' v.s. ' + target
    model = LinearRegression()
    model.fit(df[metric].values.reshape(-1,1), df[target])
    print (f'relation: {metric} v.s. {target}')
    print ('\tscore:', model.score(df[metric].values.reshape(-1,1), df[target]))
    print ('\tk: ', model.coef_)
    print ('\tb: ', model.intercept_)
    ave_err = np.mean(np.abs((df[metric] - df[target]) / df[target]))
    print ('\taverage error: ', ave_err)
    df['pred'] = df[metric]
    _, ax = plt.subplots(1,1,figsize=(6,4.5))
    df.plot.line(ax = ax, x=metric, y = 'pred', label='y=x', c = 'red', zorder = 0)
    df.plot.scatter(
        ax = ax, 
        x=metric, 
        y = target, 
        s = 50,
        title = title,
        color = 'white',
        edgecolor = 'blue',
        marker = '^',
        zorder = 1)
    plt.legend(fontsize=20)
    ax.set_xlabel(xy[0])
    ax.set_ylabel(xy[1])
    # ax.set_xlabel("")
    # ax.set_ylabel("")
    ax.set_title("")
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.xaxis.offsetText.set_fontsize(20)
    ax.get_figure().savefig(f'result-{metric}-{target}.pdf', bbox_inches = 'tight')
    ax.get_figure().savefig(f'result-{metric}-{target}.png', bbox_inches = 'tight')
    plt.close()
    # df['ratio'] = df[metric] / df[target]
    # df['dataflow'] = np.arange(len(df.index))
    # ax = df.plot.scatter(x = 'dataflow', 
    #                      y = 'ratio', 
    #                      s = 100,
    #                      xlabel = 'dataflows', 
    #                      ylabel = ylabel, 
    #                      ylim=(0,1.5), xlim=(0,131), marker='d',
    #                      color = 'pink',
    #                      edgecolors='black'
    #                      )
    # ax.hlines(y=1, xmin=0, xmax = len(df.index), color='red', linestyle='--')
    # plt.legend()
    # ax.get_figure().savefig(f'result-absolute-{metric}-{target}.png', bbox_inches = 'tight')
    
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
    os.system('mkdir -p config')
    if os.path.isfile('out.csv'):
        ret = pd.read_csv('out.csv')
    else:
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        jobs = []
        for i in range(n_proc):
            p = multiprocessing.Process(target = do_work, args = (i, return_dict))
            jobs.append(p)
            p.start()
        for proc in jobs: proc.join()
        
        ret = pd.concat(reduce(lambda x, y: x+y, return_dict.values(), []))
        ret.to_csv('out.csv')
    print ('n_data', len(ret.index))
    analyze(ret, 'tileflow-cycle', 'timeloop-cycle', 'TileFlow v.s. TimeLoop', 'Cycle', ['TimeLoop', 'TileFlow'])
    # analyze(ret, 'tileflow-cycle', 'real')
    # analyze(ret, 'timeloop-cycle', 'real')
    analyze(ret, 'timeloop-ener', 'tileflow-ener', 'TileFlow v.s. TimeLoop', 'Energy', ['TimeLoop', 'TileFlow'])
    return ret 


if __name__ == '__main__':
    with plt.style.context('seaborn-white'):
        main()