import pickle as pkl 
import re

shape = [512,512,64,64]

def parse_raw_data(raw_data):
    ret = []
    for line in raw_data.split('\n'):
        key_tuple = shape.copy()
        failed = False
        for key in ['MO', 'MM', 'KO', 'KM', 'LO', 'LM']:
            pattern = f'<{key},1,(\d+)>'
            res = re.findall(pattern, line)
            if len(res) == 0:
                print (line, key)
                failed = True
                break 
            key_tuple.append(int(res[0]))
        if failed: continue
        res = re.findall('value: (\d+)', line)
        assert len(res) 
        ret.append((tuple(key_tuple), int(res[0])))
    return ret 
     
def get_dataset():
    filename = 'data/No_Softmax/data.pkl'
    ret = dict()
    with open(filename, 'rb') as f: 
        data = pkl.load(f) 
        for line in data: 
            (M1,N1,K1,micro_M1,micro_N1,micro_K1,M2,N2,K2,micro_M2,micro_N2,micro_K2), value = line 
            MO = M1 // micro_M1
            MM = micro_M1 // 16 
            KO = K1 // micro_K1 
            KM = micro_K1 // 16 
            LO = K2 // micro_K2 
            LM = micro_K2 // 16 
            ret[(M1, N1, K1, N2, MO, MM, KO, KM, LO, LM)] = value 
    return ret 

def reconvert(shape):
    M1, N1, K1, N2, MO, MM, KO, KM, LO, LM = shape 
    micro_M1 = 16 * MM
    assert MO == (M1 // micro_M1)
    micro_N1 = N1 
    micro_K1 = 16 * KM 
    M2 = M1 
    K2 = N1 
    micro_M2 = micro_M1 
    micro_N2 = N2 
    micro_K2 = 16 * LM 
    assert LO == (K2 // micro_K2)
    return (M1,N1,K1,micro_M1,micro_N1,micro_K1,M2,N2,K2,micro_M2,micro_N2,micro_K2)

raw_data = '''
0:	<MO,1,4>,<KO,1,4>,<LO,1,1>,<KM,1,1>,<LI,1,512>,<MM,1,8>,<NO,1,1>,<LM,1,32>,<NI,1,64>,, value: 140288
1:	<MO,1,4>,<KO,1,1>,<LO,1,1>,<KM,1,4>,<LI,1,512>,<MM,1,8>,<NO,1,1>,<LM,1,32>,<NI,1,64>,, value: 140288
2:	<MO,1,4>,<KO,1,2>,<LO,1,1>,<KM,1,2>,<LI,1,512>,<MM,1,8>,<NO,1,1>,<LM,1,32>,<NI,1,64>,, value: 140288
3:	<MO,1,4>,<KO,1,4>,<LO,1,2>,<KM,1,1>,<LI,1,256>,<MM,1,8>,<NO,1,1>,<LM,1,16>,<NI,1,64>,, value: 141312
4:	<MO,1,4>,<KO,1,1>,<LO,1,2>,<KM,1,4>,<LI,1,256>,<MM,1,8>,<NO,1,1>,<LM,1,16>,<NI,1,64>,, value: 141312
5:	<MO,1,4>,<KO,1,2>,<LO,1,2>,<KM,1,2>,<LI,1,256>,<MM,1,8>,<NO,1,1>,<LM,1,16>,<NI,1,64>,, value: 141312
6:	<MO,1,4>,<KO,1,1>,<LO,1,4>,<KM,1,4>,<LI,1,128>,<MM,1,8>,<NO,1,1>,<LM,1,8>,<NI,1,64>,, value: 143360
7:	<MO,1,4>,<KO,1,4>,<LO,1,4>,<KM,1,1>,<LI,1,128>,<MM,1,8>,<NO,1,1>,<LM,1,8>,<NI,1,64>,, value: 143360
8:	<MO,1,4>,<KO,1,2>,<LO,1,4>,<KM,1,2>,<LI,1,128>,<MM,1,8>,<NO,1,1>,<LM,1,8>,<NI,1,64>,, value: 143362
9:	<MO,1,4>,<KO,1,2>,<LO,1,8>,<KM,1,2>,<LI,1,64>,<MM,1,8>,<NO,1,1>,<LM,1,4>,<NI,1,64>,, value: 147456
10:	<MO,1,4>,<KO,1,4>,<LO,1,8>,<KM,1,1>,<LI,1,64>,<MM,1,8>,<NO,1,1>,<LM,1,4>,<NI,1,64>,, value: 147456
11:	<MO,1,4>,<KO,1,1>,<LO,1,8>,<KM,1,4>,<LI,1,64>,<MM,1,8>,<NO,1,1>,<LM,1,4>,<NI,1,64>,, value: 147456
12:	<MO,1,4>,<KO,1,1>,<LO,1,1>,<KM,1,4>,<LI,1,512>,<MM,1,8>,<NO,1,2>,<LM,1,32>,<NI,1,32>,, value: 148480
13:	<MO,1,4>,<KO,1,2>,<LO,1,1>,<KM,1,2>,<LI,1,512>,<MM,1,8>,<NO,1,2>,<LM,1,32>,<NI,1,32>,, value: 148480
14:	<MO,1,4>,<KO,1,4>,<LO,1,1>,<KM,1,1>,<LI,1,512>,<MM,1,8>,<NO,1,2>,<LM,1,32>,<NI,1,32>,, value: 148480
15:	<MO,1,4>,<KO,1,4>,<LO,1,2>,<KM,1,1>,<LI,1,256>,<MM,1,8>,<NO,1,2>,<LM,1,16>,<NI,1,32>,, value: 149504
16:	<MO,1,4>,<KO,1,1>,<LO,1,2>,<KM,1,4>,<LI,1,256>,<MM,1,8>,<NO,1,2>,<LM,1,16>,<NI,1,32>,, value: 149504
17:	<MO,1,4>,<KO,1,2>,<LO,1,2>,<KM,1,2>,<LI,1,256>,<MM,1,8>,<NO,1,2>,<LM,1,16>,<NI,1,32>,, value: 149504
18:	<MO,1,4>,<KO,1,4>,<LO,1,4>,<KM,1,1>,<LI,1,128>,<MM,1,8>,<NO,1,2>,<LM,1,8>,<NI,1,32>,, value: 151552
19:	<MO,1,4>,<KO,1,1>,<LO,1,4>,<KM,1,4>,<LI,1,128>,<MM,1,8>,<NO,1,2>,<LM,1,8>,<NI,1,32>,, value: 151552'''

topk_data = parse_raw_data(raw_data)
dataset = get_dataset()
converted_dataset = []

for k, v in topk_data: 
    if k not in dataset: 
        print (k, 'is not in dataset')
        converted_dataset.append(list(reconvert(k)) + [v])
        continue 
    print (k, v, dataset[k])

import pandas as pd 
df = pd.DataFrame(data = converted_dataset, columns=
                  ['M1','N1','K1','micro_M1','micro_N1','micro_K1','M2','N2','K2','micro_M2','micro_N2','micro_K2','latency'])
df.to_csv('topk.csv')
# print (dataset)