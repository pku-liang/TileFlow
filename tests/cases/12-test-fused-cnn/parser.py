import re

def main(input_file, output_file):
    output = 'mapping:\n'
    with open(input_file, 'r') as f:
        for line in f.readlines():
            keywords = line.strip().split(' ')
            if len(keywords) == 0: break 
            keyword = keywords[0]
            keywords = keywords[1:]
            ntabs = 0
            while ntabs < len(line) and line[ntabs] in ' \t':
                ntabs += 1
            ntabs = ntabs // 2 + 1
            # print (len(re.findall('(\t*)[A-Za-z].*', line)[0]))
            if keyword == 'for' or keyword == 'pfor':
                varnames, bounds, tags = re.findall(r'p?for (.*) in [\(\[](.*)[\)\]]:\s*(#.*)\s*', line)[0]
                varnames = [x.strip().upper() for x in varnames.strip().split(',')]
                bounds = [x.strip() for x in bounds.strip().split(',')]
                assert(len(varnames) == len(bounds))
                tags = [x.strip() for x in tags.strip('# ').split(',')]
                output += '  '*(ntabs-1) + ('- ' if ntabs!= 1 else '  ') + 'node-type: Tile\n'
                output += '  '*ntabs + 'type: ' + ('Temporal' if keyword == 'for' else 'Spatial') + '\n'
                output += '  '*ntabs + 'factors: ' + ' '.join([f'{k}={v}' for k,v in zip(varnames, bounds)]) + '\n'
                output += '  '*ntabs + 'permutation: ' + ''.join(reversed(varnames)) + '\n'
                for tag in tags: 
                    k,v = [x.strip() for x in tag.split(':')]
                    output += '  '*ntabs +f'{k}: {v}\n'
                output += '\n'
                output += '  '*ntabs + 'subtree:\n'
            elif keyword == 'scope':
                t = re.findall(r'scope (\w+)', line)[0]
                output += '  '*(ntabs-1) + ('- ' if ntabs!= 1 else '  ') + 'node-type: Scope\n'
                output += '  '*ntabs + 'type: ' + t + '\n'
                output += '\n'
                output += '  '*ntabs + 'subtree:\n'
            elif keyword == 'endscope':
                pass 
            elif keyword == 'op':
                name = re.findall(r'op (\w+)', line)[0]
                output += '  '*(ntabs-1) + '- node-type: Op\n'
                output += '  '*ntabs + 'name: ' + name + '\n'
            else: 
                raise NotImplementedError
    with open(output_file, 'w') as f:
        f.write(output)
import sys

input_file = ''
output_file = 'map.yaml'

if len(sys.argv) > 1:
    input_file = sys.argv[1]
if len(sys.argv) > 2:
    output_file = sys.argv[2]

import os 

if os.path.isfile(input_file):
    main(input_file, output_file)