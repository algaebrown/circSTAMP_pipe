import re
import sys
from optparse import OptionParser
from pathlib import Path
def parse_skewer_log(f):
    data = {}
    with open(f) as fhandle:
        lines = '\n'.join(fhandle.readlines())
    
        for category in ['short read pairs filtered out',
                         'empty read pairs filtered out',
                         'read pairs available',
                         'trimmed read pairs',
                         'untrimmed read pairs'
                        ]:
            pattern = r"(\d+)(?=\s*\(.*?(\d+\.\d+%).*?\)\s*"+f"{category})"
            match = re.search(pattern, lines)
            nread, perc = int(match.group(1)), float(match.group(2)[:-1])
            print(category, nread, perc)
    
            data[f'%{category}']=perc
            data[category]=nread
    return data

if __name__ == '__main__':
    parser = option_parser()
    (options, args) = parser.parse_args()
    
    log_files = options.input.split(' ')
    all_data = []
    for f in log_files:
        name = Path(f).name.split('-trimmed')[0]
        data = pd.Series(parse_skewer_log(f))
        data.name = name
        
        all_data.append(data)
    all_data = pd.concat(data, axis = 1)