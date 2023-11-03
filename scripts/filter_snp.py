from scipy.stats import binom
import pandas as pd
import sys

if __name__ == '__main__':
    fname = sys.argv[1]
    f_handle = open(sys.argv[2],
                'w')
    chunk = 0
    for df in pd.read_csv(fname,
                sep = '\t', chunksize = 1000000, index_col = 0):
        df['total'] = df['n_ref']+df['n_alt']
        
        df['p_hetero']=df.apply(lambda row: binom.pmf(n=row['total'],
                                                        k=row['n_alt'], p = 0.5), axis = 1)
        
        df['p_homo']=df.apply(lambda row: binom.pmf(n=row['total'],
                                                        k=row['n_alt'], p = 0.99), axis = 1) # 1% seq error
        
        # remove those look like SNPs
        to_remove = df.loc[(df['p_hetero']>0.05)|(df['p_homo']>0.05)].index
        
        
        df.drop(to_remove, axis = 0, inplace = True)
        
        if chunk==0:
            df.to_csv(f_handle, sep = '\t', header = True, index = True)
        else:
            df.to_csv(f_handle, sep = '\t', header = False, index = True)
            
        chunk += 1