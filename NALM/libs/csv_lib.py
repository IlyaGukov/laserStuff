import pandas as pd
import numpy as np

def write_csv(data,cols,name):
    raw_data = {l:uu for l,uu in zip(cols,data)}
    df = pd.DataFrame(raw_data, columns = cols)
    df.to_csv(name, index=False)

def read_csv(name, cols):
    temp = pd.read_csv(name)
    res = np.array([np.array([complex(aa) for aa in temp[str(l)]]) for l in cols])
    return res