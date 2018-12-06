import numpy as np

# Read AutoLAB data file: two blank header lines and time; voltage data
def readAL(fname, I):
    data = np.loadtxt(fname,skiprows=2)
    m = np.insert(data, 2, values=I, axis=1)    # add I column
    m[:,0] = m[:,0] - m[0,0]                    # time begins at 0
    return(m)
