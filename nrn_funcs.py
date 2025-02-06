import numpy as np

def spike_counter (t,V,th) :
    V = np.array(V)
    t = np.array(t)
    ind1 = []
    for i in range(np.shape(V)[0] -1) :

        if ((V[i] - th) <= 0.0  and (V[i+1] - th) > 0.0 ) :
            ind1.append(i)
    arr = t[ind1]

    return arr
