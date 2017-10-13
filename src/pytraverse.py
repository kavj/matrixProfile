import numpy as np


def mpw(ts,m):
    n = ts.shape-m+1
    for i in range (n):
        print(i)
    return

def mpBlock(ts,m):
    n = (ts.shape[0])
    blk = 40
    step = blk - m + 1
    last = n - m + 1
    print(last)
    if last > blk:
        lastAligned = last - last%blk
    else:
        lastAligned = 0
    for i in range (0,lastAligned,blk):
        print('blk')
        for j in range(i,i+blk):
            print(j)
    for j in range(lastAligned,last):
            print(j)
    return


