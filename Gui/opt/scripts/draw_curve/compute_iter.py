import os
import sys

def load_hv(fn):
    f = open(fn,"r")
    hv = []
    for l in f:
        hv.append(float(l))

    return hv

if __name__ == "__main__":
    hv = load_hv('hyervolume.txt')

    for i in range(0, len(hv)):
        if i < 250:
            continue
        if abs(hv[i] - hv[i-250]) < 1:
            print(i)
            break
        