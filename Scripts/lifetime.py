import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def parse(file):

    print(file)
    if not os.path.exists(file) :
        print('not exists...')
        return

    a = np.loadtxt(file)
    a = np.transpose(a)
    plt.grid(color="k", linestyle=":")

    x = a[0] 
    y = a[1]  
    plt.scatter(x, y, alpha=1.0)   
    plt.xlabel('Iteration Index', fontsize=20)
    plt.ylabel('E-Node Index', fontsize=20)
    plt.savefig(file+".png")
    plt.show()
    return

def main():
    print(sys.argv[1])
    parse(sys.argv[1])

main()