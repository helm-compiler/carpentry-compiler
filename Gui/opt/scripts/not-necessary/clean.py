import os, sys
import xml.etree.ElementTree as ET

submodules = {}

def clean_egraph (egraph, path):
    tree = ET.parse(egraph)
    root = tree.getroot()
    for prog in root.findall('EClass'):
        #print(prog.attrib['ID'])
        submodules[prog.attrib['ID']] = 1

    for fnm in os.listdir(path):
        if fnm.endswith('xml'):
            fileparts = fnm.split('.')
            if fileparts[0] not in submodules:
                print(fileparts[0])
                os.remove(path+fnm)
                #TODO: if it is not existed, then delete this file!

print('hello')
#clean_egraph("../benchmarks/bookcase-random/egraph.xml", "../benchmarks/bookcase-random/programs/")
clean_egraph("../benchmarks/bench/egraph.xml", "../benchmarks/bench/programs/")