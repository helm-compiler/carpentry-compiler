#!/usr/bin/python

import xml.etree.ElementTree as ET
import sys
import os

def mk_xml (all_progs, op):
    r = ET.Element('root')
    for module in all_progs:
        (idx, progs) = module
        sid = str(idx)
        sub = ET.SubElement(r, 'sub', name = sid)
        for p in progs:
            pid = str(progs.index(p))
            ET.SubElement(sub, 'prog', name = pid).text = p
    tree = ET.ElementTree(r)
    tree.write(op)

def unify_sub_progs (d, op):
    all_progs = []
    for dnm, ss, fs in os.walk(d):
        for s in ss:
            module_progs = []
            idx = s
            for fnm in os.listdir(d + s):
                ip = d + s + "/" + fnm
                prog = open(ip).read()
                module_progs.append(prog)
            all_progs.append((idx, module_progs))
    mk_xml(all_progs, op)

#unify_sub_progs("../benchmarks/bookcase/subs/", "../benchmarks/bookcase/collapse/all_progs.xml")
unify_sub_progs("../benchmarks/bookcase/sub_sexps/", "../benchmarks/bookcase/collapse/all_progs.xml")
