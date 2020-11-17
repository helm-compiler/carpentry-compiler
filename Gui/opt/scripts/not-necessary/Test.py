# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys, io
import xml.etree.ElementTree as ET


if __name__ == "__main__":
    taskName = "bookcase-pd-c5"
    collapsedFile = "../benchmarks/" + taskName + "/collapse/all_progs.xml"
    originalFile = "../benchmarks/" + taskName + "/collapse/ori_progs.xml"
    resultFile = "../benchmarks/" + taskName + "/collapse/result.xml"
    resultFileWithProgs = "../benchmarks/" + taskName + "/collapse/resultProgs.xml"

    oriProg = ET.parse(originalFile).getroot()
    sub = oriProg.findall('sub')
    subModules = {}
    for s in sub:
        subModules[int(s.attrib['name'])] = s.findall('prog')

    inxml = ET.parse(resultFile)
    root = inxml.getroot()
    for prog in root.findall('Individual'):
        attrib = prog.attrib
        order = prog.find('Order').text
        score = prog.find('Score').text
        tree = prog.find('Tree').text
        #print(order)
        #print(score)
        tree = tree.lstrip(',').rstrip(',').split(',')

        allProg = []
        for elmt in tree:
            (a, b) = elmt.split('-')
            (a, b) = (int(a), int(b))
            allProg.extend(subModules[a][b].text.rstrip().split('\n'))
        
        orderSVec = order.rstrip().split(' ')
        #print(len(allProg))
        #print(orderSVec)
        finalProg = '\n'
        maxInd = -1
        for o in orderSVec:
            #print(allProg[int(o)])
            finalProg += allProg[int(o)] + '\n'
            if int(o) > maxInd:
                maxInd = int(o)
        
        if maxInd+1 != len(orderSVec):
            print('ERROR')
            exit(-1)

        fp = ET.SubElement(prog, 'finalProg', len = str(len(orderSVec)))
        fp.text = finalProg
    
    inxml.write(resultFileWithProgs)