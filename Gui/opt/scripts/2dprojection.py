from matplotlib.backends.backend_pdf import PdfPages

import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np

import sys

def parse_score(s):
    sc = s.split(" ")
    fc = sc[0]
    fp = sc[1]
    ft = sc[2]
    return (float(fc), float(fp), float(ft))

def dedup_scores (ss):
    scs = []
    for s in ss:
        if s not in scs:
            scs.append(s)
    return scs

def parse_results(fnm):
    scores = []
    tree = ET.parse(fnm)
    root = tree.getroot()
    for indiv in root.findall('Individual'):
        score = indiv.find('Score').text
        scores.append(parse_score(score))
        scs = dedup_scores(scores)
    return scs

def mk_fc_fp_plot(data, pp):
    fc = list(map(lambda t: t[0], data))
    fp = list(map(lambda t: t[1], data))
    f1 = plt.figure()
    plt.scatter(fc, fp, color="red")
    plt.xlabel('Material Cost')
    plt.ylabel('Precision')
    plt.xlim(min(fc), max(fc))
    plt.ylim(min(fp), max(fp))
    f1.savefig(pp)

def mk_fc_ft_plot(data, pp):
    fc = list(map(lambda t: t[0], data))
    ft = list(map(lambda t: t[1], data))
    f2 = plt.figure()
    plt.scatter(fc, ft, color="blue")
    plt.xlabel('Material Cost')
    plt.ylabel('Fabrication Time')
    plt.xlim(min(fc), max(fc))
    plt.ylim(min(ft), max(ft))
    f2.savefig(pp)

def mk_ft_fp_plot(data, pp):
    fp = list(map(lambda t: t[0], data))
    ft = list(map(lambda t: t[1], data))
    f3 = plt.figure()
    plt.scatter(fp, ft, color="green")
    plt.xlabel('Precision')
    plt.ylabel('Fabrication Time')
    plt.xlim(min(fp), max(fp))
    plt.ylim(min(ft), max(ft))
    f3.savefig(pp)

def main():
    fnm = sys.argv[1]
    nm = fnm.split(".")[0]
    scores = parse_results(fnm)
    fc_fp = list(map(lambda t:(t[0], t[1]), scores))
    fc_ft = list(map(lambda t:(t[0], t[2]), scores))
    fp_ft = list(map(lambda t:(t[1], t[2]), scores))
    mk_fc_fp_plot(fc_fp, nm + "_material_precision.pdf")
    mk_fc_ft_plot(fc_ft, nm + "_material_time.pdf")
    mk_ft_fp_plot(fp_ft, nm + "_time_precision.pdf")

main()
