import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from mpl_toolkits.mplot3d import Axes3D

import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import sys


path = 'C:/Users/wcm94/AppData/Local/Microsoft/Windows/Fonts/LinLibertine_DRah.ttf'
prop = font_manager.FontProperties(fname=path)
prop.set_weight = 'normal'
mpl.rcParams['font.family'] = prop.get_name()
mpl.rcParams['font.weight'] = 'normal'
fontSize = 22
labelSize = 14
dotSize3d = 25
dotSize2d = 15

def is_pareto_efficient(costs, return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient

import matplotlib.transforms
import sys

font = {'family' : 'sans',
        'weight' : 'normal',
        'size' : 10}

matplotlib.rc('font', **font)

def parse_score(s):
    sc = s.split(" ")
    fc = sc[0]
    fp = sc[1]
    ft = sc[2]
    return (float(fc), float(fp), float(ft))

def dedup_scores (ss):
    scs = []
    indexUnique = []
    for s in ss:
        if s not in scs:
            scs.append(s)
            indexUnique.append(ss.index(s))
    return scs, indexUnique

def parse_results(fnm):
    scores = []
    tree = ET.parse(fnm)
    root = tree.getroot()
    for indiv in root.findall('Individual'):
        score = indiv.find('Score').text
        scores.append(parse_score(score))

    scs,indexUnique = dedup_scores(scores)
    return scs,indexUnique

def mk_fc_fp_plot(data, edata,bdata, a):
    fc = list(map(lambda t: t[1], data))
    fp = list(map(lambda t: t[0], data))
    f_color = list(map(lambda t : "red", data))
    efc = list(map(lambda t: t[1], edata))
    efp = list(map(lambda t: t[0], edata))
    ef_color = list(map(lambda t : "green", edata))
    bfc = list(map(lambda t: t[1], bdata))
    bfp = list(map(lambda t: t[0], bdata))
    bf_color = list(map(lambda t : "blue", bdata))
    a.scatter(efc, efp, color=ef_color, s= dotSize2d)
    a.scatter(fc, fp, color=f_color, s= dotSize2d)
    a.scatter(bfc, bfp, color=bf_color, s= dotSize2d)
    a.set_ylabel('Material',fontproperties=prop, size = fontSize)
    a.set_xlabel('Precision',fontproperties=prop, size = fontSize)
    maxX = max([max(fc), max(efc), max(bfc)])
    minX = min([min(fc), min(efc), min(bfc)])
    maxY = max([max(fp), max(efp), max(bfp)])
    minY = min([min(fp), min(efp), min(bfp)])
    xLength =  maxX - minX
    yLength = maxY - minY
    a.set_xlim(minX - 0.1 * xLength, maxX + 0.1 * xLength)
    a.set_ylim(minY - 0.1 * yLength, maxY + 0.1 * yLength)

def mk_fc_ft_plot(data, edata, bdata, a):
    fc = list(map(lambda t: t[0], data))
    ft = list(map(lambda t: t[1], data))
    f_color = list(map(lambda t : "red", data))
    efc = list(map(lambda t: t[0], edata))
    eft = list(map(lambda t: t[1], edata))
    ef_color = list(map(lambda t : "green", edata))
    bfc = list(map(lambda t: t[0], bdata))
    bft = list(map(lambda t: t[1], bdata))
    bf_color = list(map(lambda t : "blue", bdata))
    a.scatter(fc, ft, color=f_color,   s= dotSize2d)
    a.scatter(efc, eft, color=ef_color,   s= dotSize2d)
    a.scatter(bfc, bft, color=bf_color,   s= dotSize2d)
    a.set_ylabel('Fab. Time',fontproperties=prop, size = fontSize)
    a.set_xlabel('Material',fontproperties=prop, size = fontSize)
    maxX = max([max(fc), max(efc), max(bfc)])
    minX = min([min(fc), min(efc), min(bfc)])
    maxY = max([max(ft), max(eft), max(bft)])
    minY = min([min(ft), min(eft), min(bft)])
    print(bft)
    xLength =  maxX - minX
    yLength = maxY - minY
    a.set_xlim(minX - 0.1 * xLength, maxX + 0.1 * xLength)
    a.set_ylim(minY - 0.1 * yLength, maxY + 0.1 * yLength)


def mk_ft_fp_plot(data, edata, bdata, a):
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[0], data))
    f_color = list(map(lambda t : "red", data))
    efp = list(map(lambda t: t[1], edata))
    eft = list(map(lambda t: t[0], edata))
    ef_color = list(map(lambda t : "green", edata))
    bfp = list(map(lambda t: t[1], bdata))
    bft = list(map(lambda t: t[0], bdata))
    bf_color = list(map(lambda t : "blue", bdata))
    a.scatter(fp, ft, color=f_color, s= dotSize2d)
    a.scatter(efp, eft, color=ef_color, s= dotSize2d)
    a.scatter(bfp, bft, color=bf_color, s= dotSize2d)
    a.set_ylabel('Precision',fontproperties=prop, size = fontSize)
    a.set_xlabel('Fab. Time',fontproperties=prop, size = fontSize)
    maxX = max([max(fp), max(efp), max(bfp)])
    minX = min([min(fp), min(efp), min(bfp)])
    maxY = max([max(ft), max(eft), max(bft)])
    minY = min([min(ft), min(eft), min(bft)])
    xLength =  maxX - minX
    yLength = maxY - minY
    a.set_xlim(minX - 0.1 * xLength, maxX + 0.1 * xLength)
    a.set_ylim(minY - 0.1 * yLength, maxY + 0.1 * yLength)

def mk_fc_fp_ft_plot(data, edata, bdata, ax):
    fc = list(map(lambda t: t[0], data))
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[2], data))
    f_color = list(map(lambda t : "red", data))
    efc = list(map(lambda t: t[0], edata))
    efp = list(map(lambda t: t[1], edata))
    eft = list(map(lambda t: t[2], edata))
    f_color = list(map(lambda t : "green", edata))
    bfc = list(map(lambda t: t[0], bdata))
    bfp = list(map(lambda t: t[1], bdata))
    bft = list(map(lambda t: t[2], bdata))
    f_color = list(map(lambda t : "blue", bdata))
    #fc.append(float(efc))
    #fp.append(float(efp))
    #ft.append(float(eft))
    ax.scatter(fc, fp, ft, c='r',alpha=1.0, s= dotSize3d)
    ax.scatter(efc, efp, eft, c='g',alpha=1.0, s= dotSize3d)
    ax.scatter(bfc, bfp, bft, c='b',alpha=1.0, s= dotSize3d)
    ax.set_xlabel("Material",fontproperties=prop, size = fontSize,labelpad=16)
    ax.set_ylabel("Precision",fontproperties=prop, size = fontSize, labelpad=16)
    ax.set_zlabel("Fab. time",fontproperties=prop, size = fontSize, rotation = 180, labelpad=16)

def build_xml(fnm, nnm, uniqIndex):
    scores = []
    tree = ET.parse(fnm)
    root = tree.getroot()
    cnt = 0
    minInd = [0, 0, 0]
    minVal = [1e30, 1e30, 1e30]

    for indiv in root.findall('Individual'):
        score = indiv.find('Score').text
        s = parse_score(score)
        for i in range(0, 3):
            if s[i] < minVal[i]:
                minVal[i] = s[i]
                minInd[i] = indiv.attrib['ID']

        if cnt not in uniqIndex:
            root.remove(indiv)
        cnt += 1
    tree.write(nnm)

    print(minInd)
    print(minVal)

def parse_en_time(l):
    ss = l.split(":")
    enodes = ss[0].split(",")
    enodes[0] = enodes[0].split("[")[1]
    enodes[len(enodes)-1] = enodes[len(enodes)-1].split("]")[0]
    ens = list(map(lambda e : float(e), enodes))
    times = ss[1].split(",")
    times[0] = times[0].split("[")[1]
    times[len(times)-1] = times[len(times)-1].split("]")[0]
    ts = list(map(lambda t : float(t), times))
    name = ss[2]
    return (ens, ts, name)

def parse_en_prec(l):
    ss = l.split(":")
    enodes = ss[0].split(",")
    enodes[0] = enodes[0].split("[")[1]
    enodes[len(enodes)-1] = enodes[len(enodes)-1].split("]")[0]
    ens = list(map(lambda e : float(e), enodes))
    precs = ss[1].split(",")
    precs[0] = precs[0].split("[")[1]
    precs[len(precs)-1] = precs[len(precs)-1].split("]")[0]
    ps = list(map(lambda p : float(p), precs))
    name = ss[2]
    return (ens, ps, name)

def en_vs_precision(fnm):
    fin = open(fnm, "r")
    fig = plt.figure()
    for l in fin:
        ens, ps, name = parse_en_prec(l)
        ens.sort()
        ps.sort()
        n = len(ens)
        plt.plot(ens, ps, label=name)
        #plt.annotate(name, xy=(ens[n - 1], ps[n - 1]), xytext=(ens[n - 1] + 0.5, ps[n - 1] + 0.5))
    plt.legend()
    plt.xlabel('e-nodes')
    plt.ylabel('precision')
    fig.savefig("ens-precision.png", dpi=fig.dpi, bbox_inches="tight")

def en_vs_time(fnm):
    fin = open(fnm, "r")
    fig = plt.figure()
    for l in fin:
        ens, ts, name = parse_en_prec(l)
        ens.sort()
        ts.sort()
        n = len(ens)
        plt.plot(ens, ts, label=name)
        #plt.annotate(name, xy=(ens[n - 1], ts[n - 1]), xytext=(ens[n - 1] + 0.5, ts[n - 1] + 0.5))
    plt.legend()
    plt.xlabel('e-nodes')
    plt.ylabel('time')
    fig.savefig("ens-time.png", dpi=fig.dpi, bbox_inches="tight")

def load_comb(fnm):
    f = open(fnm, "r")
    scores = []
    for line in f:
        parts = line.split(' ')
        score = (float(parts[0]), float(parts[1]), float(parts[2]))
        scores.append(score)
    return scores


def main():
    fnm = sys.argv[1]
    nm = fnm.split(".xml")[0]

    scores, indexUnique = parse_results(fnm)
    
    res = is_pareto_efficient(np.asarray(scores))
    npScore = np.asarray(scores)
    scores = npScore[res]
    fc_fp_ft = npScore[res]
    usefulIndividual = np.asarray(indexUnique)[res]

    nnm = nm + "-processed.xml"
    build_xml(fnm, nnm, usefulIndividual)

    cscores = load_comb(sys.argv[2])
    cres = is_pareto_efficient(np.asarray(cscores))
    npCscore = np.asarray(cscores)
    cscores = npCscore[cres]
    cfc_fp_ft = npCscore[cres]
    print(np.where(cres==True))
    print(np.size(cres))
    print("Value index", np.argmin(cscores, axis=0))

    bscores = load_comb(sys.argv[3])
    bres = is_pareto_efficient(np.asarray(bscores))
    npBscore = np.asarray(bscores)
    bscores = npBscore[bres]
    bfc_fp_ft = npBscore[bres]

    bfc_fp = list(map(lambda t:(t[0], t[1]), bscores))
    bfc_ft = list(map(lambda t:(t[0], t[2]), bscores))
    bfp_ft = list(map(lambda t:(t[1], t[2]), bscores))

    cfc_fp = list(map(lambda t:(t[0], t[1]), cscores))
    cfc_ft = list(map(lambda t:(t[0], t[2]), cscores))
    cfp_ft = list(map(lambda t:(t[1], t[2]), cscores))

    fc_fp = list(map(lambda t:(t[0], t[1]), scores))
    fc_ft = list(map(lambda t:(t[0], t[2]), scores))
    fp_ft = list(map(lambda t:(t[1], t[2]), scores))

    

    fig = plt.figure()

    mk_fc_fp_ft_plot(fc_fp_ft, cfc_fp_ft, bfc_fp_ft, plt.subplot2grid((30, 40), (0, 0), rowspan=30, colspan=22, projection='3d'))
    mk_fc_fp_plot(fc_fp, cfc_fp, bfc_fp, plt.subplot2grid((30, 40), (0, 25), rowspan=7, colspan=12))
    mk_fc_ft_plot(fc_ft, cfc_ft, bfc_ft, plt.subplot2grid((30, 40), (10, 25), rowspan=7, colspan=12))
    mk_ft_fp_plot(fp_ft, cfp_ft, bfp_ft, plt.subplot2grid((30, 40), (20, 25),  rowspan=7,colspan=12))

    plt.tight_layout()
    fig.set_size_inches(18, 10, forward=True)
    #fig.savefig(nm + "_plots.png", dpi=fig.dpi)
    plt.show()
    fig.savefig(nm + "_plots.svg", dpi=fig.dpi, bbox_inches="tight")

main()

#en_vs_precision("enode-precision.txt")
#en_vs_time("enode-time.txt")
