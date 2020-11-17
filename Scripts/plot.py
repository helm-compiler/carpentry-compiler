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

path = 'C:/Windows/Fonts/LeelawUI.ttf'
prop = font_manager.FontProperties(fname=path)
prop.set_weight = 'normal'
mpl.rcParams['font.family'] = prop.get_name()
mpl.rcParams['font.weight'] = 'normal'
fontSize = 22
labelSize = 14
dotSize3d = 25
dotSize2d = 45

def is_pareto(costs, return_mask = True):
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

def mk_fc_fp_plot(data, efc, efp, a):
    fc = list(map(lambda t: t[1], data))
    fp = list(map(lambda t: t[0], data))
    f_color = list(map(lambda t : "red", data))
    fc.append(float(efc))
    fp.append(float(efp))
    f_color.append("green")
    a.scatter(fc, fp, color=f_color, s= dotSize2d)
    a.set_ylabel('Material',fontproperties=prop, size = fontSize)
    a.set_xlabel('Precision',fontproperties=prop, size = fontSize)
    xLength = max(fc) - min(fc)
    yLength = max(fp) - min(fp)
    a.set_xlim(min(fc) - 0.1 * xLength, max(fc) + 0.1 * xLength)
    a.set_ylim(min(fp) - 0.1 * yLength, max(fp) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)
    a.tick_params(labelsize=labelSize)

def mk_fc_ft_plot(data, efc, eft, a):
    fc = list(map(lambda t: t[0], data))
    ft = list(map(lambda t: t[1], data))
    f_color = list(map(lambda t : "red", data))
    fc.append(float(efc))
    ft.append(float(eft))
    f_color.append("green")
    a.scatter(fc, ft, color=f_color, s= dotSize2d)
    a.set_ylabel('Fab. Time',fontproperties=prop, size = fontSize)
    a.set_xlabel('Material',fontproperties=prop, size = fontSize)
    xLength = max(fc) - min(fc)
    yLength = max(ft) - min(ft)
    a.set_xlim(min(fc) - 0.1 * xLength, max(fc) + 0.1 * xLength)
    a.set_ylim(min(ft) - 0.1 * yLength, max(ft) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)
    a.tick_params(labelsize=labelSize)


def mk_ft_fp_plot(data, eft, efp, a):
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[0], data))
    f_color = list(map(lambda t : "red", data))
    fp.append(float(efp))
    ft.append(float(eft))
    f_color.append("green")
    a.scatter(fp, ft, color=f_color, s= dotSize2d)
    a.set_ylabel('Precision',fontproperties=prop, size = fontSize)
    a.set_xlabel('Fab. Time',fontproperties=prop, size = fontSize)
    xLength = max(fp) - min(fp)
    yLength = max(ft) - min(ft)
    a.set_xlim(min(fp) - 0.1 * xLength, max(fp) + 0.1 * xLength)
    a.set_ylim(min(ft) - 0.1 * yLength, max(ft) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)


def mk_fc_fp_ft_plot(data, ex_score, ax):
    #ax.view_init(25, 0)
    fc = list(map(lambda t: t[0], data))
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[2], data))
    f_color = list(map(lambda t : "red", data))
    efc = ex_score[0]
    efp = ex_score[1]
    eft = ex_score[2]
    
    #fc.append(float(efc))
    #fp.append(float(efp))
    #ft.append(float(eft))
    ax.scatter(fc, fp, ft, c='r',alpha=1, s= dotSize3d)
    ax.scatter(efc, efp, eft, c='g', alpha=1.0, s= dotSize3d)
    ax.set_xlabel("Material",fontproperties=prop, size = fontSize,labelpad=16)
    ax.set_ylabel("Precision",fontproperties=prop, size = fontSize, labelpad=16)
    ax.set_zlabel("Fab. time",fontproperties=prop, size = fontSize, rotation = 180, labelpad=16)
    ax.tick_params(labelsize=labelSize, pad=6)


def build_xml(fnm, nnm, uniqIndex):
    scores = []
    tree = ET.parse(fnm)
    root = tree.getroot()
    cnt = 0
    minInd = [0, 0, 0]
    minVal = [1e30, 1e30, 1e30]

    for indiv in root.findall('Individual'):
        if cnt not in uniqIndex:
            root.remove(indiv)
            cnt += 1
            continue
        else:
            score = indiv.find('Score').text
            s = parse_score(score)
            for i in range(0, 3):
                if s[i] < minVal[i]:
                    minVal[i] = s[i]
                    minInd[i] = indiv.attrib['ID']
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

def main():
    print(sys.argv[1])
    fnm = sys.argv[1]
    nm = fnm.split(".xml")[0]

    ex_fc = sys.argv[2]
    ex_fp = sys.argv[3]
    ex_ft = sys.argv[4]
    ex_sc = (float(ex_fc), float(ex_fp), float(ex_ft))

    scores, indexUnique = parse_results(fnm)
    
    res = is_pareto(np.asarray(scores))

    npScore = np.asarray(scores)
    fc_fp_ft = npScore[res]
    scores = npScore[res]
    ''' For birdhouse
    scores[:, 2] += 2
    fc_fp_ft[:, 2] += 2
    '''
    usefulIndividual = np.asarray(indexUnique)[res]
    print(usefulIndividual)
    nnm = nm + "-processed.xml"
    build_xml(fnm, nnm, usefulIndividual)

    fc_fp = list(map(lambda t:(t[0], t[1]), scores))
    fc_ft = list(map(lambda t:(t[0], t[2]), scores))
    fp_ft = list(map(lambda t:(t[1], t[2]), scores))

    fig = plt.figure()

    mk_fc_fp_ft_plot(fc_fp_ft, ex_sc, plt.subplot2grid((30, 40), (0, 0), rowspan=30, colspan=22, projection='3d'))
    mk_fc_fp_plot(fc_fp, ex_sc[1], ex_sc[0], plt.subplot2grid((30, 40), (0, 25), rowspan=7, colspan=12))
    mk_fc_ft_plot(fc_ft, ex_sc[0], ex_sc[2], plt.subplot2grid((30, 40), (10, 25), rowspan=7, colspan=12))
    mk_ft_fp_plot(fp_ft, ex_sc[1], ex_sc[2], plt.subplot2grid((30, 40), (20, 25),  rowspan=7,colspan=12))

    plt.tight_layout()
    fig.set_size_inches(18, 10, forward=True)
    #fig.savefig(nm + "_plots.png", dpi=fig.dpi)
    plt.show()
    fig.savefig(nm + "_plots.svg", dpi=fig.dpi, bbox_inches="tight")

main()

#en_vs_precision("enode-precision.txt")
#en_vs_time("enode-time.txt")
