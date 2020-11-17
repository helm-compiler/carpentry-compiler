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

from matplotlib import cm
from collections import OrderedDict


path = 'C:/Windows/Fonts/LeelawUI.ttf'
prop = font_manager.FontProperties(fname=path)
prop.set_weight = 'normal'
mpl.rcParams['font.family'] = prop.get_name()
mpl.rcParams['font.weight'] = 'normal'
fontSize = 22
labelSize = 14
dotSize3d = 25
dotSize2d = 45
delta_label =0.000;

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

def mk_fc_fp_plot(data,f_color,axis_range, a, with_label = True):
    fc = list(map(lambda t: t[1], data))
    fp = list(map(lambda t: t[0], data))
    a.scatter(fc, fp, color=f_color, s= dotSize2d)
    a.set_ylabel('Material',fontproperties=prop, size = fontSize)
    a.set_xlabel('Precision',fontproperties=prop, size = fontSize)
    xLength = max(fc) - min(fc)
    yLength = max(fp) - min(fp)
    a.set_xlim(min(fc) - 0.1 * xLength, max(fc) + 0.1 * xLength)
    a.set_ylim(min(fp) - 0.1 * yLength, max(fp) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)
    a.tick_params(labelsize=labelSize)
    if with_label :
        for index in range(0,len(data)):
            a.annotate(str(index),(fc[index]+delta_label,fp[index]+delta_label))
    if len(axis_range)!=0:
        a.axis([axis_range[2],axis_range[3],axis_range[0],axis_range[1]])

def mk_fc_ft_plot(data,f_color,axis_range,  a, with_label = True):
    fc = list(map(lambda t: t[0], data))
    ft = list(map(lambda t: t[1], data))
    a.scatter(fc, ft, color=f_color, s= dotSize2d)
    a.set_ylabel('Fab. Time',fontproperties=prop, size = fontSize)
    a.set_xlabel('Material',fontproperties=prop, size = fontSize)
    xLength = max(fc) - min(fc)
    yLength = max(ft) - min(ft)
    a.set_xlim(min(fc) - 0.1 * xLength, max(fc) + 0.1 * xLength)
    a.set_ylim(min(ft) - 0.1 * yLength, max(ft) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)
    a.tick_params(labelsize=labelSize)
    if with_label :
        for index in range(0,len(data)):
            a.annotate(str(index),(fc[index]+delta_label,ft[index]+delta_label))
    if len(axis_range)!=0:
        a.axis([axis_range[0],axis_range[1],axis_range[4],axis_range[5]])


def mk_ft_fp_plot(data,f_color,axis_range, a, with_label = True):
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[0], data))
    a.scatter(fp, ft, color=f_color, s= dotSize2d)
    a.set_ylabel('Precision',fontproperties=prop, size = fontSize)
    a.set_xlabel('Fab. Time',fontproperties=prop, size = fontSize)
    xLength = max(fp) - min(fp)
    yLength = max(ft) - min(ft)
    a.set_xlim(min(fp) - 0.1 * xLength, max(fp) + 0.1 * xLength)
    a.set_ylim(min(ft) - 0.1 * yLength, max(ft) + 0.1 * yLength)
    a.tick_params(labelsize=labelSize)
    if with_label :
        for index in range(0,len(data)):
            a.annotate(str(index),(fp[index]+delta_label,ft[index]+delta_label))
    if len(axis_range)!=0:
        a.axis([axis_range[4],axis_range[5],axis_range[2],axis_range[3]])

def color_mapping(isolevel):
    r=0.0
    g=0.0
    b=0.0
    if isolevel >= 0.0 and isolevel <= 0.25:
        r = 0.0
        g = isolevel / 0.25
        b = 1
    if isolevel > 0.25 and isolevel <= 0.50:
        r = 0
        g = 1
        b = 1 - (isolevel - 0.25) / 0.25   
    if (isolevel > 0.50 and isolevel <= 0.75):
        r = (isolevel - 0.50) / 0.25
        g = 1
        b = 0
    if (isolevel > 0.75 and isolevel <= 1.0):
        r = 1
        g = 1 - (isolevel - 0.75) / 0.25
        b = 0
    if (isolevel < 0.0):
        r = 0.0
        g = 0.0
        b = 0.0
    if (isolevel > 1.0):
        r = 0.5
        g = 0.0
        b = 0.0
    return [r,g,b]

def mk_fc_fp_ft_plot(data,f_color,axis_range, ax, with_label = True):
    #ax.view_init(25, 0)
    fc = list(map(lambda t: t[0], data))
    fp = list(map(lambda t: t[1], data))
    ft = list(map(lambda t: t[2], data))
    #f_color = list(map(lambda t : [0.0,1.0,0.0], data))

    ax.scatter(fc, fp, ft, c=f_color,alpha=1, s= dotSize3d)
    ax.set_xlabel("Material",fontproperties=prop, size = fontSize,labelpad=16)
    ax.set_ylabel("Precision",fontproperties=prop, size = fontSize, labelpad=16)
    ax.set_zlabel("Fab. time",fontproperties=prop, size = fontSize, rotation = 180, labelpad=16)
    ax.tick_params(labelsize=labelSize, pad=6)

    if with_label :
        for index in range(0,len(data)):
            ax.text(fc[index],fp[index],ft[index],str(index))

    if len(axis_range)!=0:
        ax.axis([axis_range[0],axis_range[1],axis_range[2],axis_range[3]])
        ax.set_zlim(axis_range[4],axis_range[5])

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

def draw_fig(fc_fp_ft_all,fc_fp_all,fc_ft_all,fp_ft_all,f_color_all,path,axis_range, with_label = True, show = False):
    fig = plt.figure()

    mk_fc_fp_ft_plot(np.asarray(fc_fp_ft_all),f_color_all,axis_range, plt.subplot2grid((30, 40), (0, 0), rowspan=30, colspan=22, projection='3d'),with_label)
    mk_fc_fp_plot(fc_fp_all,f_color_all, axis_range, plt.subplot2grid((30, 40), (0, 25), rowspan=7, colspan=12),with_label)
    mk_ft_fp_plot(fp_ft_all,f_color_all, axis_range, plt.subplot2grid((30, 40), (10, 25),  rowspan=7,colspan=12),with_label)
    mk_fc_ft_plot(fc_ft_all,f_color_all, axis_range, plt.subplot2grid((30, 40), (20, 25), rowspan=7, colspan=12),with_label)

    #plt.tight_layout()
    fig.set_size_inches(18, 10, forward=True)
    fig.savefig(path, dpi=fig.dpi)
    if show:
        plt.show()
    #fig.savefig(path, dpi=fig.dpi, bbox_inches="tight")

def compute_axis_range(mm_scores_all):
    #compute axis range
    #material precision Fab. Time
    maxInColumns = np.amax(mm_scores_all, axis=0)
    minInColumns = np.amin(mm_scores_all, axis=0)
    axis_range=[minInColumns[0],maxInColumns[0],minInColumns[1],maxInColumns[1],minInColumns[2],maxInColumns[2]]
    for index in range(0,3):
        if maxInColumns[index]==minInColumns[index]:
            axis_range[2*index]=axis_range[2*index]-0.5
            axis_range[2*index+1]=axis_range[2*index+1]+0.5
        else:
            length = (maxInColumns[index] - minInColumns[index])*0.05
            axis_range[2*index]=axis_range[2*index]-length
            axis_range[2*index+1]=axis_range[2*index+1]+length
    return axis_range
def main():



    f_color_all =[]
    fc_fp_ft_all = []
    fc_fp_all =[]
    fc_ft_all = []
    fp_ft_all = []
    scores_all=[]
    design_all = []
    design_local_all = []
    mm_scores_all=[]

    design_nb = int(sys.argv[2])

    f = open(sys.argv[1]+'\\log\\unique_fronts.txt')
    iters_all = []
    for i in range(0,int(f.readline().split()[0])):
        iters_all.append(list(map(int,f.readline().split())))
    print("unique_fronts: ")
    for i in range(0,len(iters_all)):
        print(iters_all[i])
    f.close()

    
    for i in range(0,int(sys.argv[2])):
        print(color_mapping(float(i)/design_nb))

    for i in range(0,int(sys.argv[2])):
        fnm = sys.argv[1]+'\\design_'+str(i)+'\\egraph\\collapse\\result.xml';
        print(fnm)
        scores, indexUnique = parse_results(fnm)
        res = is_pareto(np.asarray(scores))
        for j in range(0,len(res)):
            if res[j]:
                print("Pareto Index: "+str(j))
        npScore = np.asarray(scores)
        scores = npScore[res]
        for score in scores:
            print("Front: "+str(score[0])+" "+str(score[1])+" "+str(score[2]))

        scores_all.append(scores)
        
        j=0
        for score in scores:
            mm_scores_all.append(score)
            f_color_all.append(color_mapping(float(i)/design_nb))
            fc_fp_ft_all.append([score[0],score[1],score[2]])
            fc_fp_all.append([score[0],score[1]])
            fc_ft_all.append([score[0],score[2]])
            fp_ft_all.append([score[1],score[2]])
            design_all.append(i)
            design_local_all.append(j)
            j=j+1

    #compute axis range
    axis_range = compute_axis_range(mm_scores_all)

    path=sys.argv[1] + "\\result_plots_all.png"
    draw_fig(fc_fp_ft_all,fc_fp_all,fc_ft_all,fp_ft_all,f_color_all, path, axis_range,False, False)

    #-----

    #-----

    for iter in range(0,len(iters_all)):
        res_=[]
        for i in range(0,len(fc_fp_ft_all)):
            res_.append(design_all[i] in iters_all[iter])
        res_=np.asarray(res_)
        fc_fp_ft_all_ = np.asarray(fc_fp_ft_all)[res_].tolist()
        fc_fp_all_ = np.asarray(fc_fp_all)[res_].tolist()
        fc_ft_all_ = np.asarray(fc_ft_all)[res_].tolist()
        fp_ft_all_ = np.asarray(fp_ft_all)[res_].tolist()
        f_color_all_ = np.asarray(f_color_all)[res_].tolist()
        design_all_ = np.asarray(design_all)[res_].tolist()
        path=sys.argv[1] + "\\result_plots_pareto_"+str(iter)+".png"
        draw_fig(fc_fp_ft_all_,fc_fp_all_,fc_ft_all_,fp_ft_all_,f_color_all_, path, axis_range,False, False)


    res = is_pareto(np.asarray(fc_fp_ft_all))
    fc_fp_ft_all_ = np.asarray(fc_fp_ft_all)[res].tolist()
    fc_fp_all_ = np.asarray(fc_fp_all)[res].tolist()
    fc_ft_all_ = np.asarray(fc_ft_all)[res].tolist()
    fp_ft_all_ = np.asarray(fp_ft_all)[res].tolist()
    f_color_all_ = np.asarray(f_color_all)[res].tolist()
    design_all_ = np.asarray(design_all)[res].tolist()
    design_local_all_ = np.asarray(design_local_all)[res].tolist()
    #design_local_all

    with open(sys.argv[1]+"\\pareto_fronts.txt","w") as f:
        for index in range(0,len(design_all_)):
            str_score = " : "+str(fc_fp_ft_all_[index][0])+" "+str(fc_fp_ft_all_[index][1])+" "+str(fc_fp_ft_all_[index][2]);
            f.write(str(index)+ " : "+str(design_all_[index])+"   "+str(design_local_all_[index])+str_score+"\n")
    
    path=sys.argv[1] + "\\result_plots_pareto.png"
    draw_fig(fc_fp_ft_all_,fc_fp_all_,fc_ft_all_,fp_ft_all_,f_color_all_, path, axis_range,False, False)

    for index in range(0,len(scores_all)):
        print("Plot: "+str(index))
        if index not in design_all_:
            print("dominate")
        f_color =[]
        fc_fp_ft = []
        fc_fp =[]
        fc_ft = []
        fp_ft = []
        for score in scores_all[index]:
            f_color.append(color_mapping(float(index)/design_nb)),
            fc_fp_ft.append([score[0],score[1],score[2]])
            fc_fp.append([score[0],score[1]])
            fc_ft.append([score[0],score[2]])
            fp_ft.append([score[1],score[2]])

        if index not in design_all_:
             with open(sys.argv[1]+"\\a_"+str(index)+".txt","w") as f:
                 for score_index in range(0,len(scores_all[index])):
                     str_score = str(scores_all[index][score_index][0])+" "+str(scores_all[index][score_index][1])+" "+str(scores_all[index][score_index][2]);
                     f.write(str(score_index)+":"+str_score +" : ")
                     for index_ in range(0,len(design_all_)):
                         if fc_fp_ft_all_[index_][0] <= scores_all[index][score_index][0] and fc_fp_ft_all_[index_][1] <= scores_all[index][score_index][1] and fc_fp_ft_all_[index_][2] <= scores_all[index][score_index][2]:
                             f.write(str(index_)+", ")
                     f.write("\n")

        if index not in design_all_:
            path=sys.argv[1] + "\\result_plots_"+str(index)+"_dominated_fronts_"+str(len(scores_all[index]))+".png"
        else:
            path=sys.argv[1] + "\\result_plots_"+str(index)+"_fronts_"+str(len(scores_all[index]))+".png"
        draw_fig(fc_fp_ft,fc_fp,fc_ft,fp_ft,f_color,path,axis_range,True)



main()

#en_vs_precision("enode-precision.txt")
#en_vs_time("enode-time.txt")

