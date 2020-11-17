# -*- coding: utf-8 -*-
#!/usr/bin/python

import xml.etree.ElementTree as ET
import sys
import os
import re
import itertools

uuid_rgxp = "[0-9a-fA-F]{8}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{12}"

enodes = 0

def mk_sexp (assigns):
    sequify = []
    if len(assigns) >= 3:
        #print "3 assigns"
        # for a in assigns:
        #    print(a)
        #    print("\n")
        for x in range(len(assigns) - 2):
            sequify.append ("(Seq " + assigns[x])
        final_seq = "(Seq " + assigns[len(assigns) - 2] + assigns[len(assigns) - 1]
        num_close_paren = len(assigns) - 1
        close_parens = ')' * num_close_paren
        return (''.join(sequify) + final_seq + close_parens)
    elif len(assigns) == 2:
        # print("2 assigns")
        final_seq = "(Seq " + assigns[0] + assigns[1] + ")"
        return final_seq
    elif len(assigns) == 1:
        # print("1 assign")
        return assigns[0]
    else:
        return "(Empty)"

def chopsaw (lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    # TODO: this sometimes is Lumber, sometimes is Var
    lumber = "( Lumber " + ids[0] + ")"
    face = "( Face " + ids[1] + ")"
    edge = "( Edge " + ids[2] + ")"
    rsplt = rhs.split(',')
    ang1 = "(Angle (Float " + rsplt[3].split('(')[1] + " ))"
    ang2 = "(Angle (Float " + rsplt[4] + " ))"
    lnth = "(Length (Float " + rsplt[5].split(')')[0] + " ))"
    ref = "(Refc " + ang1 + ang2 + lnth + ")"
    # TODO: everything is stackable for now
    stk = "(Bool true)"
    height = "(Height (Float " + rsplt[6].split(')')[0] + " ))"
    rhs = "(Chopsaw " + lumber + face + edge + ref + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign

def mk_xy(p):
    px = "(Float " + p.split('/')[0] + " )"
    py = "(Float " + p.split('/')[1] + " )"
    pt = "(Tup " + px + py + " )"
    return pt

def mk_path (ps):
    path = []
    p1 = ps[0].split('(')[1]
    path.append (mk_xy(p1))
    # first and last are taken care of above
    for i in range(1, len(ps) - 2):
        path.append(mk_xy(ps[i]))
    pn = ps[len(ps) - 1].split(')')[0]
    path.append (mk_xy(pn))
    spath = "(Path "
    for i in range (0, len(path)):
        spath = spath + path[i]
    return (spath + ")")

def bandsaw (lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    xprm = "(XPrime (Float " + rsplt[len(rsplt) - 2].split('(')[1] + " ))"
    lnth = "(Length (Float " + rsplt[len(rsplt) - 1].split(')')[0] + " ))"
    ps = []
    # last two list elements will always be about ref.
    for i in range(1, len(rsplt) - 2):
        ps.append(rsplt[i])
    path = mk_path(ps)
    ref = "(Refb " + xprm + lnth + ")"
    # TODO: everything is stackable for now
    stk = "(Bool true)"
    height = "(Height (Float 0.75))"
    rhs = "(Bandsaw " + lumber + path + ref + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign

def jigsaw (lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    ps = []
    # last two list elements will always be about ref.
    for i in range(1, len(rsplt) - 2):
        ps.append(rsplt[i])
    path = mk_path(ps)
    # TODO: everything is stackable for now
    stk = "(Bool true)"
    height = "(Height (Float 0.75))"
    refj = "(Refj)"
    rhs = "(Jigsaw " + lumber + path + stk + refj + height + ")"
    assign = lhs + rhs + " )"
    return assign

def drill (lhs, rhs):
    return "TODO_DRILL"

def tracksaw (lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    # TODO: this sometimes is Lumber, sometimes is Var
    lumber = "( Lumber " + ids[0] + ")"
    face = "( Face " + ids[1] + ")"
    edge = "( Edge " + ids[2] + ")"
    rsplt = rhs.split(',')
    ang1 = "(Angle (Float " + rsplt[3].split('(')[1] + " ))"
    ang2 = "(Angle (Float " + rsplt[4] + " ))"
    lnth = "(Length (Float " + rsplt[5].split(')')[0] + " ))"
    ref = "(Refk " + ang1 + ang2 + lnth + ")"
    # TODO: everything is stackable for now
    stk = "(Bool true)"
    height = "(Height (Float " + rsplt[6].split(')')[0] + " ))"
    rhs = "(Tracksaw " + lumber + face + edge + ref + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign

# TODO: right now needs to be ran in the scripts directory.
def parse_xml(i):
    global enodes
    ip = i
    tree = ET.parse(ip)
    root = tree.getroot()
    equiv_progs = []
    progId = 0
    with open(ip, "r") as in_f:
        for prog in root.findall('Program'):
            enodes += 1
            progId += 1
            hasJigOrBand = False
            equiv_assigns = []
            nEq = 0
            for equiv in prog.findall('equivalent'):
                enodes+=1
                nEq += 1
                lines = equiv.findall('line')
                line_assigns = []
                for line in lines:
                    if(line.text[0:6] == 'return'):
                        enodes-=1
                        continue
                    tool_split = line.text.split("=")
                    lhs1 = tool_split[0]
                    lhs_lst = lhs1.split(',')
                    varbs = []
                    if (len(lhs_lst) > 1):
                        fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                        lst = "(Var " + lhs_lst[len(lhs_lst) - 1].split(')')[0] + ")"
                        varbs.append(fst)
                        for i in range(1, len(lhs_lst) - 2):
                            varbs.append("(Var " + lhs_lst[i] + ")")
                        varbs.append(lst)
                        lhs = "(Assign ( Tup "
                        for v in varbs:
                            lhs = lhs + v
                        lhs = lhs + ")"
                    else:
                        # print("only one")
                        fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                        varbs.append(fst)
                        lhs = "(Assign ( Tup "
                        for v in varbs:
                            lhs = lhs + v
                        lhs = lhs #+ ")"
                    # print(lhs)
                    rhs = tool_split[1]
                    if "Chopsaw" in rhs:
                        assign = chopsaw (lhs, rhs)
                        line_assigns.append(assign)
                    elif "Bandsaw" in rhs:
                        assign = bandsaw (lhs, rhs)
                        line_assigns.append(assign)
                        hasJigOrBand = True
                    elif "Jigsaw" in rhs:
                        assign = jigsaw (lhs, rhs)
                        line_assigns.append(assign)
                        hasJigOrBand = True
                    elif "Tracksaw" in rhs:
                        assign = tracksaw(lhs, rhs)
                        line_assigns.append(assign)
                    elif "Drill" in rhs:
                        assign = drill (lhs, rhs)
                        line_assigns.append(assign)
                if len(line_assigns) != 0:
                    equiv_assigns.append(line_assigns)
            
        
            # print("Without Jig or Band ", progId - 1, nEq)
            #print("A", len(equiv_assigns))
            #print(equiv_assigns)
            #print('------------------------')
            tempSum =0
            for eqp in itertools.product(*equiv_assigns):
                tempSum+=len(eqp)
                #print(eqp)
                #enodes += 1
                equiv_progs.append(mk_sexp(eqp))
            #print("B", tempSum)
        return equiv_progs

def mk_sexps(d):
    os.chdir(d)
    for fnm in os.listdir("."):
        if fnm.endswith(".xml"):
            equiv_progs = parse_xml(fnm)
            nm = fnm.split('.')[0]
            dirnm = "../sub_sexps/" + nm
            if not os.path.exists(dirnm):
                os.mkdir(dirnm)
            os.chdir(dirnm)
            for i in range(len(equiv_progs)):
                op = nm + "_" + str(i) + ".sexp"
                with open(op, "w") as out_f:
                    out_f.write(equiv_progs[i])
            op_pick = nm + "_" + "pick.sexp"
            with open(op_pick, "w") as out_f:
                out_f.write("(Pick " + nm + ")")
            os.chdir("../../pre_processed_sub")
        else:
            continue

#mk_sexps("../benchmarks/bookcase/pre_processed_sub")

def mk_xml (all_progs, op):
    global enodes
    r = ET.Element('root')
    enodes = 0
    for module in all_progs:
        enodes += 1
        (idx, progs) = module
        sid = idx
        sub = ET.SubElement(r, 'sub', name = sid)
        for p in progs:
            pid = str(progs.index(p))
            ET.SubElement(sub, 'prog', name = pid).text = p
            enodes += 1
    tree = ET.ElementTree(r)
    tree.write(op)

def unify_sub_progs (d, op):
    module_progs = []
    for fnm in os.listdir(d):
        if fnm.endswith(".xml"):
            #TODO: does not contain the Pick programs
            equiv_progs = parse_xml(d + fnm)
            idx = fnm.split('.')[0]
            module_progs.append((idx, equiv_progs))
        else:
            continue
    mk_xml(module_progs, op)

# unify_sub_progs("../benchmarks/bookcase/193/", "../benchmarks/bookcase/collapse/all_progs.xml")
#unify_sub_progs("../benchmarks/bookcase-random/programs/", "../benchmarks/bookcase-random/collapse/all_progs.xml")
#unify_sub_progs("../benchmarks/bookcase/pre_processed_sub/", "../benchmarks/bookcase/collapse/all_progs.xml")
#unify_sub_progs("../benchmarks/bookcase/sub_sexps/", "../benchmarks/bookcase/collapse/all_progs.xml")
unify_sub_progs("../benchmarks/achair-rb-r9/programs/", "./all_progs.xml")
print(enodes)

def mk_top_prog(sub):
    return ("(Pick " + str(sub)  + ")")

# top is a number, subs is a list of numbers
def mk_tops(top, subs):
    progs = []
    op = "../benchmarks/bookcase-random/prog_sexps/" + top + ".sexp"
    with open (op, "w") as out_f:
        tops = []
        if len(subs) >= 3:
            for x in range(len(subs) - 2):
                tops.append ("(Seq " + mk_top_prog(subs[x]))
            final_seq = "(Seq " + mk_top_prog(subs[len(subs) - 2]) + mk_top_prog(subs[len(subs) - 1])
            num_close_paren = len(subs) - 1
            close_parens = ')' * num_close_paren
            out_f.write(''.join(tops) + final_seq + close_parens)
        elif len(subs) == 2:
            final_seq = "(Seq " + mk_top_prog(subs[0]) + mk_top_prog(subs[1]) + ")"
            out_f.write(final_seq)
        elif len(subs) == 1:
            out_f.write (mk_top_prog(subs[0]))
        else:
            out_f.write("(Empty)")

def parse_prog(d):
    for fnm in os.listdir(d):
        if fnm.endswith(".xml"):
            ip = "../benchmarks/bookcase-random/subs_progs_only_progs/" + fnm
            subs = open(ip).read().splitlines()
            nm = fnm.split('.')[0].split('_')[1]
            mk_tops(nm, subs)

#parse_prog("../benchmarks/bookcase/subs_progs_only_progs/")
