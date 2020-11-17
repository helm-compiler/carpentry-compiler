# -*- coding: utf-8 -*-
#!/usr/bin/python

import xml.etree.ElementTree as ET
import sys
import os
import re
import itertools
import shutil

# static variables
uuid_rgxp = r"[0-9a-fA-F]{8}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{12}"
statENodes = 0


def mk_sexp(assigns):
    sequify = []
    if len(assigns) >= 3:
        # print "3 assigns"
        # for a in assigns:
          #  print a
          #  print "\n"
        for x in range(len(assigns) - 2):
            sequify.append("(Seq " + assigns[x])
        final_seq = "(Seq " + assigns[len(assigns) -
                                      2] + assigns[len(assigns) - 1]
        num_close_paren = len(assigns) - 1
        close_parens = ')' * num_close_paren
        return (''.join(sequify) + final_seq + close_parens)
    elif len(assigns) == 2:
        #print("2 assigns")
        final_seq = "(Seq " + assigns[0] + assigns[1] + ")"
        return final_seq
    elif len(assigns) == 1:
        #print("1 assign")
        return assigns[0]
    else:
        return "(Empty)"


def chopsaw(lhs, rhs):
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


''' The version using paths as references
def mk_path(ps, pathType):
    path = []
    p1 = ps[0].split('(')[1]
    path.append(mk_xy(p1))
    # first and last are taken care of above
    for i in range(1, len(ps) - 2):
        path.append(mk_xy(ps[i]))
    pn = ps[len(ps) - 1].split(')')[0]
    path.append(mk_xy(pn))
    spath = "(" + pathType + " "
    for i in range(0, len(path)):
        spath = spath + path[i]
    return (spath + ")")

def bandsaw(lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    height = "(Height (Float "+rsplt[len(rsplt)-1].split(')')[0]+"))"
    rsplt = rhs.split('Ref((')[0]
    rsplt = rsplt.split(',')
    ps = []
    # last two list elements will always be about ref.
    for i in range(1, len(rsplt) - 1):
        ps.append(rsplt[i].strip(')'))
    path = mk_path(ps, "Refb")
    stk = "(Bool true)"
    rhs = "(Bandsaw " + lumber + path + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign


def jigsaw(lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    height = "(Height (Float "+rsplt[len(rsplt)-1].split(')')[0]+"))"
    rsplt = rhs.split('Ref((')[0]
    rsplt = rsplt.split(',')
    ps = []
    # last two list elements will always be about ref.
    for i in range(1, len(rsplt) - 1):
        ps.append(rsplt[i].strip(')'))
    path = mk_path(ps, "Refj")
    stk = "(Bool true)"
    rhs = "(Jigsaw " + lumber + path + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign
'''


def mk_ref(x, y):
    px = "(Float " + x + " )"
    py = "(Float " + y + " )"
    pt = "(Tup " + px + py + " )"
    return pt


def mk_path(ps, pathType):
    path = []
    p1 = ps[0].split('(')[1]
    path.append(mk_xy(p1))
    # first and last are taken care of above
    for i in range(1, len(ps) - 2):
        path.append(mk_xy(ps[i]))
    pn = ps[len(ps) - 1].split(')')[0]
    path.append(mk_xy(pn))
    spath = "(" + pathType + " "
    for i in range(0, len(path)):
        spath = spath + path[i]
    return (spath + ")")


def mk_refs(ps, refType):
    path = []
    for i in range(0, len(ps)):
        tmpPs = ps[i].split(',')
        path.append(mk_ref(tmpPs[1], tmpPs[3]))

    spath = "(" + refType + " "
    for i in range(0, len(path)):
        spath = spath + path[i]

    return (spath + ")")


def bandsaw(lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    height = "(Height (Float "+rsplt[len(rsplt)-1].split(')')[0]+"))"

    path = rhs.split('Ref((')[1].split('))')[0].split('), (')
    path = mk_refs(path, "Refb")

    stk = "(Bool true)"
    rhs = "(Bandsaw " + lumber + path + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign


def jigsaw(lhs, rhs):
    ids = re.findall(uuid_rgxp, rhs)
    lumber = "( Lumber " + ids[0] + ")"
    rsplt = rhs.split(',')
    height = "(Height (Float "+rsplt[len(rsplt)-1].split(')')[0]+"))"

    path = rhs.split('Ref((')[1].split('))')[0].split('), (')
    path = mk_refs(path, "Refj")

    stk = "(Bool true)"
    rhs = "(Jigsaw " + lumber + path + stk + height + ")"
    assign = lhs + rhs + " )"
    return assign


def drill(lhs, rhs):
    return "TODO_DRILL"


def tracksaw(lhs, rhs):
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


def parse_xml(i):
    ip = i
    tree = ET.parse(ip)
    root = tree.getroot()
    equiv_progs = []
    original_progs = []
    with open(ip, "r") as in_f:
        for prog in root.findall('Program'):
            equiv_assigns = []
            equiv_originals = []
            for equiv in prog.findall('equivalent'):
                lines = equiv.findall('line')
                line_assigns = []
                line_originals = []
                for line in lines:
                    if(line.text[0:6] == 'return'):
                        continue

                    tool_split = line.text.split("=")
                    lhs1 = tool_split[0]
                    lhs_lst = lhs1.split(',')
                    varbs = []
                    if (len(lhs_lst) > 1):
                        fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                        lst = "(Var " + \
                            lhs_lst[len(lhs_lst) - 1].split(')')[0] + ")"
                        varbs.append(fst)
                        for i in range(1, len(lhs_lst) - 1):
                            varbs.append("(Var " + lhs_lst[i] + ")")
                        varbs.append(lst)
                        lhs = "(Assign ( Tup "
                        for v in varbs:
                            lhs = lhs + v
                        lhs = lhs + ")"
                    else:
                        #print("only one")
                        fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                        varbs.append(fst)
                        lhs = "(Assign ( Tup "
                        for v in varbs:
                            lhs = lhs + v
                        lhs = lhs  # + ")"
                    # print(lhs)
                    rhs = tool_split[1]
                    if "Chopsaw" in rhs:
                        assign = chopsaw(lhs, rhs)
                        line_assigns.append(assign)
                        line_originals.append(line.text)
                        # print(assign)
                        # print(line.text)
                    elif "Bandsaw" in rhs:
                        assign = bandsaw(lhs, rhs)
                        line_assigns.append(assign)
                        line_originals.append(line.text)
                    elif "Jigsaw" in rhs:
                        assign = jigsaw(lhs, rhs)
                        line_assigns.append(assign)
                        line_originals.append(line.text)
                    elif "Tracksaw" in rhs:
                        assign = tracksaw(lhs, rhs)
                        line_assigns.append(assign)
                        line_originals.append(line.text)
                    elif "Drill" in rhs:
                        assign = drill(lhs, rhs)
                        line_assigns.append(assign)
                        line_originals.append(line.text)
                if len(line_assigns) != 0:
                    equiv_assigns.append(line_assigns)
                    equiv_originals.append(line_originals)
                    print(line_assigns[0])

            exit(-1)
            # for eqp in itertools.product(*equiv_assigns):
            #     equiv_progs.append(mk_sexp(eqp))

            # for eqp in itertools.product(*equiv_originals):
            #     original_progs.append(''.join(eqp))

        return equiv_progs, original_progs


def process_xml(i, o):
    ip = i
    tree = ET.parse(ip)
    root = tree.getroot()

    parsed_file = []
    original_file = []
    original_cutLineId_file = []
    prog_id = []
    prog_wlc = []
    prog_arr = []
    for prog in root.findall('Program'):
        parsed_prog = []
        original_prog = []
        original_cutLineId_prog = []
        
        for equiv in prog.findall('equivalent'):
            lines = equiv.findall('line')
            parsed_equiv = []
            original_equiv = []
            prog_equiv = []
            for line in lines:
                if(line.text[0:6] == 'return'):
                    continue
                assign = ''
                tool_split = line.text.split("=")
                lhs1 = tool_split[0]
                lhs_lst = lhs1.split(',')

                varbs = []
                if (len(lhs_lst) > 1):
                    fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                    lst = "(Var " + \
                        lhs_lst[len(lhs_lst) - 1].split(')')[0] + ")"
                    varbs.append(fst)
                    for i in range(1, len(lhs_lst) - 1):
                        varbs.append("(Var " + lhs_lst[i] + ")")
                    varbs.append(lst)
                    lhs = "(Assign ( Tup "
                    for v in varbs:
                        lhs = lhs + v
                    lhs = lhs + ")"
                else:
                    #print("only one")
                    fst = "(Var " + lhs_lst[0].split('(')[1] + ")"
                    varbs.append(fst)
                    lhs = "(Assign ( Tup "
                    for v in varbs:
                        lhs = lhs + v
                    lhs = lhs  # + ")"
                # print(lhs)
                rhs = tool_split[1]
                if "Chopsaw" in rhs:
                    assign = chopsaw(lhs, rhs)
                elif "Bandsaw" in rhs:
                    assign = bandsaw(lhs, rhs)
                elif "Jigsaw" in rhs:
                    assign = jigsaw(lhs, rhs)
                elif "Tracksaw" in rhs:
                    assign = tracksaw(lhs, rhs)
                elif "Drill" in rhs:
                    assign = drill(lhs, rhs)
                parsed_equiv.append(assign)
                original_equiv.append(line.text)
                # print(assign)
            if len(parsed_equiv) != 0:
                parsed_prog.append(parsed_equiv)
                original_prog.append(original_equiv)
                original_cutLineId_prog.append(equiv.attrib["cutLineId"])

        prog_id.append(prog.attrib["id"])
        prog_wlc.append(prog.attrib["wlc"])
        prog_arr.append(prog.attrib["arr"])
        parsed_file.append(parsed_prog)
        original_file.append(original_prog)
        original_cutLineId_file.append(original_cutLineId_prog)

    return parsed_file, original_file,original_cutLineId_file, prog_id, prog_wlc, prog_arr
    # for eqp in itertools.product(*equiv_assigns):
    #     equiv_progs.append(mk_sexp(eqp))

    # for eqp in itertools.product(*equiv_originals):
    #     original_progs.append(''.join(eqp))

def collapse_xml(all_progs, op):
    global statENodes
    r = ET.Element('root')
    statENodes = 0
    for module in all_progs:
        (idx, progs, oris, oris_cutId, prog_id,  prog_wlc, prog_arr) = module
        sid = idx
        sub = ET.SubElement(r, 'sub', eclassID=sid)
        for id in range(0, len(progs)):
            prog = ET.SubElement(sub, 'prog', enodeID=prog_id[id], wlcID = prog_wlc[id],  arrID = prog_arr[id])
            statENodes += 1
            for j in range(0, len(progs[id])):
                statENodes += 1
                equiv = ET.SubElement(prog, 'equiv', cutLineID = oris_cutId[id][j])
                for k in range(0, len(progs[id][j])):
                    statENodes += 1
                    line = ET.SubElement(equiv, 'line')
                    ET.SubElement(line, 'sexp').text = progs[id][j][k]
                    ET.SubElement(line, 'ori').text = oris[id][j][k]

    tree = ET.ElementTree(r)
    tree.write(op)


def unify_sub_progs(d, op, oop):
    module_progs = []
    module_oriprogs = []
    for fnm in os.listdir(d):
        if fnm.endswith(".xml"):
            # TODO: does not contain the Pick programs
            print(d+fnm)
            equiv_progs, ori_progs, ori_cutId_progs, prog_ids, prog_wlcs, prog_arrs = process_xml(d + fnm, '')
            idx = fnm.split('.')[0]
            print(len(equiv_progs))
            module_progs.append((idx, equiv_progs, ori_progs,ori_cutId_progs, prog_ids, prog_wlcs, prog_arrs))
            # module_oriprogs.append((idx, original_progs))
        else:
            continue
    collapse_xml(module_progs, op)
    #collapse_xml(module_oriprogs, oop)

def str_to_raw(s):
    raw_map = {8:r'\b', 7:r'\a', 12:r'\f', 10:r'\n', 13:r'\r', 9:r'\t', 11:r'\v'}
    return r''.join(i if ord(i) > 32 else raw_map.get(ord(i), i) for i in s)

def parse(folder):
    egraphGeneratedFolder = str_to_raw(folder)
    programFolder = egraphGeneratedFolder + "/programs/"
    collapsedFolder = egraphGeneratedFolder + "/collapse/"
    collapsedFile = egraphGeneratedFolder + "/collapse/all_progs.xml"
    originalFile = egraphGeneratedFolder + "/collapse/ori_progs.xml"

    if not os.path.exists(programFolder):
        os.makedirs(programFolder)
    else:
        print('error, folders already exist')
        exit(-2)
    
    if not os.path.exists(collapsedFolder):
        os.makedirs(collapsedFolder)
    else:
        print('error, folders already exist')
        exit(-2)

    fileIndices = []
    for fnm in os.listdir(egraphGeneratedFolder):
        if fnm.endswith(".xml"):
            if not fnm.endswith("egraph.xml"):
                fileIndices.append(int(fnm.split('.')[0]))
                shutil.move(egraphGeneratedFolder+"/"+fnm, programFolder+fnm)
    if not os.path.exists(programFolder):
        os.makedirs(programFolder)
    if not os.path.exists(collapsedFolder):
        os.makedirs(collapsedFolder)

    unify_sub_progs(programFolder, collapsedFile, originalFile) 
    print("#E-Nodes = " + str(statENodes))

def main():
    print(sys.argv[1])
    parse(sys.argv[1])

main()
