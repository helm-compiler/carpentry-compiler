import xml.etree.ElementTree as ET
import numpy as np
import sys
import uuid
import re
import os
import itertools
import shutil

uuid_re = "[0-9a-fA-F]{8}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{12}"
uuid_cre = re.compile(uuid_re)

stock_library = {
"b024b206-900c-4e9a-a934-00293b4bd186" : "lumber_2x8x96",
"aa4af7be-294c-4d35-b5ce-9a75861d60cd" : "lumber_2x8x48", 
"652f5cd5-9e6e-48fa-b3f8-176f10fcdb56" : "lumber_2x8x24",
"764917a9-4b5a-4e52-a402-5cdfc7b0d055" : "lumber_4x4x24",
"906eb3e2-1df2-4764-a150-e736a3a37483" : "lumber_4x4x48",
"7eaa6d88-3d87-4afa-a5f5-2d0cae939480" : "lumber_4x4x96",
"b054fd41-d023-4cde-b5ca-cf6356ec6b14" : "lumber_2x4x24",
"32b733c8-9bb2-4c2b-b3af-4844992d5807" : "lumber_2x4x48",
"805da52d-c8e4-4ec6-8cfe-8a304f52ff24" : "lumber_2x4x96",
"8538f8db-adb6-4b5d-8ec0-6d6bad325ac4" : "wood_0.75x12x20",
"786c393f-b7fa-4910-9fc7-1492021f1a77" : "wood_0.75x24x20",
"69a58769-b39e-4267-a28e-aa8a09a93669" : "wood_0.75x48x96",
"af330a8f-db6b-46ee-9427-2af4ac227167" : "wood_0.5x12x20",
"44335d53-2916-4f71-81b3-2a5d76d71fff" : "wood_0.5x24x20",
"9d786f94-62a6-4e98-9c9c-0a2a2c9c27dd" : "wood_0.5x48x96",
# don't know the meaning
"15411ad4-b2b1-4dbe-b5f3-f6e5df128d31" : "0.3",
"d01ecf45-d768-404c-8518-2dae6444cdff": "0.55",
"1df1a8cb-d30e-4c7f-93ab-7de601a610ad": "1.0",
"e49f6c42-16fe-4b9c-8198-0469f33182da" : "0.55",
"aac77b14-f1cf-4950-b48c-c4d9e012fddb": "1.0"
}


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
    ids = re.findall(uuid_cre, rhs)
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
    ids = re.findall(uuid_cre, rhs)
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
    ids = re.findall(uuid_cre, rhs)
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
    ids = re.findall(uuid_cre, rhs)
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

def process_xml(prog):
    parsed_file = []
    original_file = []
    prog_id = []
    lines = prog.split('\n')
    for line in lines:
        if(line[0:6] == 'return'):
            continue
        assign = ''
        tool_split = line.split("=")
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
        parsed_file.append(assign)
        original_file.append(line)


    return parsed_file, original_file

def collapse_xml(all_progs, op):
    global statENodes
    r = ET.Element('root')
    statENodes = 0
    for module in all_progs:
        (idx, progs, oris, prog_id) = module
        sid = idx
        sub = ET.SubElement(r, 'sub', name=sid)
        for id in range(0, len(progs)):
            prog = ET.SubElement(sub, 'prog', name=prog_id[id])
            statENodes += 1
            for j in range(0, len(progs[id])):
                statENodes += 1
                equiv = ET.SubElement(prog, 'equiv')
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
            equiv_progs, ori_progs, prog_ids = process_xml(d + fnm, '')
            idx = fnm.split('.')[0]
            print(len(equiv_progs))
            module_progs.append((idx, equiv_progs, ori_progs, prog_ids))
            # module_oriprogs.append((idx, original_progs))
        else:
            continue
    collapse_xml(module_progs, op)
    #collapse_xml(module_oriprogs, oop)

def dedup_part_uuids (uuids):
    return (list(dict.fromkeys(uuids)))

def get_all_part_uuids(l):
    uuids = []
    ss = l.split("=")
    if (len(ss) == 2):
        lhs = ss[0]
        rhs = ss[1]
        ids = re.findall(uuid_cre, lhs)
        for i in ids:
            uuids.append(i)
        id_from_rhs = re.search(uuid_cre, rhs).group(0)
        uuids.append(id_from_rhs)
        return uuids
    else:
        lhs = ss[0]
        ids = re.findall(uuid_cre, lhs)
        for i in ids:
            uuids.append(i)
        return uuids

def parse_results(fnm):
    outputs = []
    tree = ET.parse(fnm)
    root = tree.getroot()
    for indiv in root.findall('Individual'):
        output = indiv.find('Output').text
        outputs.append(output)
    return outputs

def mk_dict_from_uuids (uuids):
    dd_uuids = dedup_part_uuids (uuids)
    uuid_d = dict()
    for i in dd_uuids:
        uuid_d[i] = uuid.uuid4().__str__()
    return uuid_d

def update_parts(prog, id_dict):
    s = ""
    for l in prog:
        ids = re.findall(uuid_cre, l)
        for i in ids:
            if i in stock_library:
                l = l
            elif i in id_dict:
                l = l.replace(i, id_dict[i])
            else:
                continue
        s = s + l + "\n"
    return s[0:-1]

def update_prog(prog):
    ids = []
    lines = prog.lstrip().rstrip().split('\n')
    for l in lines:
        ids = ids + get_all_part_uuids(l)
    id_dict = mk_dict_from_uuids(ids)
    s = update_parts(lines, id_dict)
    return s


if __name__ == "__main__":
    schairs = parse_results("schair-result-processed.xml")
    coffeetables = parse_results("coffeetable-result-processed.xml")
    cnt = 0
    fileId = 0
    r = ET.Element('root')
    nSchairs = len(schairs)
    print(nSchairs)
    for p1 in range(0, nSchairs):
        for p2 in range(0, nSchairs):
            if p2 < p1:
                continue
            for p3 in range(0, nSchairs):
                if p3 < p2:
                    continue
                for p7 in coffeetables:
                    s = ""
                    s = update_prog(schairs[p1]) + "\n" + update_prog(schairs[p2]) + "\n" +  update_prog(schairs[p3]) + "\n" +  update_prog(p7)
                    parsed_file, original_file = process_xml(s)
                    sub = ET.SubElement(r, 'Program', name=str(cnt))
                    for line in parsed_file:
                        ET.SubElement(sub, 'sexp').text = line

                    #fn = open('test.txt', 'w+')
                    # fn.write(s)
                    # fn.close()
                    # print(s[0])
                    cnt += 1
                    print(cnt)

    tree = ET.ElementTree(r)
    tree.write("result"+str(fileId)+".xml")
    fileId += 1
    r.clear()
    exit(-1)


    for p1 in schairs:
        for p2 in schairs:
            for p3 in schairs:
                for p4 in schairs:
                    for p5 in schairs:
                        for p6 in schairs:
                            for p7 in coffeetables:
                                s = ""
                                s = update_prog(p1) + "\n" + update_prog(p2) + "\n" +  update_prog(p3) + "\n" +  update_prog(p4) + "\n" +  update_prog(p5) + "\n" +  update_prog(p6) + "\n" +  update_prog(p7)
                                parsed_file, original_file = process_xml(s)
                                sub = ET.SubElement(r, 'Program', name=str(cnt))
                                for line in parsed_file:
                                    ET.SubElement(sub, 'sexp').text = line

                                if cnt % 20000 == 0:
                                    tree = ET.ElementTree(r)
                                    tree.write("result"+str(fileId)+".xml")
                                    fileId += 1
                                    r.clear()
                                #fn = open('test.txt', 'w+')
                                # fn.write(s)
                                # fn.close()
                                # print(s[0])
                                cnt += 1
                                print(cnt)