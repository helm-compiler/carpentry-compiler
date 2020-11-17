import re
import os
import sys
uuid_re = "[0-9a-fA-F]{8}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{4}\-[0-9a-fA-F]{12}"
uuid_cre = re.compile(uuid_re)

stock_library = {
    "b024b206-900c-4e9a-a934-00293b4bd186": "lumber_2x8x96",
    "aa4af7be-294c-4d35-b5ce-9a75861d60cd": "lumber_2x8x48",
    "652f5cd5-9e6e-48fa-b3f8-176f10fcdb56": "lumber_2x8x24",
    "764917a9-4b5a-4e52-a402-5cdfc7b0d055": "lumber_4x4x24",
    "906eb3e2-1df2-4764-a150-e736a3a37483": "lumber_4x4x48",
    "7eaa6d88-3d87-4afa-a5f5-2d0cae939480": "lumber_4x4x96",
    "b054fd41-d023-4cde-b5ca-cf6356ec6b14": "lumber_2x4x24",
    "32b733c8-9bb2-4c2b-b3af-4844992d5807": "lumber_2x4x48",
    "805da52d-c8e4-4ec6-8cfe-8a304f52ff24": "lumber_2x4x96",
    "8538f8db-adb6-4b5d-8ec0-6d6bad325ac4": "wood_0.75x12x20",
    "786c393f-b7fa-4910-9fc7-1492021f1a77": "wood_0.75x24x20",
    "69a58769-b39e-4267-a28e-aa8a09a93669": "wood_0.75x48x96",
    "af330a8f-db6b-46ee-9427-2af4ac227167": "wood_0.5x12x20",
    "44335d53-2916-4f71-81b3-2a5d76d71fff": "wood_0.5x24x20",
    "9d786f94-62a6-4e98-9c9c-0a2a2c9c27dd": "wood_0.5x48x96",
    "15411ad4-b2b1-4dbe-b5f3-f6e5df128d31": "lumber_2x2x24",
    "d01ecf45-d768-404c-8518-2dae6444cdff": "lumber_2x2x48",
    "1df1a8cb-d30e-4c7f-93ab-7de601a610ad": "lumber_2x2x96"
}


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


def dedup_part_uuids(uuids):
    return (list(dict.fromkeys(uuids)))


def mk_dict_from_uuids(uuids):
    dd_uuids = dedup_part_uuids(uuids)
    uuid_d = dict()
    ctr = 0
    for i in dd_uuids:
        uuid_d[i] = "a" + str(ctr)
        ctr = ctr + 1
    return uuid_d


def rename_parts(fnm, s, id_dict):
    fin = open(fnm, "r")
    for l in fin:
        ids = re.findall(uuid_cre, l)
        for i in ids:
            if i in stock_library:
                l = l.replace(i, stock_library[i])
            elif i in id_dict:
                l = l.replace(i, id_dict[i])
            else:
                continue
        s = s + l + "\n"
    return s


def renmame_parts_from_string(fnm, s, id_dict):
    ids = re.findall(uuid_cre, fnm)
    for i in ids:
        if i in stock_library:
            fnm = fnm.replace(i, stock_library[i])
        elif i in id_dict:
            fnm = fnm.replace(i, id_dict[i])
        else:
            continue
    s = s + fnm + "\n"
    return s


def get_return_parts(fnm, id_dict):
    fin = open(fnm, "r")
    stock_used = []
    for l in fin:
        l = l.split("=")[1]
        # print(cur)
        ids = re.findall(uuid_cre, l)
        for i in ids:
            if i in stock_library:
                l = l.replace(i, stock_library[i])
            elif i in id_dict:
                stock_used.append(i)
                l = l.replace(i, id_dict[i])
            else:
                continue
    stock_all = list(id_dict.keys())

    for s in stock_used:
        if s in stock_all:
            stock_all.remove(s)

    for s in stock_library.keys():
        if s in stock_all:
            stock_all.remove(s)

    s = "Return("
    for i in range(0, len(stock_all)):
        if i != len(stock_all) - 1:
            s += str(stock_all[i]) + ','
        else:
            s += str(stock_all[i]) + ')'
    return s


def rename_faces(s):
    s = s.rstrip()
    res = ""
    ss = s.split("\n")
    ctr = 0
    for l in ss:
        if (not ("Bandsaw" in l)) and (not ("Jigsaw" in l)):
            r = re.search(uuid_cre, l)
            if r != None:
                fid = r.group(0)
                l1 = l.replace(fid, "face_" + str(ctr))
                ctr = ctr + 1
                res = res + l1 + "\n"
            else:
                res = res + l + "\n"
        else:
            res = res + l + "\n"
    return res


def rename_edges(s):
    s = s.rstrip()
    res = ""
    ss = s.split("\n")
    ctr = 0
    for l in ss:
        rs = re.findall(uuid_cre, l)
        print (len(rs))
        if rs != []:
            for r in rs:
                l = l.replace(r, "edge_" + str(ctr))
                ctr = ctr + 1
            res = res + l + "\n"
        else:
            res = res + l + "\n"
    return res


def remove_last_arg(st):
    st1 = st.rstrip()
    ls = st1.split("\n")
    res = ""
    for l in ls:
        eq = l.split("=")
        if (len(eq) > 1):
            rhs = eq[1]
            if ("Chopsaw" in rhs) or ("Tracksaw" in rhs) or ("Drill" in rhs) or ("Jigsaw" in rhs) or ("Bandsaw" in rhs):
                cs = rhs.split(",")
                idx = len(cs) - 1
                l1 = l.replace("," + cs[idx], ")")
                res = res + l1 + "\n"
            else:
                res = res + l + "\n"
        else:
            res = res + l
    return res


def get_setup(l):
    ref_idx = l.find("Ref")
    ref_ss = l[ref_idx: len(l) - 1]
    return ref_ss

def get_band_setup(l):
    ref_idx = l.find("Path")
    ref_ss = l[ref_idx: len(l) - 1]
    return ref_ss

def get_new_chop_setup(l):
    ref_idx = l.find("Ref")
    ref_ss = l[ref_idx: len(l) - 1]
    rs = ref_ss.split("f")
    return "Setup_Chopsaw" + rs[1]

def get_new_drill_setup(l):
    ref_idx = l.find("Ref")
    ref_ss = l[ref_idx: len(l) - 1]
    rs = ref_ss.split("f")
    return "Setup_Drill" + rs[1]


def get_new_jig_setup(l):
    ref_idx = l.find("Ref")
    ref_ss = l[ref_idx: len(l) - 1]
    rs = ref_ss.split("f")
    return "Setup_Jigsaw" + rs[1]


def get_new_band_setup(l):
    ref_idx = l.find("Path")
    ref_ss = l[ref_idx: len(l) - 1]
    return "Setup_Bandsaw" + "(" + ref_ss + ")"


def get_new_track_setup(l):
    ref_idx = l.find("Ref")
    ref_ss = l[ref_idx: len(l) - 1]
    rs = ref_ss.split("f")
    return "Setup_Tracksaw" + rs[1]


def separate_setup(s):
    s = s.rstrip()
    res = ""
    ls = s.split("\n")
    for l in ls:
        if "Chopsaw" in l:
            setup = get_setup(l)
            n_setup = get_new_chop_setup(l)
            l1 = n_setup + "\n" + l.replace(", " + setup, "")
            res = res + l1 + "\n"
#        elif "Bandsaw" in l:
#            setup = get_band_setup(l)
#            n_setup = get_new_band_setup(l)
#            l1 = n_setup + "\n" + l.replace(", " + setup, "")
#            res = res + l1 + "\n"
#        elif "Jigsaw" in l:
#            setup = get_setup(l)
#            n_setup = get_new_jig_setup(l)
#            l1 = n_setup + "\n" + l.replace(", " + setup, "")
#            res = res + l1 + "\n"
        elif "Tracksaw" in l:
            setup = get_setup(l)
            n_setup = get_new_track_setup(l)
            l1 = n_setup + "\n" + l.replace(", " + setup, "")
            res = res + l1 + "\n"
        elif "Drill" in l:
            setup = get_setup(l)
            n_setup = get_new_track_setup(l)
            l1 = n_setup + "\n" + l.replace(", " + setup, "")
            res = res + l1 + "\n"
        else:
            res = res + l + "\n"
    return res

def remove_path(s):
    s = s.rstrip()
    res = ""
    ls = s.split("\n")
    for l in ls:
        if ("Bandsaw" in l) or ("Jigsaw" in l):
            path_idx = l.index("Path")
            ref_idx = l.index("Ref")
            remove = l[path_idx : ref_idx]
            l = l.replace(remove, "")
            res = res + l + "\n"
        else:
            res = res + l + "\n"
    return res

def dedup_setup(s):
    s = s.rstrip()
    res = ""
    curr = ""
    prev = ""
    ls = s.split("\n")
    for l in ls:
        if "Setup_" in l:
            curr = l
            if curr == prev:
                prev = curr
                res = res + ""
            else:
                prev = curr
                res = res + l + "\n"
        else:
            res = res + l + "\n"
    return res


def sugar_llhelm(fnm):
    fin = open(fnm, "r")
    fout = open("sugared-" + fnm, "w")
    ids = []
    for l in fin:
        ids = ids + get_all_part_uuids(l)
    id_dict = mk_dict_from_uuids(ids)
    rnm_ps = rename_parts(fnm, "", id_dict)
    ret_ps = get_return_parts(fnm, id_dict)
    rnm_fs = rename_faces(rnm_ps)
    rnm_es = rename_edges(rnm_fs)
    rmv_pt = remove_path(rnm_es)
    # the line below removes last argument from saw and drill
    clean_last = remove_last_arg(rmv_pt)
    #sep_set = separate_setup(rnm_es)
    sep_set = separate_setup(clean_last)
    dedup_set = dedup_setup(sep_set)

    fout.write(dedup_set)
    ret_ps = renmame_parts_from_string(ret_ps, "", id_dict)
    fout.write(ret_ps)
    fout.close()


if __name__ == "__main__":
    sugar_llhelm(sys.argv[1])
