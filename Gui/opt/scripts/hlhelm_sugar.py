import re
import sys

indent = 2

def space_after_comma (l):
    if l.group(0) == ',':
        return ', '
    else:
        return ''

def remove_keep (l):
    if (l.group(0) == ', Keep=1' or l.group(0) == ', Keep=0'):
        return ''
    else:
        return ''

def rename_query (l):
    if l.group(0) == 'Query_Face_By_Closest_Point':
        return 'Query'
    else:
        return ''

def space (n):
    return " " * n

def indent_sketch_line (l):
    if "Make_Sketch" in l:
        idx_sketch = l.find("Make_Sketch")
        idx_query = l.find("Query")
        idx_geom = l.find("Geom")
        idx_constr = l.find("Constraint")

        ss_sketch = l[idx_sketch : idx_query]
        ss_query = l[idx_query : idx_geom - 1]
        ss_geom = l[idx_geom : idx_constr - 1]
        ss_constr = l[idx_constr : ]
        ss_rest = l[0 : idx_sketch]

        i1 = ss_rest + ss_sketch
        idx1_sketch = i1.find("Make_Sketch")
        i2 = i1 + "\n" + space(indent + idx1_sketch) + ss_query
        i3 = i2 + "\n" + space(indent + idx1_sketch) + ss_geom
        i4 = i3 + "\n" + space(indent + idx1_sketch) + ss_constr
        return i4
    else:
        return l


def sugarify (l):
    l1 = re.sub(',', space_after_comma, l)
    l2 = re.sub('Query_Face_By_Closest_Point', rename_query, l1)
    l3 = re.sub(', Keep=0', rename_query, l2)
    l4 = re.sub(', Keep=1', rename_query, l3)
    l5 = indent_sketch_line (l4)
    return l5

def sugar_hlhelm (fnm):
    fin = open(fnm, "r")
    fout = open("sugared-" + fnm, "w+")
    sugared = ""
    for l in fin:
        new_l = sugarify(l)
        sugared = sugared + new_l
    fout.write(sugared)
    fout.close()

sugar_hlhelm (sys.argv[1])
