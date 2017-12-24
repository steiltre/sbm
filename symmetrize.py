#! /usr/bin/env python

import sys

def symmetrize(pairs):
    lst = []
    for i in pairs:
        row = i[0]
        col = i[1]
        lst.append( (row,col) )
        lst.append( (col,row) )

    return set(lst)

def read_tsv(fname):
    lst = []
    with open(fname) as file:
        for line in file:
            vals = line.split("\t")
            row = vals[0]
            col = vals[1]
            lst.append( (row,col) )
    return lst

def write_tsv(fname, graph):
    file = open(fname, "w")
    for edge in graph:
        #file.write("{}\t{}\t1\n".format(edge[0]), edge[1])
        file.write(edge[0] + "\t" + edge[1] + "\t1\n")

if __name__ == "__main__":
    ifname = sys.argv[1]
    ofname = sys.argv[2]

    dir_graph = read_tsv(ifname)
    udir_graph = symmetrize(dir_graph)
    write_tsv(ofname, udir_graph)
