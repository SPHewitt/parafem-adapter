#!/usr/bin/env python

import os
import sys
import argparse


def parse_arguments():
    '''
    Argument parser
    :return:
    '''
    parser = argparse.ArgumentParser(description='Script to write an nset to an interface file (.int)')    
    parser.add_argument('-i',metavar='--input',type=str,nargs=1,dest='infile',help='input file name')
    parser.add_argument('-n',metavar='--nset',type=str,nargs=1,dest='nsetname',help='input nset name')
    return parser


def file_exists(filename):
    '''
    Routine checks if input filename exist in cwd
    :param filename:
    :return: bool
    '''
    file = True
    for i in os.listdir(os.getcwd()):
        if filename == i:
            file = False
            continue
    return file


def nset_exists(filename,nsetname):
    '''
    Check nsetname is present in file
    :param filename: input file name
    :param nsetname: nset name
    :return: bool
    '''
    nset = True
    with open(filename) as file:
        for line in file:
            line = line.split()
            if line[0] == "*NSET" and len(line) == 3:
                if line[2] == nsetname:
                    nset = False
                    break
    return nset


def get_node_list(filename, nsetname):
    '''
    Return a list of node for the setname
    :param filename:
    :param nsetname:
    :return:
    '''
    node_list = []
    nset_current = False

    with open(filename) as file:
        for line in file:
            line = line.split()

            if nset_current and (len(line) == 1):
                node_list.append(int(line[0]))
            elif nset_current and (len(line) > 1):
                break
            elif line[0] == "*NSET" and len(line) == 3:
                    if line[2] == nsetname:
                        nset_current = True
            else:
                continue

    return node_list


def get_node_coords(filename, nodelist):
    '''
    Get coordinates for the nodelist
    Please not this routine only works if the nodelist is ordered
    :param filename:
    :param nodelist:
    :return:
    '''
    xcoord = []
    ycoord = []
    zcoord = []
    node = []

    count = 0
    maxcount = len(nodelist)-1
    with open(filename) as file:
        for line in file:
            line = line.split()
            if count == maxcount:
                break
            if line[0] == str(nodelist[count]) and (len(line) == 4):
                node.append(int(line[0]))
                xcoord.append(float(line[1]))
                ycoord.append(float(line[2]))
                zcoord.append(float(line[3]))
                count = count + 1
            else:
                continue

    coords_ = {'node':node, 'x': xcoord, 'y': ycoord, 'z': zcoord}
    return coords_


if __name__ == "__main__":

    print("\n  -------nset2int.py-------\n")

    # Parse arguments
    args = parse_arguments().parse_args()

    # Error checking
    if len(sys.argv) < 4:
        print("  Too few command line arguments have been set\n")
        exit(1)

    if file_exists(args.infile[0] + '.nset'):
        print("  " + args.infile[0] + ".nset not found\n")
        exit(1)

    nsetfile = str(args.infile[0] + '.nset')
    dfile = str(args.infile[0] + '.d')
    intfile = str(args.infile[0] + '.int')

    if nset_exists(nsetfile, args.nsetname[0]):
        print("  Set name, " + args.nsetname[0] + ", not found in " + nsetfile + "\n")
        exit(1)

    nsetname = str(args.nsetname[0])

    print("  Error checks have been passed, creating"+intfile +" file\n")

    # Create lists of data for interface file:
    # nodes x-coordinate y-coordinate z-coordinate
    nodelist = get_node_list(nsetfile, nsetname)
    data = get_node_coords(dfile, nodelist)

    file = open(intfile, "w+")
    for i in range(len(data['node'])):
        file.write("{0}\t{1}\t{2}\t{3}\n".format(data['node'][i], data['x'][i], data['y'][i], data['z'][i]))
    file.close()

    print("  nset2int.py completed successfully\n")

