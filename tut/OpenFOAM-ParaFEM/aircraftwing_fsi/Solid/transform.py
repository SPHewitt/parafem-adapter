#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

def parse_arguments():
    '''
    Argument parser
    :return: '''
    parser = argparse.ArgumentParser(description='Script to transform the nodal coordinates')    
    parser.add_argument('-i',metavar='--input',type=str,nargs=1,dest='infile',help='input file name')
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
    if len(sys.argv) < 3:
        print("  Too few command line arguments have been set\n")
        exit(1)

    if file_exists(args.infile[0] + '.d'):
        print("  " + args.infile[0] + ".d not found\n")
        exit(1)

    dfile = str(args.infile[0] + '.d')

    print("  Error checks have been passed, transforming coordinates\n")

    # Tranformation matrix
    cos90 = np.cos(np.pi/2.0)
    sin90 = np.sin(np.pi/2.0)
    
    # Rotation about X
    rx = np.array([[1, 0, 0], [0, cos90, -sin90], [0, sin90, cos90]])
   
    # Read d file line by line
    data = list()
    with open(dfile) as file:
        for line in file:
            line = line.split()
            data.append(line)

    # transform node
    # Find Nodes
    nodeFlag = False
    for i in range(len(data)):
        if data[i] == ['*NODES']:
            nodeFlag = True
        elif data[i] == ["*ELEMENTS"]:
            nodeFlag = False
            break

        if (nodeFlag == True) and (data[i] != ['*NODES']):
            # Create temporary vector from data
            tmp = [float(data[i][1]),  float(data[i][2]), float(data[i][3])]
            tmp = np.matmul(rx, tmp)
            data[i][1] = str(tmp[0])
            data[i][2] = str(tmp[1])
            data[i][3] = str(tmp[2])

    # rewrite dfile
    file = open(dfile, "w+")
    for i in range(len(data)):
        string = str()
        for j in range(len(data[i])):
            string = string + ' ' + data[i][j]
        file.write(string+'\n')
    file.close()

    print("  Writing new file completed successfully\n")

