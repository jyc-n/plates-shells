#############################################
import os
import sys
import subprocess
import multiprocessing as mp
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#############################################

label_node    = False
label_element = False
MAX_THREADS = 24

# read connectivity
def readConn(case_folder, conn, nel):
    fconn = case_folder+"connectivity.txt"
    with open(fconn, 'r') as fin:
        for i in range(nel):
            line = fin.readline()
            iel = np.fromstring(line, dtype=int, sep=' ')
            conn[i,:] = iel[1:4]

# read coordinate at a given frame
def readCoord(case_folder, frame_file, coord, nn):
    fcoord = case_folder + frame_file
    with open(fcoord, 'r') as fin:
        for i in range(nn):
            line = fin.readline()
            ixyz = np.fromstring(line, dtype=float, sep=' ')
            coord[i,:] = ixyz

# plot a given frame
def plotShape(case_folder, coord, conn, nn, nel, iframe, frame_name):
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    # ax.set_aspect(1)

    x = coord[:,0]
    y = coord[:,1]
    z = coord[:,2]

    for i in range(nel):
        iel = conn[i]
        index = np.array([iel[0]-1, iel[1]-1, iel[2]-1, iel[0]-1])
        el_x = x[index]
        el_y = y[index]
        el_z = z[index]

        # labeling elements
        if label_element:
            l_x = (x[iel[0]-1] + x[iel[1]-1] + x[iel[2]-1]) / 3.0
            l_y = (y[iel[0]-1] + y[iel[1]-1] + y[iel[2]-1]) / 3.0
            l_z = (z[iel[0]-1] + z[iel[1]-1] + z[iel[2]-1]) / 3.0
            label = str(i+1)
            ax.text(l_x, l_y, l_z, label, color='b')

        ax.plot(el_x, el_y, el_z, lw=0.5, c='k', ls='-')

    # labeling nodes
    if label_node:
        for i in range(nn):
            inn = coord[i]
            label = str(i+1)
            ax.text(inn[0], inn[1], inn[2], label, color='r')

    # First remove fill
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Now set color to white 
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    # plt.show()
    plt.title(frame_name)
    plt.savefig(case_folder+"frame{:0>5d}.png".format(iframe))
    plt.close()

# get the file id of the last file
def findLastFrame(case_folder):
    iframe = 0
    while iframe < 99999:
        frame_name = "{:0>5d}".format(iframe)
        frame_file = case_folder+"result"+frame_name+".txt"
        if not os.path.exists(frame_file):
            # iframe is the total number of frames, -1 to get the last file id
            return iframe-1
        else:
            iframe = iframe+1

# plot a given case
def plotCase(case_folder, nSide, last_frame_num, step_size=1):
    nn = nSide * nSide
    nel = (nSide-1)*2 * (nSide-1) 

    #! for clamped case, 1 extra row of elements
    # nn = (nSide+1) * nSide
    # nel = (nSide)*2 * (nSide-1)         

    coord = np.empty((nn,3), dtype=float)
    conn = np.empty((nel,3), dtype=int)

    readConn(case_folder, conn, nel)

    iframe = 1
    for istep in range(0, last_frame_num+1, step_size):
        frame_name = "{:0>5d}".format(istep)
        frame_file = "result"+frame_name+".txt"
        readCoord(case_folder, frame_file, coord, nn)
        plotShape(case_folder, coord, conn, nn, nel, iframe, frame_name)
        # print("frame {:d} plotted".format(iframe))
        iframe = iframe+1
    # print("All figures plotted!")

    # create videos
    subprocess.call(["ffmpeg", "-i", case_folder+"frame%05d.png", case_folder+"out.mp4"])
    # print("Video generated")

# plot multiple cases
def multiCasePlot(path, nSide, n_list1, n_list2):
    #* the following operations are in parallel
    pool = mp.Pool(processes=MAX_THREADS)
    for i in range(n_list1):
        for j in range(n_list2):
            jobName = "Job-"+str(i+1)+str(j+1)
            case_folder = path + jobName + "/"
            last_frame_num = findLastFrame(case_folder)

            last_frame_file = case_folder+"result"+"{:0>5d}".format(last_frame_num)+".txt"
            subprocess.call(["cp", last_frame_file, "./"+jobName+"_output.txt"])

            pool.apply_async(plotCase, args=(case_folder, nSide, last_frame_num,))

    pool.close()
    pool.join()
    #* --------------------------- end parallel


# project nodes onto the given bases
def projNodes(nvec1, nvec2, nvec3, node):
    t_node = np.array([np.dot(node,nvec1), np.dot(node,nvec2), np.dot(node,nvec3)], dtype=float)
    return t_node

# get projected coordinates
def projCoords(jobName, nside):
    myfile = jobName+"_output.txt"
    nn = nside * nside
    # read and save all coordinates
    coord = np.zeros([nn, 3], dtype = float)
    with open(myfile, 'r') as fin:
        for i in range(nn):
            line = fin.readline()
            xyz = np.fromstring(line, dtype=float, sep=' ')
            coord[i] = xyz

    # get coordinates of four corners
    #      3               4
    #      -----------------
    #      |               |
    #      |               |
    #      |               |
    #      -----------------
    #      1 (fixed)        2
    node1 = coord[0]
    node2 = coord[nside-1]
    node3 = coord[-nside]
    node4 = coord[-1]

    # dimensions vector measured between 1 and 4
    vecDist1 = node4 - node1
    nvec3 = np.array([0,0,1], dtype=float)

    # dimensions vector measured between 2 and 3
    vecDist2 = node2 - node3

    nvec2 = np.array([vecDist2[0], vecDist2[1], 0])
    nvec2 = nvec2 / np.linalg.norm(nvec2)
    
    nvec1 = np.cross(nvec3, nvec2)

    newfile = jobName+"_new.txt"
    with open(newfile, 'w') as fout:
        for i in range(nn):
            coord[i] = projNodes(nvec1, nvec2, nvec3, coord[i])
            line = "{:e} {:e} {:e}".format(coord[i,0], coord[i,1], coord[i,2])
            fout.write(line+"\n")


# main
def main():
    args = sys.argv
    if len(args) != 2:
        sys.exit("Must give an argument!")

    try:

        if args[-1] == "single":
            nSide = 30
            case_folder = os.getcwd()+"/results/hanging_param/Job-21/"
            last_frame_num = findLastFrame(case_folder)

            last_frame_file = case_folder+"result"+"{:0>5d}".format(last_frame_num)+".txt"
            # subprocess.call(["cp", last_frame_file, "./Job-00_output.txt"])

            plotCase(case_folder, nSide, last_frame_num)

        elif args[-1] == "all":
            nSide = 30
            E = list( range(1,10) ) 
            E.extend( list( range(10,100+1,10) ) )
            E = [i*1e5 for i in E]
            T = [0.3]
            path = os.getcwd()+"/results/param/"
            multiCasePlot(path, nSide, len(E), len(T))

        elif args[-1] == "last":
            # nSide = list( range(2,50+1) )
            E = list( range(1,101) ) 
            # E.extend( list( range(10,100+1,10) ) )
            # E.extend( list( range(100,1000,100) ) )
            # E.extend( list( range(1000,10000+1,1000) ) )
            E = [i*1e6 for i in E]
            path = os.getcwd()+"/results/param/"
            for i in range(len(E)):
                # jobName = "Job-"+str(i+1)+str(i+1)
                jobName = "Job-"+str(i+1)+"1"
                case_folder = path + jobName + "/"
                last_frame_num = findLastFrame(case_folder)

                last_frame_file = case_folder+"result"+"{:0>5d}".format(last_frame_num)+".txt"
                subprocess.call(["cp", last_frame_file, "./"+jobName+"_output.txt"])

        elif args[-1] == "proj":
            nSide = 30
            E = list( range(1,101) )
            E = [i*1e6 for i in E]
            path = os.getcwd()+"/"
            for i in range(len(E)):
                jobName = path+"Job-"+str(i+1)+"1"
                projCoords(jobName, nSide)

    except Exception as err:
        print(err.args[0])
        sys.exit("Exception caught, existing...")

if __name__ == "__main__":
    main()