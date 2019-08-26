#############################################
import os
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#############################################

label_node    = False
label_element = False

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

# plot a given case
def plotCase(icase, last_step_num):
    nn = icase * icase
    nel = (icase-1)*2 * (icase-1) 

    #! for clamped case, 1 extra row of elements
    # nn = (icase+1) * icase
    # nel = (icase)*2 * (icase-1)         

    case_folder = path+str(icase)+"_"+str(icase)+"/"

    coord = np.empty((nn,3), dtype=float)
    conn = np.empty((nel,3), dtype=int)

    readConn(case_folder, conn, nel)

    step_size = 10
    iframe = 1
    for istep in range(step_size, last_step_num+1, step_size):
        frame_name = "{:0>5d}".format(istep)
        frame_file = "result"+frame_name+".txt"
        readCoord(case_folder, frame_file, coord, nn)
        plotShape(case_folder, coord, conn, nn, nel, iframe, frame_name)
        iframe = iframe+1
    
    # create videos
    subprocess.call(["ffmpeg", "-i", case_folder+"frame%05d.png", case_folder+"out.mp4"])
    print("video generated")


# main
if __name__ == "__main__":
    path = os.getcwd()+"/"

    # nSide = np.arange(25, 75+1, 5, dtype=int)
    nSide = [10]
    last_step_num = 500
    for icase in nSide:
        plotCase(icase, last_step_num)
