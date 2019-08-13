import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

label_node    = False
label_element = False


def readGeo(case_folder, coord, conn, nn, nel):
    fcoord = case_folder+"result03000.txt"
    with open(fcoord, 'r') as fin:
        for i in range(nn):
            line = fin.readline()
            ixyz = np.fromstring(line, dtype=float, sep=' ')
            coord[i,:] = ixyz

    fconn = case_folder+"connectivity.txt"
    with open(fconn, 'r') as fin:
        for i in range(nel):
            line = fin.readline()
            iel = np.fromstring(line, dtype=int, sep=' ')
            conn[i,:] = iel[1:4]


if __name__ == "__main__":

    path = "/Users/chenjingyu/Downloads/refine/"

    nSide = np.arange(25, 75+1, 5, dtype=int)
    nSide = [25]
    for icase in nSide:
        nn = (icase+1) * icase
        nel = (icase)*2 * (icase-1)         #! for clamped case, 1 extra row of elements

        case_folder = path+str(icase)+"_"+str(icase)+"/"

        coord = np.empty((nn,3), dtype=float)
        conn = np.empty((nel,3), dtype=int)

        readGeo(case_folder, coord, conn, nn, nel)

        fig = plt.figure(facecolor='white')
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect(1)

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

        plt.show()
        # plt.savefig('mesh.png')