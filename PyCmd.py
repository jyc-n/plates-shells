'''
    Python Command File to Run Multiple Simulations
'''

#######################################################
import os
import sys
import time
import subprocess
import multiprocessing as mp
#######################################################

def pre():
    print("compiling c++ code")
    subprocess.call(["cmake", "CMakeLists.txt"])
    subprocess.call(["make"])
    

def createPath(len, wid):
    pwd = os.getcwd()
    case_folder = str(int(len)) + '_' + str(int(wid))
    path = pwd + '/results/' + case_folder
    if not os.path.exists(path):
        os.mkdir(path)
    src = pwd + '/plates_shells'
    dst = path + '/plates_shells'
    if not os.path.islink(dst):
        os.symlink(src, dst)
    return path


def run(path, cmd1, cmd2):
    pid = mp.current_process().name
    cmd = path + ' ' + str(cmd1) + ' ' + str(cmd2)
    print(pid)
    print(cmd)
    SimProcess = subprocess.Popen(cmd, shell=True)
    SimProcess.communicate()
    time.sleep(1)


def runProcess():
    nLen = [5, 10, 15, 20]
    nWid = [15, 15, 15, 15]
    nProcesses = len(nLen)

    pool = mp.Pool(processes=nProcesses)
    for i in range(nProcesses):
        cmd1 = nLen[i]
        cmd2 = nWid[i]
        path = createPath(cmd1, cmd2)
        exe_file = path + '/plates_shells'
        pool.apply_async(run, args=(exe_file, cmd1, cmd2,))

    pool.close()
    pool.join()    
        

def main():
    if sys.argv[1] == 'pre':
        pre()
    elif sys.argv[1] == 'run':
        runProcess()
    else:
        print('wrong arguments')


if __name__ == "__main__":
    main()