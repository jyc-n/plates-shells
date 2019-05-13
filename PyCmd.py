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
    # check if result folder exists
    result_folder = pwd + '/results/'
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    # check if case folder exists
    case_folder = str(int(len)) + '_' + str(int(wid))
    case_path = result_folder + case_folder
    if not os.path.exists(case_path):
        os.mkdir(case_path)
    # create soft link
    src = pwd + '/plates_shells'
    dst = case_path + '/plates_shells'
    if not os.path.islink(dst):
        os.symlink(src, dst)
    return case_path


def run(path, cmd1, cmd2):
    pid = mp.current_process().name
    cmd = path + ' ' + str(cmd1) + ' ' + str(cmd2)
    print(pid)
    print(cmd)
    SimProcess = subprocess.Popen(cmd, shell=True)
    SimProcess.communicate()
    time.sleep(1)

def runSingleProcess(nLen=30, nWid=30):
    path = createPath(nLen, nWid)
    exe_file = 'plates_shells'
    subprocess.call([path+'/'+exe_file])

def runMultiProcesses():
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

def runRefine():
    nLen = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
    nWid = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
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
    elif sys.argv[1] == 'srun':
        runSingleProcess()
    elif sys.argv[1] == 'mrun':
        runMultiProcesses()
    elif sys.argv[1] == 'refine':
        runRefine()
    else:
        print('wrong arguments')


if __name__ == "__main__":
    main()