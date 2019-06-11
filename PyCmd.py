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

MAX_THREADS = 24

def pre():
    print("compiling c++ code")
    subprocess.call(["cmake", "CMakeLists.txt"])
    subprocess.call(["make"])


def createPath(case_id, opt_name=''):
    # get current work directory
    pwd = os.getcwd()

    # check if result folder exists
    result_folder = pwd + '/results/'
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)

    # check if operation folder exists
    opt_folder = result_folder
    if opt_name is not '':
        opt_folder = opt_folder + opt_name + '/'
        if not os.path.exists(opt_folder):
            os.mkdir(opt_folder)

    # check if case folder exists
    case_path = opt_folder + case_id + '/'
    if not os.path.exists(case_path):
        os.mkdir(case_path)

    # create soft link
    src = pwd + '/plates_shells'
    dst = case_path + '/plates_shells'
    if not os.path.islink(dst):
        os.symlink(src, dst)
    return case_path


def runSim(path, cmd1, cmd2):
    pid = mp.current_process().name
    cmd = path + ' ' + str(cmd1) + ' ' + str(cmd2)
    print(pid)
    print(cmd)
    SimProcess = subprocess.Popen(cmd, shell=True)
    SimProcess.communicate()


def runSingleProcess(nLen, nWid):
    case_id = str(int(nLen)) + '_' + str(int(nWid))
    path = createPath(case_id)
    exe_file = path + 'plates_shells'
    runSim(exe_file, nLen, nWid)


def runMultiProcesses():
    list_nLen = [5, 10, 15, 20]
    list_nWid = [15, 15, 15, 15]
    nProcesses = len(list_nLen)

    #* the following operations are in parallel
    pool = mp.Pool(processes=MAX_THREADS)
    for i in range(nProcesses):
        nLen = list_nLen[i]
        nWid = list_nWid[i]
        case_id = str(int(nLen)) + '_' + str(int(nWid))
        path = createPath(case_id, 'multi')
        exe_file = path + '/plates_shells'
        pool.apply_async(runSim, args=(exe_file, nLen, nWid,))

    pool.close()
    pool.join()
    #* --------------------------- end parallel


def runRefine(list_nLen, list_nWid):
    if len(list_nLen) != len(list_nWid):
        raise Exception("Target lists for refinement tests has different length")

    nCases = len(list_nLen)

    #* the following operations are in parallel
    pool = mp.Pool(processes=MAX_THREADS)

    for i in range(nCases):
        nLen = list_nLen[i]
        nWid = list_nWid[i]
        case_id = str(int(nLen)) + '_' + str(int(nWid))
        path = createPath(case_id, 'refine')
        exe_file = path + '/plates_shells'
        pool.apply_async(runSim, args=(exe_file, nLen, nWid,))

    pool.close()
    pool.join()
    #* --------------------------- end parallel


if __name__ == "__main__":
    args = sys.argv
    if len(args) != 2:
        sys.exit("Must give an argument!")

    try:
        cmd = args[-1]

        # compile code
        if cmd == 'pre':
            pre()

        # run single simulation
        elif cmd == 'srun':
            nLen = 20
            nWid = 20
            runSingleProcess(nLen, nWid)

        # run multiple simulations
        elif cmd == 'mrun':
            runMultiProcesses()

        # perform mesh refinement test
        elif cmd == 'refine':
            # list_nSide = [10]
            # list_nSide = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
            list_nSide = [21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76]
            runRefine(list_nSide, list_nSide)

        else:
            raise Exception("Wrong Commands")

    except Exception as err:
        print(err.args[0])
        sys.exit("Exception caught, existing...")