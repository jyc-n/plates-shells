"""
    Python Command File to Run Multiple Simulations
"""

#######################################################
import os
import sys
import time
import subprocess
import multiprocessing as mp
#######################################################

MAX_THREADS = 24

def createPath(case_id, opt_name=""):
    # get current work directory
    pwd = os.getcwd()

    # check if result folder exists
    result_folder = pwd + "/results/"
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)

    # check if operation folder exists
    opt_folder = result_folder
    if opt_name is not "":
        opt_folder = opt_folder + opt_name + "/"
        if not os.path.exists(opt_folder):
            os.mkdir(opt_folder)

    # check if case folder exists
    case_path = opt_folder + "Job-" + case_id + "/"
    if not os.path.exists(case_path):
        os.mkdir(case_path)

    # create soft link
    src = pwd + "/plates_shells"
    dst = case_path + "/plates_shells"
    if not os.path.islink(dst):
        os.symlink(src, dst)
    return case_path


def runSim(path, flag, jobName, cmd1, cmd2):
    pid = mp.current_process().name
    cmd = path + " " + flag + " " + jobName + " " + str(cmd1) + " " + str(cmd2)
    # print(pid)
    print(cmd)
    SimProcess = subprocess.Popen(cmd, shell=True)
    SimProcess.communicate()


def runMultiProcesses(flag, list1, list2):

    #* the following operations are in parallel
    pool = mp.Pool(processes=MAX_THREADS)
    for i, var1 in enumerate(list1):
        for j, var2 in enumerate(list2):
            case_id = str(i+1) + str(j+1)
            if (flag == "-p"):
                sim_type = "param"
            elif (flag == "-r"):
                sim_type = "refine"

            path = createPath(case_id, sim_type)
            jobName = "Job-" + case_id
            exe_file = path + "plates_shells"
            pool.apply_async(runSim, args=(exe_file, flag, jobName, var1, var2,))

    pool.close()
    pool.join()
    #* --------------------------- end parallel


def runRefine(flag, list_nodes):

    #* the following operations are in parallel
    pool = mp.Pool(processes=MAX_THREADS)

    for i, var1 in enumerate(list_nodes):
        var2 = var1
        case_id = str(i+1) + str(i+1)
        
        path = createPath(case_id, "refine")
        jobName = "Job-" + case_id
        exe_file = path + "plates_shells"
        pool.apply_async(runSim, args=(exe_file, flag, jobName, var1, var2,))

    pool.close()
    pool.join()
    #* --------------------------- end parallel

def main():
    args = sys.argv
    if len(args) != 2:
        sys.exit("Must give an argument!")

    try:
        cmd = args[-1]

        # check if compiled
        exefile = os.getcwd() + "/plates_shells"
        if not os.path.exists(exefile):
            raise Exception("Haven't compiled yet!!")

        # run single simulation
        if cmd == "srun":
            flag = "-p"
            jobName = "Job-00"
            var1 = "1e6"
            var2 = "0.3"
            case_id = "00"

            path = createPath(case_id, "param")
            exe_file = path + "plates_shells"
            runSim(exe_file, flag, jobName, var1, var2)

        # run multiple simulations
        elif cmd == "param":
            flag = "-p"
            E = list( range(1,101) ) 
            # E.extend( list( range(10,100+1,10) ) )
            # E.extend( list( range(100,1000,100) ) )
            # E.extend( list( range(1000,10000+1,1000) ) )
            E = [i/3.0*1e6 for i in E]
            T = [0.03]
            runMultiProcesses(flag, E, T)

        # perform mesh refinement test
        elif cmd == "refine":
            flag = "-r"
            nSide = list( range(2,50+1) )
            runRefine(flag, nSide)

        else:
            raise Exception("Wrong Commands")

    except Exception as err:
        print(err.args[0])
        sys.exit("Exception caught, existing...")


if __name__ == "__main__":
    main()