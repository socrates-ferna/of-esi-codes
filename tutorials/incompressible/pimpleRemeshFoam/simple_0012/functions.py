import os
import subprocess
from glob import glob
import shutil
import errno

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def get_latestTime(casePath,parallel):
    latestTime = 0
    if parallel:
        casePath = os.path.join(casePath,'processor0')
    for entry in os.listdir(casePath):
        try:
            entry = float(entry)
            print("Time dir: ",entry)
            latestTime = max(latestTime,entry)
        except: continue
    print("latestTime is: ",latestTime,'s')

    return latestTime


def run_solver(solver,parallel,casePath,logFileName,crescent):
    #logFileName = 'log.'+solver
    if casePath != os.getcwd():
        os.chdir(casePath)
        print("Changed dir to :", casePath)
        return ValueError
    
    if parallel: ## aquí debe ir también el flag de si estamos en crescent
        if os.path.exists('processor0'):
            print("case already decomposed. RESTARTING")
            for procname in glob('processor*'):
                shutil.rmtree(procname)
                subprocess.run("decomposePar",shell=True)
        else:
            subprocess.run("decomposePar",shell=True)
        #decompeParDictPath = os.path.join(os.getcwd(),'system/decomposeParDict')
        with open('system/decomposeParDict','r') as f:
            lines = f.readlines()
            nProcs = searchOFScalar(lines,'numberOfSubdomains')
            nProcs = int(nProcs)
        if crescent:
            print("CRESCENT RUN")
            crescentStr = ("mpirun -x LD_LIBRARY_PATH -x FOAM_USER_LIBBIN -x WM_PROJECT_USER_DIR -x PATH"
             " -x WM_PROJECT_DIR -machinefile $PBS_NODEFILE -np "
             + str(nProcs) + " " + solver + " -parallel >> log." + solver)
            print("Running the command: "+crescentStr)
            subprocess.run(crescentStr,shell=True)
            #en crescent no sé si controlas dónde se escribe el output
        else:
            logFile = open(logFileName,'w')
            subprocess.run(' '.join(['mpirun','-np',str(nProcs),solver,'-parallel']),shell=True,stdout=logFile)
            logFile.close()
    else:
        nProcs = 1
        logFile = open(logFileName,'w')
        subprocess.run(solver,shell=True,stdout=logFile)
        logFile.close()

    return nProcs

def searchOFScalar(iterable,searchStr,unique=True):
    if not unique:
        iterable = reversed(iterable)
    
    for line in iterable:
        if searchStr in line:
            return float(line.split()[-1].replace(';',''))
    return 0

def replace_line(filename,new,pos):
    with open(filename,'r') as f:
        filelines = f.readlines()
    
    filelines[pos] = new

    with open(filename,'w') as f:
        f.writelines(filelines)
    


def replaceOFScalars(filename,valuesDict,endstr=';\n',sep=' ',avoid=''):
    with open(filename,'r') as f:
        filelines = f.readlines()
    
    for pos,line in enumerate(filelines):
        for key,val in valuesDict.items():
            if key in line:
                if avoid and avoid in line:
                    continue
                else:
                    newline = line.replace(line.split(sep)[-1],str(val)+endstr)
                    filelines[pos] = newline

    with open(filename,'w') as f:
        f.writelines(filelines)




def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]
