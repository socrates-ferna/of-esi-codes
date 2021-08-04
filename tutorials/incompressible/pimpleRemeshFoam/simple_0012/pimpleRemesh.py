import os
#import sys
import subprocess
import shutil
import parameters as p
from glob import glob

from numpy.lib.function_base import angle
import functions
from numpy import pi

def main():
    
    """
    0. Estaría bien que hiciera la parte de inicialización steady
    1. Lanza pimpleFoam. Se para solo
    2. Coger tail log. Toda la info: error,I,D,omega,angle
    3. Caso nuevo: copia sólo carpetas necesarias + el 0 de otro lao
    4. sustituciones: controlDelay, error,I,D,omega,angle, (coger windowValuesy windowTimes tb)
    5. mapFields (o antes del 4?)
    6. mv unmapped. rm CellDisplacement
    """
    
    #print("REMEMBER TO SOURCE OPENFOAM BEFORE USING THIS SCRIPT")
    #dummy = input("Press any key to continue")
    namesList = ["Cl","omega","angle","error","errorIntegral","errorDifferential"]
    currentTime = 0.0
    PIDreplNameList = namesList[3:]
    PIDreplNameList.append("controlDelay")
    valueList = [0]*len(namesList) #EL TIEMPO LO VAS A COGER DE LAS CARPETAS ESCRITAS
    searchList = ["Controlled Variable","Saturated omega","Saturated angle","error  ","errorIntegral","errorDifferential"]
    valueDict = dict(zip(namesList,valueList))
    searchDict = dict(zip(namesList,searchList))
    parentPath = os.getcwd()
    print("Parent Path: ", parentPath)
    caseName = p.caseName
    origCaseName = p.caseName
    casePath = os.path.join(parentPath,caseName)
    print("Startup case: ", casePath)
    solver = p.solver
    parallel = p.parallel
    crescent = p.crescent
    meshType = p.meshType
    mapFieldsCmd = "mapFields -consistent -sourceTime latestTime"
    if parallel:
        mapFieldsCmd = mapFieldsCmd + " -parallelSource ../"
    else:
        mapFieldsCmd = mapFieldsCmd + " ../"

    if meshType == "airfoil":
        remeshCommand = p.AIRFOILREMESH
    elif meshType == "flap":
        remeshCommand = p.FLAPREMESH
        rotateCenterstr = ' '.join([str(i) for i in p.rotateCenter])
        print(rotateCenterstr, "ROTATION CENTER")
        functions.replaceOFScalars(origCaseName+'parameters.dat',{"rotateCenter": rotateCenterstr},endstr='\n',sep='=')
    else:
        ValueError
        return "INVALID MESH TYPE SCRIPT STOPPED"



    if not crescent:
        logFileName = 'log.' + solver
    else:
        logFileName = "dummy" ## pendiente. REQUIERE DE UN GLOB O ALGO QUE JUNTE EL NOMBRE DEL ACCESO DIRECTO CON EL NOMBRE DEL CASO

    #Get end Time
    with open(caseName+'/system/controlDict','r') as f:
        lines = f.readlines()
        endTime = functions.searchOFScalar(lines,'endTime',unique=False)
    
    ### AQUÍ DEBERÍA IR UN WHILE QUE FUERA CREANDO NUEVOS CASOS HASTA QUE SE TERMINA EL RUN
    
    #FIRST RUN
    os.chdir(casePath)
    subprocess.run(['./Allclean'],shell=True)
    caseFileList = os.listdir(casePath)
    with open('0/pointDisplacement','r') as f:
        lines = f.readlines()
        angle0 = functions.searchOFScalar(lines,'angle0')

    nProcs = functions.run_solver(solver,parallel,casePath,logFileName,crescent)

    with open('system/PIDcontrolDict','r') as f:
        lines = f.readlines()
        controlDelay = functions.searchOFScalar(lines,'controlDelay')
    #miro las carpetas de tiempo escritas
    currentTime = functions.get_latestTime(casePath,parallel)

    if not p.testing:
        if currentTime >= endTime:
            return "FINISHED RUN"
    else:
        endTime += 0.1
        functions.replaceOFScalars('system/controlDict',{'endTime': str(endTime)},avoid='stopAt')

    #start looping
    while currentTime < endTime:
        prevCasePath = casePath
        prevCaseName = caseName
        os.chdir(prevCasePath)

        #1. catch all info from previous run
        with open('system/PIDcontrolDict','r') as f:
                lines = f.readlines()
                controlDelay = functions.searchOFScalar(lines,'controlDelay')

        controlDelay = max(0,controlDelay-currentTime)

        f = open(logFileName)
        lastLines = functions.tail(f,lines=200) ## LOOK HOW MANY LINES YOU REALLY NEED
        f.close()
        
        for key,searchStr in searchDict.items():
            valueDict[key] = functions.searchOFScalar(lastLines,searchStr,unique=False)
        
        #2. Prepare "empty" case for next run
        caseName = '_'.join([origCaseName,"{:2f}".format(currentTime)])
        casePath = os.path.join(parentPath,caseName)
        os.chdir(parentPath)
        if os.path.exists(casePath):
            shutil.rmtree(casePath)

        os.mkdir(casePath)
        
        #Copy case without results
        for fname in caseFileList:
            tempPath = os.path.join(origCaseName,fname)
            tempDestPath = os.path.join(casePath,fname)
            if os.path.isdir(tempPath):
                shutil.copytree(tempPath,tempDestPath)
            elif os.path.isfile(tempPath):
                shutil.copyfile(tempPath,tempDestPath)            
            else:
                print("WARNING, NEITHER FILE NOR DIRECTORY")
            print("Copied: ",fname, " to destination: ", tempDestPath)
        
        #Copy uninitialized 0
        zeroPath = os.path.join(casePath,'0')
        templateZeroPath = os.path.join(parentPath,'0')
        shutil.rmtree(zeroPath)
        shutil.copytree(templateZeroPath,zeroPath)


        #3. Prepare startup, remesh, mapFields
        os.chdir(casePath)

        #PID controller to take up run from previous conditions
        PIDreplValuelist = [valueDict["error"],valueDict["errorIntegral"],valueDict["errorDifferential"],controlDelay]
        PIDreplDict = dict(zip(PIDreplNameList,PIDreplValuelist))
        functions.replaceOFScalars('system/PIDcontrolDict',PIDreplDict)


        angle0Dict = {"angle0": valueDict["angle"],"omega":valueDict["omega"]}
        parameterDict = {"rotateAngle": valueDict["angle"]*180/pi}
        functions.replaceOFScalars('parameters.dat',parameterDict,endstr='\n',sep='=')
        
        #Remesh call. Remember to have licence active
        subprocess.run(remeshCommand,shell=True)

        #mapFields with parallel support
        subprocess.run(mapFieldsCmd+prevCaseName,shell=True)
        shutil.copyfile('pointDisplacement','0/pointDisplacement')
        functions.silentremove('0/cellDisplacement')
        #os.rename("0/pointDisplacement.unmapped","0/pointDisplacement")
        os.chdir('0')
        subprocess.run("sed -i 's/nan/0.0/g' *",shell=True)
        os.chdir('../')
        subprocess.run("renumberMesh -overwrite",shell=True)
        """if parallel:
            for dir in glob('processor*'):
                functions.replaceOFScalars(dir+'/0/pointDisplacement')
        else:"""
        functions.replaceOFScalars('0/pointDisplacement',angle0Dict)
        #correr solver de nuevo
        if p.testing:
            endTime += 0.1
            functions.replaceOFScalars('system/controlDict',{"endTime": str(endTime)},avoid='stopAt')
            functions.replaceOFScalars('system/controlDict',{"writeInterval": 0.02})
        
        os.rename('0',str(currentTime))
        functions.run_solver(solver,parallel,os.path.join(parentPath,caseName),logFileName,crescent)

        currentTime = functions.get_latestTime(casePath,parallel)
        endTime += 0.1
        if currentTime >= endTime:
            print("FINISHED RUN")
        else:
            print("BAD QUALITY. BACK TO REMESHING")
    



    return 0

    """
    if parallel: ## aquí debe ir también el flag de si estamos en crescent
        subprocess.run(['decomposePar'])
        #decompeParDictPath = os.path.join(os.getcwd(),'system/decomposeParDict')
        with open('system/decomposeParDict','r') as f:
            lines = f.readlines()
            nProcs = functions.searchOFScalar(lines,'numberOfSubdomains')
            nProcs = int(nProcs)
        if crescent:
            print("CRESCENT RUN")
        logFile = open(logFileName,'w')
        subprocess.run(['mpirun','-np',nProcs,solver,'-parallel'],shell=True,stdout=logFile)
        logFile.close()
    else:
        logFile = open(logFileName,'w')
        subprocess.run([solver],shell=True,stdout=logFile)
        logFile.close()
    """

if __name__ == "__main__":
    main()