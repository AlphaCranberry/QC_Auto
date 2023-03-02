#!/usr/bin/env python3
#native imports that don't need installed
import os # operation system functions
import sys # system functions
import subprocess #subproccess systm package. similar to os functions
import time #native import capable of keeping track of time for me
import random as rand

#start tracking the time the program takes to run
timerStart = time.perf_counter()

#try to import all packages, install if not
isPass = False
while isPass == False:
    try:
        #import bunch of programs needed
        import re
        import multiprocessing #multiprocess package
        from multiprocessing.connection import wait
        import pandas as pd # gives excel like functions
        import seaborn as sns # pretty plot package
        from matplotlib import pyplot as plt # basic plotting package
        import glob
        import re
        import csv
        import gc
        import math
        from collections import Counter
        import argparse

        # configure where plots are displayed
        #%matplotlib inline 
        import numpy as np # numerical functions
        import scipy.stats as stats #statistical engine
        import matplotlib.lines as mlines # for midlines
        sys.setrecursionlimit(3000)
        
        # exit the loop
        isPass = True
    #split up the error data and grab just the module name, then install that module
    #with mamba... hopefully
    except ModuleNotFoundError:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        splitObjects = str(exception_object).split("'")
        
        #now try to install it
        os.system("mamba install "+splitObjects[1])


## here will be the functions block. define all functions to be used in the program

# homer trimming function
def trimFunc(trimFile):
    os.system(f'homerTools trim -3 AGATCGGAAGAGCACACGTCT -mis 2 -minMatchLength 4 -min 20 {trimFile}')

# a function for alinging to multiprocess
def hisat2Func(file, builtGenome, folder):
    output = file + '.sam'
    mapping_stats = output + '_mappingstats.txt'
    ## p.s. the built genome ref should be a hisat2 built genome directory pointing to the name of the
    ## file with the prefix of the hisat2 files.
    ## i.e. /data/lab/duttke/genome/Hisat2/chlamydomonas_reinhardtii/Chlamy_V55
    cmd = f' hisat2 -p 30 --rna-strandness RF -x {builtGenome} -U {folder}/{file} -S {folder}/{output} --dta 2> {folder}/{mapping_stats}'
    print(f'running: {cmd}')
    os.system(cmd)

#troubleshooting stop program
def stopProgram():
    repeat = True
    while repeat == True:
        contFunc = input("Continue running program? (y/n):\n")
        if contFunc.casefold() == "n":
            sys.exit("\nExiting...\n")
        elif contFunc.casefold() == "y":
            print("Continuing...\n")
            
def findPeaks(csRNA, csInput, totalRNA, ntagsCount):
    #run the findcsRNATSS.pl
    print(f'Peak finding with samples csRNA: {csRNA} csInput: {csInput} totalRNA: {totalRNA}\n')
    csDir = f'{path}/tagDirs/csRNA'
    inputDir = f'{path}/tagDirs/csInput'
    os.system(f'findcsRNATSS.pl {csDir}/{csRNA} -o ./analysis/peakFiles/{csRNA}.peaks -i {inputDir}/{csInput} -genome hg38 -size 150 -ntagThreshold {ntagsCount}')
    
def starFunc(starDir,file):
    starCmd = f'STAR --genomeDir {starDir} --runThreadN 8 --readFilesIn {file} --outFileNamePrefix {file}. --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterMultimapNmax 10000 --limitOutSAMoneReadBytes 10000000 '
    os.system(starCmd)

#basic zipping function to allow multiprocessing of the gzip command    
def zipFunc(zipFile):
    cmd = f'gzip (file)'
    os.system(cmd)

## information parser for running
parser=argparse.ArgumentParser(description = 'Auttomated QC function. this code will take files in either the current directory, or designated directory, and analyze all fastq.gz for QC analysis')

# define user inputs to program through linux command line
parser.add_argument("-p",
                    "--species", 
                    help="species you are working with. If left blank, generic naming schemes will be used instead")
parser.add_argument("-i",
                    "--hisat2Align", 
                    help="use the Hisat2 aligner for analysis. Must include this or STAR for program to run",
                    action="store_true", 
                    default=False)
parser.add_argument('-o',
                    '--homerGenome',
                    help="use specified Homer genome for tagDir creation",
                    required=True)
parser.add_argument("-s",
                    "--starAlign", 
                    help="use the STAR aligner for analysis. Must include this or Hisat2 for program to run",
                    action="store_true", 
                    default=False)
parser.add_argument("-r",
                    "--raw", 
                    help="location of the folder with raw fastq.gz files. If blank, program will use pwd instead")
parser.add_argument("-g",
                    "--genomeDir", 
                    required=True, 
                    help="path to the genome diretory used for the given species. program will try to detect hisat prefix if Hisat used")
parser.add_argument("-b",
                    "--rebuildTagDir",
                    help= "Optional rebuilding the tagdirs according to the scripts needs. Script will attempt to autodetect tagdir presence",
                    action ="store_true",
                    default=False)
args=parser.parse_args()
    
# pull from parser
useHisat = args.hisat2Align
useStar = args.starAlign
rawPath = args.raw
genomeDir = args.genomeDir
rebuildTD = args.rebuildTagDir
homerGenome = args.homerGenome

# check to see if either Hisat or Star were used. if neither, the program will stop
if (useHisat == False) and (useStar == False):
    exit('\n\nNeither Hisat or Star was specified for use.\nPlease choose one, the other, or both\nExiting...') 

#define the genome directories
if useHisat == True:
    if genomeDir.endswith('/'):
        hisatGenomeDir = genomeDir+'index'
    else:
        hisatGenomeDir = genomeDir+'/index'

if useStar == True:
    starGenomeDir = genomeDir

#no rawPath was declared, so use cwd
if rawPath == None:
    rawPath = os.getcwd()
if rawPath == './':
    rawPath == os.getcwd()

# begin defining some variables to use for the rest of the program, i.e. parent directories, species to use genome for, etc.
lastDir = rawPath.split('/')[-1]
if len(lastDir) == 0:
    lastDir = rawPath.split('/')[-2]
parentDir = rawPath.split(('/'+lastDir))[0]

# also take the time now to build up the other directories needed
# How it needs to look:
#                
#                         ParentDir________
#                        /    |   \        |
#               ____analysis  QC  tagDirs  raw
#              /     |     \         |    /   \
#       keyFiles peakFiles diffExp   | csRNA  csInput                         
#                                   / \                
#                               csRNA  csInput
#
# loading |\--/|\--/|

maybe = True
while maybe == True:

    ## try moving to the main parent dir file. this should be the raw file -1
    if os.path.exists(parentDir):
        os.chdir(parentDir)
    else:
        print(f'{parentDir} does not exist')
        sys.exit()
    
    ## check the rest of the dir based user defined variables
    HisatPrefix = genomeDir.split('/')[-1]
    justDir = genomeDir.split(HisatPrefix)[0]
    if not os.path.exists(justDir):
        print(f'{justDir} does not exist')
        sys.exit()

    mainDirs = ['analysis','QC','tagDirs','raw']
    analysisDirs = ['keyFiles','peakFiles','diffExp']
    triDirs = ['csRNA','csInput','totalRNA','other']

    print('creating needed directories if not present')
    for eachDir in mainDirs:
        os.chdir(parentDir)
        if not os.path.exists(f'{parentDir}/{eachDir}'):
            os.system(f'mkdir {eachDir}')

        if eachDir == 'analysis':
            os.chdir(f'{parentDir}/{eachDir}')
            for subDir in analysisDirs:
                if not os.path.exists(f'{parentDir}/{eachDir}/{subDir}'):
                    os.system(f'mkdir {subDir}')
        elif (eachDir == 'tagDirs') or (eachDir == 'raw'):
            os.chdir(f'{parentDir}/{eachDir}')
            for subDir in triDirs:
                if not os.path.exists(f'{parentDir}/{eachDir}/{subDir}'):
                    os.system(f'mkdir {subDir}')

    # now start moving files in the raw folder 
    print('\tZipping fastq files to save space...')
    processes = []
    for eachFile in os.listdir(rawPath):
        #do some zipping first
        if eachFile[-len('.fastq'):] == '.fastq':
            print(f'\t\tZipping {eachFile}...')
            zipPro = multiprocessing.Process(target = zipFunc, args = (eachFile,), name = "zipPro_"+eachFile.split('.')[0])
            zipPro.start()
            #add to list to join later
            processes.append(zipPro)

            # pause to let processes clear
            if len(processes) >= 4:
                print('\tWaiting for processes to clear')
                multiprocessing.connection.wait(p.sentinel for p in processes)
        else:
            continue

    #wait for the rest to finish
    if len(processes) != 0:
        print('\twaiting for zipping to finish')
        multiprocessing.connection.wait(p.sentinel for p in processes)

    print('Adjusting file locations')
    for eachFile in os.listdir(rawPath):
        if os.path.isdir(f'{rawPath}/{eachFile}'):
            continue
        elif ('csRNA' in eachFile) or ('GRO' in eachFile):
            cmd = f'mv {rawPath}/{eachFile} {rawPath}/csRNA/'
            os.system(cmd)
        elif ('csInp' in eachFile) or ('_sRNA' in eachFile):
            cmd = f'mv {rawPath}/{eachFile} {rawPath}/csInput/'
            os.system(cmd)
        elif ('totalRNA' in eachFile) or ('_RNAseq' in eachFile):
            cmd = f'mv {rawPath}/{eachFile} {rawPath}/totalRNA/'
            os.system(cmd)
        else:
            cmd = f'mv {rawPath}/{eachFile} {rawPath}/other/'
            os.system(cmd)

    ## now with files moved to where we want them, 
    ## look to see if tagdirs already exists. if so, skip trimming, alinging, and tagdir creation
    
    # set my variable for tagdir creation or not
    buildTagDir = True

    #start looping all the raw files
    dirCount = 0
    rawCount = 0
    if os.path.exists(f'{parentDir}/tagDirs/'):
        for eachRaw in os.listdir(f'{parentDir}/raw/'):
            if os.path.isdir(f'{parentDir}/raw/{eachRaw}'):
                os.chdir(f'{parentDir}/raw/{eachRaw}')

                #loop all the raw files and check if has a tagdir:
                for eachFile in os.listdir():
                    if (eachFile[-len('.fastq.gz'):] == '.fastq.gz') or (eachFile[-len('.fastq'):] == '.fastq'):
                        rawCount += 1
                        #convert the raw file name into what the tagDir would be
                        if ('csRNA' in eachFile) or ('GRO' in eachFile):
                            ## create a file path for this tagdir to be from this file name
                            if '_L0' in eachFile:
                                modVal = eachFile.split('_L0')[0]
                            else:
                                modVal = eachFile.split('.')[0]
                            
                            dirName = f'{parentDir}/tagDirs/csRNA/{modVal}'
                            #check if that directory exists
                            if os.path.exists(dirName):
                                dirCount += 1

                        elif ('csInp' in eachFile) or ('_sRNA' in eachFile):
                            ## same here
                            if '_L0' in eachFile:
                                modVal = eachFile.split('_L0')[0]
                            else:
                                modVal = eachFile.split('.')[0]

                            dirName = f'{parentDir}/tagDirs/csInput/{modVal}'
                            if os.path.exists(dirName):
                                dirCount += 1

                        elif ('totalRNA' in eachFile) or ('_RNAseq' in eachFile):
                            ## same here
                            if '_L0' in eachFile:
                                modVal = eachFile.split('_L0')[0]
                            else:
                                modVal = eachFile.split('.')[0]

                            dirName = f'{parentDir}/tagDirs/totalRNA/{modVal}'
                            if os.path.exists(dirName):
                                dirCount += 1


    #print(f'\n\n\ndircount: {dirCount}; rawcount: {rawCount}\n\n\n')
    if dirCount == rawCount:
        buildTagDir = False
        
    ## check user submitted requirements of rebuilding tagdirs
    if rebuildTD == True:
        buildTagDir = True 
    print(f'build tag directories: {buildTagDir}')    
    if buildTagDir == True: 
        ## begin trimming and alinging file in each subfolder
        for eachDir in os.listdir(rawPath):
            if os.path.isdir(f'{rawPath}/{eachDir}'):
                os.chdir(f'{rawPath}/{eachDir}')
                processes = []
                print(f'Begining trimming in {eachDir}...')

                if len(os.listdir()) == 0:
                    #skip empty directories
                    continue

                for eachFile in os.listdir():
                    if eachFile[-len('.fastq.gz'):] == '.fastq.gz' or eachFile[-len('.fastq'):] == '.fastq':
                        # multiprocess the file
                        print(f'Began trimming {eachFile}...')
                        trimPro = multiprocessing.Process(target = trimFunc, args = (eachFile,), name = f'trimPro_{eachFile}')
                        trimPro.start()
                        processes.append(trimPro)

                #check if length of processes is too long and wait for active ones to finish
                    if len(processes) >= 10:
                        print('\nWaiting for processes to clear')
                        timerStart = time.perf_counter()
                        multiprocessing.connection.wait(p.sentinel for p in processes)
                        timerFinish = time.perf_counter()
                        print(f"elapsed time: {round(timerFinish-timerStart, 2)} seconds(s)")

                # waiting for all process to finish
                if len(processes) != 0:
                    multiprocessing.connection.wait(p.sentinel for p in processes)
                print(f'Finished trimming {eachDir}.')
            else:
                continue

            ## determine wether or not to use the Hisat2 aligner or the STAR aligner:
            # align some files with Hisat2!
            if useHisat == True and useStar == True:
                print('Using both aligners... waiting 10 seconds to quit if other option wanted')
                for i in range(10):
                    time.sleep(1)
                    print(i)                

            if useHisat == True:
                processes = []
                print(f'Begining aligning in {eachDir}')
                for file in os.listdir():
                    if file[-8:] == '.trimmed':
                        ## run the subprocess with the function above. p.s. don't run more than 10 at once
                        print(f'Began Hisat2 alinging {file}...')
                        filePath =f'{parentDir}/raw/{eachDir}'
                        hisatPro = multiprocessing.Process(target = hisat2Func, args = (file,genomeDir,filePath,), name = "hisatPro_"+file.split('.')[0])
                        hisatPro.start()
                        #add to list to join later
                        processes.append(hisatPro)

                    # pause to let processes clear
                    if len(processes) >= 10:
                        print('\nWaiting for processes to clear')
                        timerStart = time.perf_counter()
                        multiprocessing.connection.wait(p.sentinel for p in processes)
                        timerFinish = time.perf_counter()
                        print(f"elapsed time: {round(timerFinish-timerStart, 2)} seconds(s)")
                #again, wait for all processes to clear
                multiprocessing.connection.wait(p.sentinel for p in processes)

            # run the star Alinger
            if useStar == True:
                processes = []
                print(f'Begining aligning in {eachDir}')
                for file in os.listdir():
                    if file[-8:] == '.trimmed':
                        ## run the subprocess with the function above. p.s. don't run more than 10 at once
                        print(f'Began STAR alinging {file}...')
                        filePath =f'{parentDir}/raw/{eachDir}'
                        starPro = multiprocessing.Process(target = starFunc, args = (genomeDir,file,), name = "starPro_"+file.split('.')[0])
                        starPro.start()
                        #add to list to join later
                        processes.append(starPro)

                    # pause to let processes clear
                    if len(processes) >= 10:
                        print('\nWaiting for processes to clear')
                        timerStart = time.perf_counter()
                        multiprocessing.connection.wait(p.sentinel for p in processes)
                        timerFinish = time.perf_counter()
                        print(f"elapsed time: {round(timerFinish-timerStart, 2)} seconds(s)")
                #again, wait for all processes to clear
                multiprocessing.connection.wait(p.sentinel for p in processes)


        ## files built and bidoodled, aka aligned, time to make some tag directories
        # first, return to parent dir
        os.chdir(parentDir)
        pwd = parentDir
        samLst = []
        for folder in os.listdir(f'{parentDir}/raw/'):
            os.chdir(f'{parentDir}/raw/{folder}')

            if len(os.listdir()) == 0:
                #skip the empty directories
                continue

            for eachFile in os.listdir():
                #print(eachFile[-len('.sam'):])
                if eachFile[-len('.sam'):] == '.sam':
                    ## append to a list for writing to a file
                    samLst.append(eachFile)
        os.chdir(pwd)

        #create a file with just the sam names in it
        with open(f'{parentDir}/QC/samNames.txt','w') as samFiles:
            for eachSam in samLst:
                samFiles.write(f'{eachSam}\n')

        ## write some naming scheme and append to a new file with the QC samNames and such
        with open(f'{parentDir}/QC/sampleInfoFile.txt','w') as namesFile:
            for eachVal in samLst:
                if ('csRNA' in eachVal) or ('GRO' in eachVal):
                    ## create a file path for this tagdir to be from this file name
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]
                    dirName = f'{parentDir}/tagDirs/csRNA/{modVal}'
                    line = f'{dirName}\traw/csRNA/{eachVal}'
                    namesFile.write(f'{line}\n')

                elif ('csInp' in eachVal) or ('_sRNA' in eachVal):
                    ## same here
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]

                    dirName = f'{parentDir}/tagDirs/csInput/{modVal}'
                    line = f'{dirName}\traw/csInput/{eachVal}'
                    namesFile.write(f'{line}\n')

                elif ('totalRNA' in eachVal) or ('_RNAseq' in eachVal):
                    ## same here
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]

                    dirName = f'{parentDir}/tagDirs/totalRNA/{modVal}'
                    line = f'{dirName}\traw/totalRNA/{eachVal}'
                    namesFile.write(f'{line}\n')


        ## now with the path and file locations written, lets use the batchmaketagdirs command  to make all tagdirs simultaneously
        if 'Hisat2' in genomeDir:
            tagDirGenome = genomeDir.split('/Hisat2')[0] + 'genome.fa'
            
        cmdBatch = f'batchMakeTagDirectory.pl QC/sampleInfoFile.txt -cpu 4 -genome {homerGenome}  -omitSN -checkGC -fragLength 150 -single'
        os.system(cmdBatch)
    else:
        #create the sample info files
        os.chdir(parentDir)
        pwd = parentDir
        samLst = []
        for folder in os.listdir(f'{parentDir}/raw/'):
            os.chdir(f'{parentDir}/raw/{folder}')

            if len(os.listdir()) == 0:
                #skip the empty directories
                continue

            for eachFile in os.listdir():
                #print(eachFile[-len('.sam'):])
                if eachFile[-len('.sam'):] == '.sam':
                    ## append to a list for writing to a file
                    samLst.append(eachFile)
        os.chdir(pwd)

        #create a file with just the sam names in it
        with open(f'{parentDir}/QC/samNames.txt','w') as samFiles:
            for eachSam in samLst:
                samFiles.write(f'{eachSam}\n')

        ## write some naming scheme and append to a new file with the QC samNames and such
        with open(f'{parentDir}/QC/sampleInfoFile.txt','w') as namesFile:
            for eachVal in samLst:
                if ('csRNA' in eachVal) or ('GRO' in eachVal):
                    ## create a file path for this tagdir to be from this file name
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]

                    dirName = f'{parentDir}/tagDirs/csRNA/{modVal}'
                    line = f'{dirName}\traw/csRNA/{eachVal}'
                    namesFile.write(f'{line}\n')

                if ('csInput' in eachVal) or ('_sRNA' in eachVal):
                    ## same here
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]

                    dirName = f'{parentDir}/tagDirs/csInput/{modVal}'
                    line = f'{dirName}\traw/csInput/{eachVal}'
                    namesFile.write(f'{line}\n')

                if ('totalRNA' in eachVal) or ('_RNAseq' in eachVal):
                    ## same here
                    if '_L0' in eachVal:
                        modVal = eachVal.split('_L0')[0]
                    else:
                        modVal =eachVal.split('.')[0]

                    dirName = f'{parentDir}/tagDirs/totalRNA/{modVal}'
                    line = f'{dirName}\traw/totalRNA/{eachVal}'
                    namesFile.write(f'{line}\n')

    ## this section is what starts the actual QC graphs, mapping stats, pcr clonality, etc, ad will save these files 
    ## with appropriate names to the QC file.

    #############################################################################################################################################################################
    ################################################################ Mapping Stats ##############################################################################################
    #############################################################################################################################################################################
    os.chdir(parentDir)

    #start by creating copies of the .sam files where the only difference is that they have .Log.final.out at the end as opposed to trimmed.sam
    os.chdir(f'{parentDir}/raw/')
    for eachDir in os.listdir(f'{parentDir}/raw/'):
        if len(os.listdir(f'{parentDir}/raw/{eachDir}')) == 0:
            continue
        else:
            os.chdir(f'{parentDir}/raw/{eachDir}')
            for eachFile in os.listdir():
                if eachFile[-len('trimmed.sam'):] == 'trimmed.sam':
                    name = eachFile.split('.sam')[0]
                    newName = name+'.Log.final.out'
                    copyCmd = f'cp {eachFile} {newName}'
                else:
                    continue
    os.chdir(parentDir)
    
    ## getHomerMappingStats.pl doesn't completely work with hisat mapping stat items, so i need to manualy pull them if the user 
    ## wants to use Hisat alinging
    # this will be done after the hisat alinger is finished
    statsLst = []
    #print(f"Didn't Align\tUnique Reads\tAligned Multiple\tOverall Alignment Rate\tCount Didn't Align\tCount Unique Reads\tCount Aligned Multiple\tCount Overall Alignment Rate")
    # Write this information to the file below
    with open(f'{parentDir}/QC/stats.txt', 'w') as statstxt:
        statstxt.write("Samp Name\tDidn't Align\tUnique Reads\tAligned Multiple\tOverall Alignment Rate\tCount Didn't Align\tCount Unique Reads\tCount Aligned Multiple\tCount Overall Alignment Rate\n")
        for eachDir in os.listdir(rawPath):
            if os.path.isdir(f'{rawPath}/{eachDir}'):
                os.chdir(f'{rawPath}/{eachDir}')
                for eachFile in sorted(os.listdir()):
                    if eachFile[-len('mappingstats.txt'):] == 'mappingstats.txt':

                        # append this to a list for later use if need be
                        statsLst.append(eachFile)

                        #load open the file on a per line basis
                        mapStatsFile = open(eachFile)

                        #get the name of the file
                        name = re.split('_S[0-9][0-9]_',eachFile)[0]

                        ##start by creating a string seperated by tabular deliminted values
                        sumString = ""
                        otherString = ""
                        for eachLine in mapStatsFile:
                            if 'Warning' in eachLine:
                                continue

                            elif 'overall' in eachLine:
                                percentValue = eachLine.split('%')[0]
                                percentValue = percentValue+'%'
                                sumString += percentValue+"\t"
                                #print(f'overall;  whole:{wholeValue} eachLine: {eachLine}')

                            elif 'reads; of these' in eachLine:
                                wholeValue = eachLine.split(' ')[0]
                                otherString += wholeValue+'\t'
                                #print(f'of these;  whole:{wholeValue} eachLine: {eachLine}')

                            else:
                                percentValue = eachLine.split(' (')[-1].split(')')[0]
                                sumString += percentValue+"\t"
                                wholeValue = eachLine.split(' (')[0]
                                otherString += wholeValue+'\t'
                                #print(f'else;  whole:{wholeValue} eachLine: {eachLine}')

                        ## combine the two strings
                        compleat = sumString + otherString
                        statstxt.write(name+'\t'+compleat+'\n')
                        #print(name+'\t'+compleat)
                        mapStatsFile.close()                        
    #############################################################################################################################################################################
    ################################################################ UCSD Files for viewing #####################################################################################
    #############################################################################################################################################################################

    ## loop through all of the tagDir files for the csRNA ones and create ucsc files. Only need the combo files
    os.chdir(f'{parentDir}/tagDirs/')
    pwd = os.getcwd()

    for eachDir in os.listdir():
        if eachDir == 'csRNA':
            # jump into the file
            os.chdir(f'{pwd}/{eachDir}')

            # loop through the files again
            for eachTagDir in  os.listdir():
                #make a ucsc file of the dir
                cmd = f'makeUCSCfile {eachTagDir} -strand separate -style tss > {eachTagDir}/{eachTagDir}.ucsc.bedGraph'
                os.system(cmd)

        # don't do anything fot the csInp or totalRNA tagdirs if available
        else:
            continue

    #############################################################################################################################################################################
    ##################################################################### PCR and Clonality #####################################################################################
    #############################################################################################################################################################################

    output = f'{parentDir}/QC'
    os.chdir(f'{parentDir}/tagDirs')
    pwd = os.getcwd()
    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            tagdirs.append(eachLine)

    patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
    patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Val'})
    patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]

    my_dict = {"Library":[],"Median tags per tag position (should be =1)":[]}

    for f in tagdirs:
        file = f +'/tagCountDistribution.txt'
        #read in  file
        read_tagCountDist_file = pd.read_csv(file, sep='\t')

        name = list(read_tagCountDist_file.columns)
        median_val = str(name).split('=')[1].split(',')[0]

        plotName = f.split('/')[-1]
        my_dict["Library"].append(plotName)
        my_dict["Median tags per tag position (should be =1)"].append(median_val)

    Median_frame = pd.DataFrame(my_dict)
    Median_frame["Median tags per tag position (should be =1)"] = pd.to_numeric(Median_frame["Median tags per tag position (should be =1)"]) #convert column to numeric values for plotting

    graph = sns.barplot(data=Median_frame,  y = 'Library', x="Median tags per tag position (should be =1)", capsize=.4, errcolor=".5", linewidth=2, edgecolor=".5", facecolor=(0, 0, 0, 0),)
    graph.axvline(1, color='g', ls='--')
    graph.axvline(1.2, color='r')

    sns.set(rc={'figure.figsize':(6,6)})
    plt.savefig(f'{output}/medianTagsPerPosition.png')


    path = f'{parentDir}/tagDirs/'
    output_dir = f'{parentDir}/QC/'
    os.chdir(f'{parentDir}/raw')

    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]

            eachLine = '/'+eachLine
            tagdirs.append(eachLine)

    patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
    patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Tags per tag position'})
    patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]

    for f in tagdirs:
        file = f +'/tagCountDistribution.txt'
        #read in  file
        read_tagCountDist_file = pd.read_csv(file, sep='\t')
        #define variables for names
        median_value = str(list(read_tagCountDist_file)).split('=')[1].split(',')[0].split(' ')[1] #gets the median value out of the columnname
        modName = f.split('/')[-1]
        column_name = modName + ' (' + median_value + ')'

        #rename columns for merge
        read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[0]: 'Tags per tag position'})
        read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[1]: column_name})
        #cat all tagDir values together
        merged_frame = pd.merge(patens_tagCountDist_start, read_tagCountDist_file, left_on='Tags per tag position', right_on='Tags per tag position', how='left')
        patens_tagCountDist_start = merged_frame.copy()

    patens_tagCountDist_start = patens_tagCountDist_start.set_index('Tags per tag position')
    patens_tagCountDist_start.to_csv(output_dir + 'sum_TagCountDist.txt', sep = '\t')

    #convert into a df for funky scatter plots
    combined_frames = patens_tagCountDist_start.replace(['0', 0], np.nan)  #make 0 to NaN 
    combined_frames_stacked = combined_frames.stack().reset_index() #stack and reset index
    combined_frames_stacked.columns = ['Tags per tag position','Library','Fraction of Positions'] #rename columns
    combined_frames_stacked_log = combined_frames_stacked
    combined_frames_stacked_log['Tags per tag position'] = np.log(combined_frames_stacked_log['Tags per tag position']) #log 
    combined_frames_stacked_log['Fraction of Positions'] = np.log(combined_frames_stacked_log['Fraction of Positions']) #log 

    g = sns.FacetGrid(combined_frames_stacked_log, col='Library',height=8, aspect=1)
    plots = g.map(sns.scatterplot, "Tags per tag position", "Fraction of Positions", alpha=.5)
    plt.savefig(f'{output}/tagsPer_Vs_FracofPos.png')

    #############################################################################################################################################################################
    ##################################################################### Read Length Dist. #####################################################################################
    #############################################################################################################################################################################
    plt.close()
    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            tagdirs.append(eachLine)

    tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
    tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
    tagLength_start = tagLength_start.iloc[:,[0]]

    my_dict = {"Library":[],"Average tag length":[]}

    for f in tagdirs:
        file = f +'/tagLengthDistribution.txt'
        #read in  file
        read_tagLength_file = pd.read_csv(file, sep='\t')

        name = list(read_tagLength_file.columns)
        tagLength_val = str(name).split('= ')[1].split(')')[0]

        modName = f.split('/')[-1]
        my_dict["Library"].append(modName)
        my_dict["Average tag length"].append(tagLength_val)

    AvrgLength_frame = pd.DataFrame(my_dict)
    AvrgLength_frame["Average tag length"] = pd.to_numeric(AvrgLength_frame["Average tag length"]) #convert column to numeric values for plotting

    graph = sns.barplot(data=AvrgLength_frame,  y = 'Library', x="Average tag length", capsize=.4, errcolor=".5", linewidth=2, edgecolor=".5", facecolor=(0, 0, 0, 0),)
    #graph.axvline(1, color='g', ls='--')
    #graph.axvline(1.2, color='r')
    sns.set(rc={'figure.figsize':(6,6)})
    plt.tight_layout()
    plt.savefig(f'{output}/read_Length_Dist.png')

    #############################################################################################################################################################################
    ##################################################################### RNA Length plots ######################################################################################
    #############################################################################################################################################################################
    plt.close()
    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            tagdirs.append(eachLine)

    tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
    tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
    tagLength_start = tagLength_start.iloc[:,[0]]

    for f in tagdirs:
        file = f +'/tagLengthDistribution.txt'
        #read in  file
        read_tagLength_file = pd.read_csv(file, sep='\t')
        #define variables for names
        column_name = f.split('/')[-1] 

        #rename columns for merge
        read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[0]: 'Length (nt)'})
        read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[1]: column_name})
        #cat all tagDir values together
        merged_frame = pd.merge(tagLength_start, read_tagLength_file, left_on='Length (nt)', right_on='Length (nt)', how='left')
        tagLength_start = merged_frame.copy()


    tagLength_start = tagLength_start.iloc[1: , :] #delete 0 cause not a read after trimming
    tagLength_start = tagLength_start.set_index('Length (nt)')
    tagLength_start.to_csv(output_dir + 'sum_TagLengthDist.txt', sep = '\t')

    combined_frames = tagLength_start.replace(['0', 0], np.nan)  #make 0 to NaN 
    combined_frames_stacked = combined_frames.stack().reset_index() #stack and reset index
    combined_frames_stacked.columns = ['Length (nt)','Library','Fraction of Reads'] #rename columns

    g = sns.FacetGrid(combined_frames_stacked, col='Library', height=8, aspect=1)
    g.map(sns.lineplot, 'Length (nt)','Fraction of Reads', alpha=.5)
    plt.tight_layout()
    plt.savefig(f'{output}/RNA_Length_Plots.png')

    #############################################################################################################################################################################
    ##################################################################### NT Dist. Plots ########################################################################################
    #############################################################################################################################################################################
    plt.close()
    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            tagdirs.append(eachLine)

    #complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
    nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
    nt_freq_file_start = nt_freq_file_start.iloc[:,[0]]

    for f in tagdirs:
        nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
        plot_me = nt_freq_file.set_index('Offset')
        plot_me_stack = plot_me.stack().reset_index() #stack and reset index
        axisName = f.split('/')[-1]
        plot_me_stack.columns = ['Distance from TSS','nt',axisName + ' - %'] #rename columns

        sns.set(rc={'figure.figsize':(20,15)})
        plt.figure()
        g = sns.lineplot(data=plot_me_stack, x="Distance from TSS", y= axisName + ' - %', hue="nt")
        plt.tight_layout()
        plt.savefig(f'{output}/{axisName}_nt_Preference.jpeg')


    #############################################################################################################################################################################
    ################################################################### percent nt A Plots ######################################################################################
    #############################################################################################################################################################################
    plt.close()
    tagdirs = []
    with open('../QC/sampleInfoFile.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            tagdirs.append(eachLine)

    #select key parameter
    nucleotide_to_plot = 'A'
    selected_files = 'csRNA'

    #complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
    nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
    nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
    nt_freq_file_start = nt_freq_file_start.set_index('Offset')
    nt_freq_file_start = nt_freq_file_start.stack().reset_index() #stack and reset index
    nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
    nt_freq_file_start = nt_freq_file_start.iloc[: , :2] #select both offset and nt info

    unused = []
    for f in tagdirs:
        if selected_files in f:
            nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
            nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
            plot_me = nt_freq_file.set_index('Offset')
            plot_me_stack = plot_me.stack().reset_index() #stack and reset index
            modName = f.split('/')[-1]
            plot_me_stack.columns = ['Distance from TSS','nt', modName] #rename columns

            del plot_me_stack['Distance from TSS']
            del plot_me_stack['nt']

            merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
            nt_freq_file_start = merged_frame.copy()
        else:
            #append unused values to a list
            unused.append(f)



    AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
    del AntPlot_frame['nt']
    AntPlot_frame = AntPlot_frame.set_index('Offset')
    AntPlot_frame_stack = AntPlot_frame.stack().reset_index() #stack and reset index
    AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] #rename columns

    #plot with legend next to plot
    plt.figure()
    ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(f'{output}/csRNA_Aplot.png')

    ################################################################### redo with those libs not plotted in first plot ###################################################################

    if not len(unused) == 0:
        #complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
        nt_freq_file_start = nt_freq_file_start.set_index('Offset')
        nt_freq_file_start = nt_freq_file_start.stack().reset_index() #stack and reset index
        nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
        nt_freq_file_start = nt_freq_file_start.iloc[: , :2] #select both offset and nt info

        for f in unused:
            nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
            nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
            plot_me = nt_freq_file.set_index('Offset')
            plot_me_stack = plot_me.stack().reset_index() #stack and reset index
            modName = f.split('/')[-1]
            plot_me_stack.columns = ['Distance from TSS','nt', modName] #rename columns

            del plot_me_stack['Distance from TSS']
            del plot_me_stack['nt']

            merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
            nt_freq_file_start = merged_frame.copy()

        AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
        del AntPlot_frame['nt']
        AntPlot_frame = AntPlot_frame.set_index('Offset')
        AntPlot_frame_stack = AntPlot_frame.stack().reset_index() #stack and reset index
        AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] #rename columns

        #plot with legend next to plot
        plt.figure()
        ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(f'{output}/csInput_Aplot.png')    
    
    ## stop the program
    maybe = False
