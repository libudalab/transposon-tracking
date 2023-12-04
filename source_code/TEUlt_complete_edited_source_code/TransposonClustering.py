from __future__ import print_function
############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import os
import shutil

# Methods
def getQuery(records, queryChrom, start, end, seq):
    return SeqFeature(FeatureLocation(start, end), type="CDS").extract(seq)

def findOptimal_FactorForClusteringBatchSize(start, end, total): # start = batchSizeStart, end = batchSizeEnd, total = counter returned line 99 (total number of transposons found from gff3 candidates file)
    limit = start*1.0
    lastND = 1  # what is ND?
    lastNU = total # what is NU?
    lastTimes = 0
    n = lastND+(lastNU-lastND)/2 # n = (ND + (NU-ND)) /2
    while True:
        F = (end/start)**(1.0/n) # start = batchSizeStart, end = batchSizeEnd
        times = 0
        target = 0
        counter = 0
        last = 1
        counter = 1
        for i in range(0,total): # for i in number of total transposons
            counter += 1 # add 1 to counter
            if(counter>limit): # if the counter value is larger than the limit
                times += 1 # add 1 to times
        #        print(i-last,"\t",i,"\t",limit)
                target = i-last # target = transposon index value minus the 'last' value
                last = i # set 'last' equal to the transposon index number
                counter = 0 # reset counter to 0
                limit = limit * F # set limit equal to the limit times F
                if(limit < end): # if the limit value is less than the end value
                    limit = end*1.0 # set limit equal to end value time 1
    #    print(target)
        if(target-2>end):
            lastNU = n
        else:
            if(lastTimes<times):
                lastND = n
                lastTimes = times
            else:
                break
        n = int(lastND+(lastNU-lastND)/2)
#        print(n,"\t",lastND,"\t",lastNU,"\t",target,"\t",times)
        if(lastNU-lastND<=2):
            break
    return n, F

def createTransposonSequences(folderIn, folderOut, fastaFile):
    MAXLEN=100000
    # Determine latest iteration
    files = os.listdir(folderIn)
    filesN = list()
    for f in files:
        if(f.startswith("CandidatesB_it")):
            filesN.append(int(f.replace("CandidatesB_it","").replace(".gff3","")))
    filesN.sort()
    fileIn = os.path.join(folderIn,"CandidatesB_it"+str(filesN[-1])+".gff3")
    # Load Annotations, sort by size
    transposons = {}
    counter = 0
    f = open(fileIn,"r")
    line = f.readline()
    while line!="": #if not a blank line
        if(not line.startswith("#")): # if not a header line
            parts = line.split("\t") # parts = line split at tabs
            if(parts[2]=="transposon"): # if transposon
                if(parts[0] not in transposons): # if chromosome isn't in transposons dict keys
                    transposons[parts[0]] = list() # add chromosome to transposon dict as a key with an empty list as value
                transposons[parts[0]].append([parts[0],int(parts[3]),int(parts[4]),int(parts[8]),int(parts[4])-int(parts[3])]) # append transposon info to chromosome key's list of tes in transposons dictt
                counter+=1 # add one to counter
        line = f.readline() # next line
    f.close()
    # Create Fasta File of sequences
    fileOut = os.path.join(folderOut,"TransposonSequences.fasta") # initialize file
    f = open(fileOut,"w+")
    records = SeqIO.parse(fastaFile, "fasta") # use seqio to parse input fasta file (fastaFile = "parsedAnnotations/sequence.fasta")
    for r in records: # for each seq record from input fasta
        for key in list(transposons.keys()): # for each chromosome key in the transposons dict
            if(r.id==key): # if the seq record id is the chromosome key from the transposons dict
                seq = r.seq # assign variable seq to the sequence of that te from the input fasta
                for idx in range(0,len(transposons[key])): # for each te
                    query = transposons[key][idx]
                    queryChrom = query[0]
                    queryStart = query[1]
                    queryEnd   = query[2]
                    queryID    = query[3]
                    if(abs(queryEnd-queryStart)<MAXLEN):
                        querySequence   = getQuery(r, queryChrom, queryStart, queryEnd, seq) # get sequence of TE from the input fasta
                        f.write(">transposon"+str(queryID)+"\n"+str(querySequence)+"\n") # write the TE id and TE sequence to the TransposonSequences.fasta
    f.close()
    return counter # return number of transposons that were found and put into transposons dict value lists

# Parameters
def doTransposonClustering(folderIn,folderOut,fastaFile,batchSizeStart=200,batchSizeEnd=10):
#    folderIn  = "transposonCandB"
#    folderOut = "transposonCandC"
#    fastaFile = "parsedAnnotations/sequence.fasta"
#    batchSizeStart = 200
#    batchSizeEnd   = 10
#    counter = 21498

    # Create Transposon Sequences
    counter = createTransposonSequences(folderIn, folderOut, fastaFile) # line 56
    print("Transposon Sequences extracted...")

    # Sort Transposon Sequences by size
    os.system("seqkit sort -l "+os.path.join(folderOut,"TransposonSequences.fasta")+" > "+os.path.join(folderOut,"TransposonSequences_sort.fasta"))
    print("Transposon Sequences sorted...")

    # Create folder with parts of sequences for clustering
    if(os.path.isdir(os.path.join(folderOut,"clusterIt1"))):
        shutil.rmtree(os.path.join(folderOut,"clusterIt1"))
    os.mkdir(os.path.join(folderOut,"clusterIt1"))

    # Start to cluster small portions
    n,F = findOptimal_FactorForClusteringBatchSize(batchSizeStart,batchSizeEnd,counter) # line 17, n = , F =
    batchSize = batchSizeStart
    counterZ = 0
    fileCounter = 1
    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences_sort.fasta"), "fasta") # parse seq records from transposon sequences fasta sorted by size
    selection = list() # initialize selection list
    for r in records: # for every transposon sequence
        if(counterZ>batchSize): # if the counter is larger than the batchsize
            print(counterZ, batchSize) # print the counter and batchsize
            SeqIO.write(selection, os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta"), "fasta") # make a new candidates b file
            os.system("cd-hit -i "+os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta")+" -o "+os.path.join(folderOut,"clusterIt1","Results"+str(fileCounter)+".txt")) # run cd hit on new file
            counterZ = 0
            selection = list()
            fileCounter += 1
            batchSize = batchSize*F
        else:
            counterZ += 1
            selection.append(r)
    print("Clustering finished successfully...")
    SeqIO.write(selection, os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta"), "fasta")
    os.system("cd-hit -i "+os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta")+" -o "+os.path.join(folderOut,"clusterIt1","Results"+str(fileCounter)+".txt"))
    fileCounter += 1
    counterZ = 0
    selection = list()

    # Collect cluster results
    clusterData = {}
    clusterFileAss = {}
    collected = list()
    previousClusterNo = -1
    previousClusterB = 0
    for i in range(1,fileCounter): # for range i to number of files
        previousClusterNo = previousClusterB # previousClusterNo equals previousClusterB
        f = open(os.path.join(folderOut,"clusterIt1","Results"+str(i)+".txt.clstr"),"r") #open cluster results file i
        line = f.readline()
        while line!="": # while not blank line
            if(line.startswith(">")): # if line with cluster number ('>Cluster #')
                if(len(collected)==0): # if the collected list has no entries pass the line
                    pass
                else:
                    clusterData[previousClusterNo] = collected # set the key value in clusterData of previousClusterNo to the value of collected
                    clusterFileAss[previousClusterNo] = i # set the clusterFileAss key of previousClusterNo to the value of i
                    collected = list() # initialize collected as a blank list
                previousClusterNo = int(line.replace(">Cluster ","").replace("\n","").replace("",""))+previousClusterB # set previousClusterNo equal to the cluster number from that line of the .txt.clstr file plus the previousClusterB number
            else: # if not a line that has cluster number ('0	107aa, >transposon7735... at 91.59%')
                collected.append(int(line.split(">transposon")[1].split(".")[0])) # append transposon number to collected list
            line = f.readline() # next line
        f.close()
        clusterData[previousClusterNo] = collected # set key value previousClusterNo in clusterData equal to list of collected transposon numbers from previous cluster file
        clusterFileAss[previousClusterNo] = i # set key value previousClusterNo in clusterFileAss equal to i (file number)
        collected = list() # set collected equal to blank list
        previousClusterB = previousClusterNo+1 # go to next cluster number
    print("Collecting clustering results finished...")

    # Sort Transposon Sequences by size
    os.system("seqkit sort -n "+os.path.join(folderOut,"TransposonSequences.fasta")+" > "+os.path.join(folderOut,"TransposonSequences_sort2.fasta")) # sort the clustered transposon sequences by size
    print("Transposon Sequences sorted...")
    print("cluster data->")
    print(clusterData)

    # Create a fasta file and save cluster sequences
    transpIDs = list() # initialize transpIDs as blank list
    for key in list(clusterData.keys()): # for every cluster number key in the list of clusterData keys
        if len(clusterData[key]) != 0:
            cluster = clusterData[key] # set cluster equal to the collected list value from the clusterData dict key
            transpIDs.append(["transposon"+str(cluster[0]),key]) # ERROR # append the transpIDs list with 'transposon transposon# cluster#'
        else:
            pass
    transpIDs.sort() # sort the transposon ids low to high
    fileOut = os.path.join(folderOut,"ClusterSequences.fasta")
    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences_sort2.fasta"), "fasta")
    f = open(fileOut,"w+")
    for t in transpIDs:
        print(t)
        found = False
        for r in records:
            if(r.id==t[0]):
                seq = r.seq
                f.write(">cluster"+str(t[1])+"\n"+str(seq)+"\n")
                found = True
                break
        if(not found):
            print("ERROR",t)
    f.close()

    # Create a Table with cluster and transposon assignment
    fileOut = os.path.join(folderOut,"ClusterSequencesTransposons.txt")
    f = open(fileOut,"w+")
    for key in clusterData:
        f.write(">cluster"+str(key)+"\n"+str(clusterData[key])+"\n")
    f.close()

    ## Create a Fasta with a sequence for each cluster at begin and end of file to check if there is potential in clustering these together
    #fileCluAss = {}
    #for i in range(1,fileCounter):
    #    fileCluAss[i] = list()
    #    for ass in clusterFileAss:
    #        if(clusterFileAss[ass]==i):
    #            fileCluAss[i].append(ass)
    #transpIDs = list()
    #for i in range(1,fileCounter):
    #    if(len(fileCluAss[i])>1):
    #        firstCluster = fileCluAss[i][0]
    #        lastCluster  = fileCluAss[i][-1]
    #        transpIDs.append([firstCluster,clusterData[firstCluster][0]])
    #        transpIDs.append([lastCluster,clusterData[lastCluster][0]])
    #    else:
    #        firstCluster = fileCluAss[i][0]
    #        transpIDs.append([firstCluster,clusterData[firstCluster][0]])
    #transpIDs.sort()
    #fileOut = os.path.join(folderOut,"TransposonSequences_ClusterSeqs.fasta")
    #f = open(fileOut,"w+")
    #for t in transpIDs:
    #    print(t)
    #    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences.fasta"), "fasta")
    #    for r in records:
    #        if(r.id=="transposon"+str(t[1])):
    #            seq = r.seq
    #            f.write(">cluster"+str(t[0])+"\n"+str(seq)+"\n")
    #f.close()

    # Cluster again on them
    #os.system("cd-hit -i "+os.path.join(folderOut,"TransposonSequences_ClusterSeqs.fasta")+" -o "+os.path.join(folderOut,"Results_TransposonSequences_ClusterSeqs.txt"))

    ## Collect cluster results  2
    #clusterData2 = list()
    #collected = list()
    #f = open(os.path.join(folderOut,"Results_TransposonSequences_ClusterSeqs.txt.clstr"),"r")
    #line = f.readline()
    #while line!="":
    #    if(line.startswith(">")):
    #        if(len(collected)==0):
    #            pass
    #        else:
    #            clusterData2.append(collected)
    #            collected = list()
    #    else:
    #        collected.append(int(line.split(">cluster")[1].split(".")[0]))
    #    line = f.readline()
    #f.close()
    #clusterData2.append(collected)
    #collected = list()
    #print("Collecting clustering results finished...")
