# LTR Retriever integration
The following changes were made to the TEUlt source code to integrate results from LTR Retriever

---

# LTR Retriever parsing function
- [ ] add to AnnotationParsing.py

```python
# pathResDir, fastaFile, targetGFFFile, targetFastaFile
def parseLtrRetriever(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse ltrRetriever outputs...")
    transposonAnnotations = {}
    folder = os.listdir(pathResDir)[0]
    files  = os.listdir(os.path.join(pathResDir,folder))
    selFile = ""
    for f in files:
        if(f.endswith(".gff")):
            selFile = os.path.join(pathResDir,folder,f)
    f = open(selFile,"r")
    line = f.readline()
    counter = 0
    while line!="":

        if line.startswith("#") == False:
            transposons = line.split("\t")
            chrom = transposons[0]
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()

            # print(counter)
            annoSoftware = "ltrRetriever"
            features     = "transposon"
            strand = transposons[6]
            start  = int(transposons[3])
            end    = int(transposons[4])
            start,end = correctPositions(start, end)
            score  = "."
            phase  = "."
            transpNr = str(counter)
            comments = transposons[8].split(";")
#             print(comments)

            component_id = comments[0].split('=')[1]

            if(component_id[0:13] == 'repeat_region'):
                features     = "repeat_region"
                attributes = "Repeat region of transposon "+str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)

            elif(component_id[0:5] == 'LTRRT'):
                features     = "transposon"
                attributes = str(comments)
                counter += 1
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)

            elif(component_id[0:4] == 'lLTR'):
                features     = "ltr"
                attributes = "Left LTR of transposon "+str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)

            elif(component_id[0:4] == 'rLTR'):
                features     = "ltr"
                attributes = "Right LTR of transposon "+str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)

            elif( component_id[0:4] == 'lTSD'):
                features     = "tsd"
                attributes = "Left TSD of transposon "+str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)

            elif(component_id[0:4] == 'rTSD'):
                features     = "tsd"
                attributes = "Right TSD of transposon "+str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)


        line = f.readline()
    f.close()

    # Print results to Annotation file

    f = open(targetGFFFile,"w+")
    f.write("# ltrRetriever Annotation")
    # f.write('\n')
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()

#     Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

```

# LTR Retriever check function
- [ ] add to AnnotationChecker.py

```python
def checkLtrRetriever(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "ltrRetriever")):
        if(getFile(os.path.join(projectFolderPath, "ltrRetriever"),".gff3")!=""):
            return True
#        if(path.isfile(os.path.join(projectFolderPath, "ltrHarvest", "result.txt"))):
#            return True
    return False
```

# adding LTR Retriever to lists and import statements
- [ ] append 'ltrRetriever' to list in line 94 in DuplicateFilterA.py
- [ ] append 'ltrRetriever' to list in line 58 in DuplicateFilterB.py
- [ ] append 'ltrRetriever' to validSoftwares list in line 72
- [ ] append checkLtrRetriever to import statement in line 11 in AnnotationParser.py