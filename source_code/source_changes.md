# Source Code Changes
The following were other changes made to the TEUlt source code

## transposonclustering.py changes
transposonclustering.py changes to avoid this error:

```python
Traceback (most recent call last):
  File "/home/calbers/.conda/envs/transposon_annotation_reasonaTE/share/TransposonAnnotator_reasonaTE/TransposonAnnotator.py", line 160, in <module>
    doTransposonClustering(os.path.join(arg1,arg2,"transposonCandB"),os.path.join(arg1,arg2,"transposonCandC"),os.path.join(arg1,arg2,"sequence.fasta"))
  File "/gpfs/home/calbers/.conda/envs/transposon_annotation_reasonaTE/share/TransposonAnnotator_reasonaTE/TransposonClustering.py", line 187, in doTransposonClustering
    transpIDs.append(["transposon"+str(cluster[0]),key]) # ERROR # append the transpIDs list with 'transposon transposon# cluster#'
IndexError: list index out of range 
```

original (lines 184-188):

```python
# Create a fasta file and save cluster sequences
transpIDs = list() # initialize transpIDs as blank list
for key in list(clusterData.keys()): # for every cluster number key in the list of clusterData keys
    cluster = clusterData[key] # set cluster equal to the collected list value from the clusterData dict key
    transpIDs.append(["transposon"+str(cluster[0]),key]) # ERROR # append the transpIDs list with 'transposon transposon# cluster#'

```

changed (lines 184-191):

```python
# Create a fasta file and save cluster sequences
transpIDs = list() # initialize transpIDs as blank list
for key in list(clusterData.keys()): # for every cluster number key in the list of clusterData keys
      if len(clusterData[key]) != 0:
          cluster = clusterData[key] # set cluster equal to the collected list value from the clusterData dict key
          transpIDs.append(["transposon"+str(cluster[0]),key]) # ERROR # append the transpIDs list with 'transposon transposon# cluster#'
      else:
          pass

```