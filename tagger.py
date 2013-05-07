#!/usr/bin/env python2.7

#Written by Tristan Chong on 10/23/12 for UW LING 570

import re
import math

trainFile = open("korean-training.txt","r")
trainFileRead = trainFile.read()
trainFile.close()

morphTagPairList = []
morphTagTupleList = []
morphList = []
tagList = []
tagAsFirstElementInBigramList = []
tagBigramList = []

# Make a list of all substrings matching the format "morpheme/TAG"
morphTagPairList = re.findall(r"(\^EOS|[^\s\+]+/[^\s\+]+)",trainFileRead)

# Make a list of tuples where the 1st element is "morpheme" & the 2nd is "TAG"
for pair in morphTagPairList:
	if "/" in pair:
		morphTagTupleList.append(pair.split("/"))
	else:
		morphTagTupleList.append((pair,"<s>"))

# Make 4 lists: morphemes, tags, bigrams ("TAG1,TAG2"), and tags that are the 
# first element in a bigram
prevtag = "<s>"
for morph,tag in morphTagTupleList:
	morphList.append(morph)
	tagList.append(tag)
	tagBigramString = str(prevtag + "," + tag)
	tagBigramList.append(tagBigramString)
	tagAsFirstElementInBigramList.append(prevtag)
	prevtag = tag

# Make 2 lists: unique morphemes and unique tags
morphUniqList = []
for morpheme in set(morphList):
	morphUniqList.append(morpheme)
tagUniqList = []
for tag in set(tagList):
	tagUniqList.append(tag)

# Tally all tag bigram sequences
tagBigramTally = {}
for bigram in tagBigramList:
	tagBigramTally[bigram] = tagBigramTally.setdefault(bigram,0) + 1

# Tally all tags that occur as the first element in a bigram
tagAsFirstElementInBigramTally = {}
for tag in tagAsFirstElementInBigramList:
	tagAsFirstElementInBigramTally[tag] = tagAsFirstElementInBigramTally.setdefault(tag,0) + 1

# Create transition probability matrix
tProbMatrix = [[0 for x in range(len(tagUniqList))] for y in range(len(tagAsFirstElementInBigramTally))]
rowNum = 0
for firstTag in tagUniqList:
	firstTagCount = tagAsFirstElementInBigramTally[firstTag]
	colNum = 0
	for secondTag in tagUniqList:
		bigramString = firstTag + "," + secondTag
		if tagBigramTally.has_key(bigramString):
			bigramCount = tagBigramTally[bigramString]
		else:
			bigramCount = 0
		tProb = float(bigramCount + 1) / float(firstTagCount + len(tagAsFirstElementInBigramTally))
		tProbMatrix[rowNum][colNum] = "%0.6f" % math.log(tProb,2)
		colNum += 1
	rowNum += 1

# Output transition probability matrix
tProbMatrixFile = open("t_prob_matrix.txt","w+")
colLabelString = "\t" + "\t".join(secondTag for secondTag in tagUniqList) + "\n"
tProbMatrixFile.write(colLabelString)
rowNum = 0
for row in tProbMatrix:
	rowLabelString = tagUniqList[rowNum] + "\t"
	tProbMatrixFile.write(rowLabelString)
	for col in row:
		colString = str(col) + "\t"
		tProbMatrixFile.write(colString)
	tProbMatrixFile.write("\n")
	rowNum += 1
tProbMatrixFile.close()

# Tally all morphemes
morphTally = {}
for morpheme in morphList:
	morphTally[morpheme] = morphTally.setdefault(morpheme,0) + 1

# Replace once-occurring morphemes with UNK, reflect these changes in unique 
# morpheme list, and tally all morpheme/tag pairs
pairUnkTally = {}
for pair in morphTagPairList:
	if "^EOS" not in pair:
		morpheme = re.match(r"([^\s\+]+)/",pair)
		if morphTally[morpheme.group(1)] < 2:
			pairUnk = re.sub(r"[^\s\+]+/","UNK/",pair)
			pairUnkTally[pairUnk] = pairUnkTally.setdefault(pairUnk,0) + 1
			morphUniqList.remove(morpheme.group(1))
		else:
			pairUnkTally[pair] = pairUnkTally.setdefault(pair,0) + 1
morphUniqList.append("UNK")

# Tally all tags
tagNoEosTally = {}
for tag in tagList:
	if "<s>" not in tag:
		tagNoEosTally[tag] = tagNoEosTally.setdefault(tag,0) + 1

# Create emission probability matrix
eProbMatrix = [[0 for x in range(len(tagNoEosTally))] for y in range(len(morphUniqList))]
rowNum = 0
for morpheme in morphUniqList:
	colNum = 0
	for tag in tagNoEosTally:
		morphTagPair = str(morpheme) + "/" + str(tag)
		if morphTagPair in pairUnkTally:
			eProb = float(pairUnkTally[morphTagPair]) / float(tagNoEosTally[tag])
			eProbMatrix[rowNum][colNum] = "%0.6f" % math.log(eProb,2)
		else:
			eProbMatrix[rowNum][colNum] = "-inf"
		colNum += 1
	rowNum += 1

# Output emission probability matrix
eProbMatrixFile = open("e_prob_matrix.txt","w+")
colLabelString = "\t" + "\t".join(tag for tag in tagNoEosTally) + "\n"
eProbMatrixFile.write(colLabelString)
rowNum = 0
for row in eProbMatrix:
	rowLabelString = morphUniqList[rowNum] + "\t"
	eProbMatrixFile.write(rowLabelString)
	for col in row:
		colString = str(col) + "\t"
		eProbMatrixFile.write(colString)
	eProbMatrixFile.write("\n")
	rowNum += 1
eProbMatrixFile.close()

# Format testing file
testFile = open("korean-testing.txt","r")
formattedTestFile = open("korean-testing-formatted.txt","w+")
for line in testFile.readlines():
	line = line.rstrip()
	if "^EOS" in line:
		formattedTestFile.write("\n")
	else:
		if line == "+":
			string = "+ "
		elif "+" in line:
			morphemeList = line.split("+")
			string = " ".join(morpheme for morpheme in morphemeList) + " "
		else:
			string = line + " "
		formattedTestFile.write(string)
testFile.close()
formattedTestFile.close()

# Import and run Viterbi algorithm, create list of morphemes in test file
from Viterbi import Viterbi
v = Viterbi()
formattedTestFile = open("korean-testing-formatted.txt","r")
viterbiOutputFile = open("viterbi-out.txt","w+")
testMorphemeList = []

for line in formattedTestFile.readlines():
	line = line.rstrip()
	string = v.tag(line)
	viterbiOutputFile.write(string)
	morphemes = line.split(" ")
	for morpheme in morphemes:
		testMorphemeList.append(morpheme)
formattedTestFile.close()
viterbiOutputFile.close()

# Create list of tags in Viterbi test output
viterbiOutputFile = open("viterbi-out.txt","r")
testTagList = []
for line in viterbiOutputFile.readlines():
	line = line.rstrip()
	tags = line.split(" ")
	for tag in tags:
		testTagList.append(tag)
viterbiOutputFile.close()

# Zip the test results into a list of tuples and create a new list separating 
# morphemes and tags with "/"
testPairTupleList = zip(testMorphemeList,testTagList)
testPairList = []
for morpheme,tag in testPairTupleList:
	testPairListString = morpheme + "/" + tag
	testPairList.append(testPairListString)

# Create a list of tuples containing the morphemes with their correct tags 
# (gold standard)
taggedTestFile = open("korean-testing-tagged.txt","r")
taggedPairList = re.findall(r"([^\s\+]*[^\s]/[^\s\+]*)",taggedTestFile.read())
taggedTestFile.close()

# Zip the gold standard pairs with the test result pairs
taggedTestZipList = zip(taggedPairList,testPairList)

# Count the number of correctly-predicted tags
correctCount = 0
for tagged,test in taggedTestZipList:
	if tagged == test:
		correctCount += 1

# Generate the output file, formatted like the morpheme/tag+morpheme/tag 
# columns in korean-testing-tagged.txt, noting discrepancies with **TAG**
taggedTestFile = open("korean-testing-tagged.txt","r")
outputFile = open("output.txt","w+")
for line in taggedTestFile.readlines():
	line = line.rstrip()
	segments = line.split("\t")
	for segment in segments:
		if "/" in segment:
			outputString = ""
			if "+/" in segment:
				if taggedTestZipList[0][0] in segment:
						if taggedTestZipList[0][0] == taggedTestZipList[0][1]:
							outputString += taggedTestZipList[0][1]
						else:
							tag = re.search(r"([^\s]*)/([^\s]*)",taggedTestZipList[0][1])
							outputString += tag.group(1) + "/**" + tag.group(2) + "**"
						taggedTestZipList.pop(0)
			elif "+" in segment:
				splitPlus = segment.split("+")
				pairs = []
				for item in splitPlus:
					pairs.append(item)
					pairs.append("+")
				pairs.pop()
				for pair in pairs:
					if pair == "+":
						outputString += "+"
					else:
						if taggedTestZipList[0][0] in pair:
							if taggedTestZipList[0][0] == taggedTestZipList[0][1]:
								outputString += taggedTestZipList[0][1]
							else:
								tag = re.search(r"([^\s]*)/([^\s]*)",taggedTestZipList[0][1])
								outputString += tag.group(1) + "/**" + tag.group(2) + "**"
							taggedTestZipList.pop(0)
			else:
				if taggedTestZipList[0][0] in segment:
					if taggedTestZipList[0][0] == taggedTestZipList[0][1]:
						outputString += taggedTestZipList[0][1]
					else:
						tag = re.search(r"([^\s]*)/([^\s]*)",taggedTestZipList[0][1])
						outputString += tag.group(1) + "/**" + tag.group(2) + "**"
					taggedTestZipList.pop(0)
			outputString += "\n"
			outputFile.write(outputString)

# Add morpheme count and % accuracy to the output file
morphemeCount = len(testMorphemeList)
percentAccuracy = 100*float(correctCount)/float(morphemeCount)
outputString = "\nTotal # of morphemes evaluated:\t%i\nAccuracy:\t%0.2f" % (morphemeCount,percentAccuracy) + "%"
outputFile.write(outputString)

taggedTestFile.close()
outputFile.close()
