### Protein Dimer Genomic Binding Site Search ###

### Search Settings

## Forward and Reverse Binding Sites for search
FW_1 = "TCCCG"
FW_2 = "CGGGA"
RV_1 = "AGGGC"
RV_2 = "GCCCT"

## Annotated GenBank File Name
genomeFile = "NA1000.gb"

## Search Interval (in BP)
searchInterval = 100

# Program Constants, do not change!

GENESTART = 0
GENEEND = 1
DIRECTION = 2
SEARCHSTART = 3
SEARCHEND = 4
FOUND = 5
FOUNDINDEX = 6

FORWARD = True


dataMatrix = []
genomeLibrary = []
genomeLength = 0

def getLocii():
	text_file = open(genomeFile, "r")
	lines = text_file.readlines()

	for line in lines:
		if line.find("  gene  ") != -1:
			a = 8 + line.find("  gene  ")
			newLine = line[a:].strip()
			dataMatrix.append(readAnnotations(newLine))
	text_file.close()
	return

def buildGenomeLibrary():
	global genomeLibrary
	global genomeLength

	currString = ""

	text_file = open(genomeFile, "r")
	lines = text_file.readlines()

	buildLibrary = False

	for line in lines:
		if buildLibrary:
			DNAline = line[10:-1].split(" ")

			for segment in DNAline:
				currString = currString + segment.upper()
				genomeLength += 10

			if len(currString) == 600000:
				genomeLibrary.append(currString)
				currString = ""

		elif line.find("ORIGIN  ") != -1:
			buildLibrary = True
	genomeLibrary.append(currString)

	genomeLibrary.append(currString)

	text_file.close()
	return


def readAnnotations(line):
	annotation = []

	if line.find("complement") != -1:
		line = line[line.find('(')+1:line.find(')')]
		fw = True
	else:
		fw = False
	indicies = line.split("..")
	annotation.append(int(indicies[0]))
	annotation.append(int(indicies[1]))
	annotation.append(fw)
	return annotation

def processLocii():
	for gene in dataMatrix:
		if gene[DIRECTION] == FORWARD:
			gene.append(gene[GENEEND])
			gene.append(gene[GENEEND] + searchInterval)
		else:
			gene.append(gene[GENESTART] - searchInterval)
			gene.append(gene[GENESTART])


def getSequence(start, end):

	sequence = genomeLibrary[start//600000]
	index = start % 600000
	index2 = index + (end - start)

	if index2 > 600000:
		print "problem"

	sequence = sequence[index:index2]

	return sequence


def search():
	for gene in dataMatrix:
		searchFragment = getSequence(gene[SEARCHSTART],gene[SEARCHEND])

		if gene[DIRECTION] == FORWARD:
			firstSearch = FW_1
			secondSearch = FW_2
		else:
			firstSearch = RV_1
			secondSearch = RV_2

		index = searchFragment.find(firstSearch)
		if index == -1:
			gene.append(False)
			gene.append(-1)

		else:
			# print gene
			# print genome[gene[SEARCHSTART] + index:gene[SEARCHSTART] + index + 5]

			if searchFragment[index + 7: index + 12].find(secondSearch) != -1:
				gene.append(True)
				gene.append(gene[SEARCHSTART] + index)
			else:
				gene.append(False)
				gene.append(-1)
	return

def saveResultsGB():
	resultsFile = open("results.gb", "w")
	resultsFile.write("LOCUS       CP001340 features                                       11-AUG-2014 \nUNIMARK     CP001340 features \nFEATURES             Location/Qualifiers\n")
	currNum = 0

	for gene in dataMatrix:
		if gene[FOUND]:
			resultsFile.write("     bindsite        " + str(gene[FOUNDINDEX]) + ".." + str(gene[FOUNDINDEX]+12) + "\n")
			resultsFile.write('''                     /ugene_name="name''' + str(currNum) + '''"\n''')
			resultsFile.write('''                     /ugene_group="bindgroup"\n''')
			currNum = currNum + 1
	resultsFile.write("//")		
	resultsFile.close()
	return

def quote(string):
	return ('''"''' + string + '''"''')

def saveResultsCSV():
	currNum = 0

	resultsFile = open("results.csv", "w")
	resultsFile.write('''"Group","Name","Start","End","Length","Complementary"\n''')
	for gene in dataMatrix:
		if gene[FOUND]:
			resultsFile.write("" + quote("bindSite") + "," + quote("site" + str(currNum)) + "," + quote(str(gene[FOUNDINDEX])) + "," + quote(str(gene[FOUNDINDEX] + 12)) + "," + quote(str(12)) + "," + quote("no") + "\n")
			currNum = currNum + 1
	resultsFile.close()
	return


buildGenomeLibrary()
getLocii()
processLocii()
search()
saveResultsCSV()



