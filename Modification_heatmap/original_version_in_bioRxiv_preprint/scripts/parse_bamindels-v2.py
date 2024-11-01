import sys

if __name__ == '__main__':
	args = sys.argv[1:]
	indelsFile = open(args[0], 'r')
	lines = indelsFile.readlines()
	outLines = '\n'

	for line in lines:
		indels = line.split('\t')
		Icount = 0
		Dcount = 0

		if len(indels[0]) > 1:
			for indel in indels:
				data = indel.split(':')
				if data[0][0] == '+':
					Icount += int(data[1])
				else:
					Dcount += int(data[1])
		outLines += str(Icount) + '\t' + str(Dcount) + '\n'


destFile = open(args[1], 'w')
destFile.write(outLines)
destFile.close()
