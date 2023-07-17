param = {}
with open('data/param.csv', 'r') as fin:
	while True:
		line = fin.readline()
		if not line: break
		param[line.split(',')[1].strip()] = [
		    float(line.split(',')[2]),
		    float(line.split(',')[3]),
		]

color = {}
with open('data/color.csv', 'r') as fin:
	while True:
		line = fin.readline()
		if not line: break
		color[line.split(',')[0].strip()] = [
		    float(line.split(',')[1]),
		    float(line.split(',')[2]),
		    float(line.split(',')[3]),
		]
