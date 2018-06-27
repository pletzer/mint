import re
import sys

for line in open(sys.argv[1]).readlines():
	m = re.match(r'Grid: no of cells (\S+)', line)
	if m:
		print m.group(1)
		break

