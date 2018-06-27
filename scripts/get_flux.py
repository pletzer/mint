import re
import sys

for line in open(sys.argv[1]).readlines():
	m = re.match(r'Flux across the line: (\S+)', line)
	if m:
		print m.group(1)
		break

