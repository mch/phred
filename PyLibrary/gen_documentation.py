# Generates documentation when run from Phred
import sys

try:
    import Phred

except ImportError:
    print "This script must be executed by Phred."
    sys.exit()

try:
    import pydoc

except ImportError:
    print "Pydoc must be installed on your system to generate Phred's documentation. "
    sys.exit()

textoutput = pydoc.text.document(Phred)
fp = open('phred_python.txt', 'w')
fp.write(textoutput)
fp.close()

htmloutput = pydoc.html.document(Phred)
fp = open('phred_python.html', 'w')
fp.write(htmloutput)
fp.close()

print "Documentation for Phred's python interface has been written to "
print "phred_python.txt and phred_python.html."
