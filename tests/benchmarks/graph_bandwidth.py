from pyx import *
import sys

g = graph.graphxy(width=10, key=graph.key.key(pos="tl"),
                  x=graph.axis.log(title="Message Size (bytes)"),
                  y=graph.axis.linear(title="Bandwidth (MB/s)"))

colours = [color.rgb.red, color.rgb.green, color.rgb.blue]
titles = ["xy plane", "yz plane", "zx plane"]
i=0

for f in sys.argv[1:]:
    print "graphing %s..." % f
    data = graph.data.file(f, title=titles[i],
                           x="$1**2", y="$3/1e6")
    g.plot(data, [graph.style.line([ colours[i] ])])
    i = i + 1

g.writeEPSfile("graph.eps")
