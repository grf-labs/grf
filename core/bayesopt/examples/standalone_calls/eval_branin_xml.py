#!/usr/bin/env python

# Imports
try:
    import sys
    import traceback
    import xml.dom.minidom as minidom
    from time import sleep
    from random import randint
    from math import sin, cos, pi
except ImportError, err:
    print "ERROR: Cannot load module: %s" % (err)
    sys.exit(2)
    
def sqr(x):
    return x*x

def branin(x, y):
    # sqr(y-(5.1/(4*sqr(M_PI)))*sqr(x)
    #       +5*x/M_PI-6)+10*(1-1/(8*M_PI))*cos(x)+10;
    return sqr(y-(5.1/(4*sqr(pi)))*sqr(x)+5*x/pi-6)+10*(1-1/(8*pi))*cos(x)+10
#} end branin

# main function implementation
def main(argv):
    if len(argv) != 2:
        print "ERROR: 2 arguments required"
        return

    inputFilename = argv[0];
    outputFilename = argv[1];

    # Read parameters from xml document
    try:
        doc = minidom.parse(inputFilename)
    except:
        print "Unable to open and parse input xml file: " + inputFilename
        return

    node = doc.getElementsByTagName("node")[0]
    x1 = float(node.getAttribute("par1"))
    x2 = float(node.getAttribute("par2"))
    
    # Some random delay
    ms = randint(0,500)/1000.0
    sleep(ms)
    print "Python: delay = " + str(ms) + " ms."
    
    # Output results into file
    # TODO (Javier): Change into XML
    y = branin(float(x1)*15-5,float(x2)*15)
    f1=open(outputFilename, "w+")
    f1.write("y="+str(y))
#}end main

# Execution   
try:
    if __name__ == '__main__' :
        main(sys.argv[1:])
except Exception, e:
    tb = sys.exc_info()[2]
    traceback.print_exception(e.__class__, e, tb)
