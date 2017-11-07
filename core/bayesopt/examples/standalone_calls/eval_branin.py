#!/usr/bin/env python

# Imports
try:
    import sys
    import traceback
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
    y = branin(float(argv[0])*15-5,float(argv[1])*15)
    f1=open("./results.txt", "w+")
    f1.write("y="+str(y))
#}end main

# Execution   
try:
    if __name__ == '__main__' :
        main(sys.argv[1:])
except Exception, e:
    tb = sys.exc_info()[2]
    traceback.print_exception(e.__class__, e, tb)
