#!/usr/bin/python

import random, sys, math

nFiles = 1
if len(sys.argv) < 4:
    print 'usage: python generateRandomData.py filename d N [nFiles]'
    sys.exit(1)
else:
    filename = sys.argv[1]
    d = int(sys.argv[2])
    N = int(sys.argv[3])
    if len(sys.argv) > 4:
      nFiles = int(sys.argv[4])

MIN = -100
MAX = 100

valuesPerFile = int(math.ceil(float(N) / nFiles))

coeff = [random.uniform(-10, 10) for i in xrange(d)]
for n in xrange(nFiles):
  if nFiles > 1:
    name = '%s_%d.dat' % (filename, n)
  else:
    name = '%s.dat' % filename
    
  f = open(name, 'w')
  
  if nFiles == 1 or n < nFiles - 1 or N % valuesPerFile == 0:
    nValues = valuesPerFile
  else:
    nValues = N % valuesPerFile
  
  f.write('%d %d\n' % (d, nValues))
  
  for i in xrange(nValues):
      val = [random.randint(MIN, MAX) for i in xrange(d)]
      label = sum([coeff[i]*val[i] for i in xrange(d)])
      label += random.gauss(0,100);
      
      for j in xrange(d):
        f.write('%d ' % val[j])
      f.write('%d\n' % label)
  f.close()
