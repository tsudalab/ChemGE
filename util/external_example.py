# coding: utf-8
# This is an example of a script called from RbtExternalTransform,
# and only queries the vectors csv

import sys
import numpy as np

class RDock():
  '''Interface to RbtExternalTransform
  '''
  def __init__(self):
    self.out = sys.stdout
    self.vec = None
    self.score = None

  def __del__(self):
    self.out.write('\0')
    self.out.close()

  def load(self):
    '''Load the population data
  
    The csv file output by rDock is loaded, then
    vectors and scores are returned separately.
    '''
    if (self.vec is None):
      csv = np.loadtxt(sys.argv[1], skiprows=1, delimiter=',')
      self.vec = csv[:,:-1]
      self.score = csv[:,-1]
    return self.vec, self.score

  def query(self, v):
    '''Query an evaluation of vector

    Given a feature vector, a query is sent to rDock process
    and its score is returned.
    '''
    sys.stdout.write(','.join(['%.6f'%x for x in v]))
    sys.stdout.write('\n')
    sys.stdout.flush()
    score = sys.stdin.readline()
    return np.float32(score)

  def log(self, s, bl=True):
    '''Logging

    Because stdin/stdout are used to communicate with rDock,
    the logging messages are displayed through stderr.
    '''
    sys.stderr.write(str(s))
    if bl:
      sys.stderr.write('\n')

# entry point
rdock = RDock()
vec, _ = rdock.load()
for v in vec:
  score = rdock.query(v)
  rdock.log('%.6f'%score)
