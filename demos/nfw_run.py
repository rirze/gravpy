from gravpy import gravlens
from gravpy.models import nfw

pmodelargs = [nfw(2,0,0,0.2,0,0.5)]
ppolargs = [[(0,0),0.9,10,42]]
pcarargs = [[-2.5,2.5,0.5],[-2.5,2.5,0.5]] # lower bound, upper bound, initial spacing (two sets--for x and y axes)

example = gravlens(pcarargs, ppolargs, pmodelargs)
example.run()
