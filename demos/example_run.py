from gravpy import gravlens
from pycallgraph import PyCallGraph,Config,GlobbingFilter
from pycallgraph.output import GraphvizOutput
from gravpy.models import SIE,alpha,nfw
############################# Actual parameters needed for evaluation
pmodelargs = [
    SIE(0.75,-0.45,-0.5,0.1,20,0.1),
    SIE(0.5,1.5,1.5,0,0.2,0.0),
    SIE(0.1,0.25,0.25,0.1,0,0),
    SIE(0.3,-0.4,0.6,0.5,45,0)
] #models' args, - look in models.py for parameter order (simliar to gravlens param order minus the two shear parameters)
ppolargs = [[(-0.45,-0.5),0.9,10,42],
            [(1.5,1.5),0.5,10,42]] # center position (coordinate pair), outer radius, number of divisions in radius, number of divisions in angle (for 360 degrees)
pmodelargs2 = [alpha(1,0,0,0.2,0,0.,0.75)]
pmodelargs3 = [nfw(2,0,0,0.1,0,0.5)]
ppolargs2 = [[(0,0),0.9,10,42]]
pcarargs = [[-2.5,2.5,0.5],[-2.5,2.5,0.5]] # lower bound, upper bound, initial spacing (two sets--for x and y axes)
pimage = [0.25,0.25] #image location -if- we want to specify
##############################

# If we want an image with runtime stats, set bool to true, otherwise runs the statement in the else branch.
callgraph = False
filepath = 'runs/alphafintelparallel.png' #where we want to save the output image with the runtime breakdown
config = Config()
config.trace_filter = GlobbingFilter(exclude=[
    'num*',
    'scipy*',
    
    'pycallgraph*'
])
if callgraph:
    with PyCallGraph(output=GraphvizOutput(output_file=filepath),config=config):
        example = gravlens(pcarargs,ppolargs2,pmodelargs2, image=pimage, show_plot=False, include_caustics=False)
        example.run()
else:
    example = gravlens(pcarargs,ppolargs2,pmodelargs3,image=pimage,include_caustics=True,logging_level='info' )
    example.run()
    

    
