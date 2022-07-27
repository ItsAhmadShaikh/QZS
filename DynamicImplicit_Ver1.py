from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqus import getInput
from odbAccess import *
import numpy as np
import math
import regionToolset
from scipy.optimize import fsolve
import numpy as np
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


t = 0.5e-3;
h = -0.67e-3;
# h = -1*t*math.sqrt(2);
r1 = 11e-3;
r2 = 17.25e-3;

L1 = [r1-r2, -h, (r2-r1)*r1];
L2 = [r1-r2, -h, h*h + (r2-r1)*r2];
L3 = [h, r1-r2, (-h*r1)+0.5*t*math.sqrt(h*h*r1*r1+(r2-r1)*(r2-r1))];
L4 = [h, r1-r2, (-h*r1)-0.5*t*math.sqrt(h*h*r1*r1+(r2-r1)*(r2-r1))];

p1 = [ (L1[1]*L3[2] - L3[1]*L1[2]) / (L1[0]*L3[1]-L3[0]*L1[1]) , (L3[0]*L1[2]-L1[0]*L3[2])/(L1[0]*L3[1]-L3[0]*L1[1]) ]
p2 = [ (L1[1]*L4[2] - L4[1]*L1[2]) / (L1[0]*L4[1]-L4[0]*L1[1]) , (L4[0]*L1[2]-L1[0]*L4[2])/(L1[0]*L4[1]-L4[0]*L1[1]) ]
p3 = [ (L2[1]*L4[2] - L4[1]*L2[2]) / (L2[0]*L4[1]-L4[0]*L2[1]) , (L4[0]*L2[2]-L2[0]*L4[2])/(L2[0]*L4[1]-L4[0]*L2[1]) ]
p4 = [ (L2[1]*L3[2] - L3[1]*L2[2]) / (L2[0]*L3[1]-L3[0]*L2[1]) , (L3[0]*L2[2]-L2[0]*L3[2])/(L2[0]*L3[1]-L3[0]*L2[1]) ]



mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((206000000000.0, 
    0.3), ))
mdb.models['Model-1'].materials['Material-1'].Density(table=((7900.0, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-1', name=
    'Section-1', thickness=None)


#Geometry

mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.03)
mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues(
    decimalPlaces=4, viewStyle=AXISYM)
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -0.015), point2=(0.0, 0.015))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.0), 
    ))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(p1[0], p1[1]), 
    point2=(p2[0], p2[1]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(p2[0], p2[1]), 
    point2=(p3[0], p3[1]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(p3[0], p3[1]), 
    point2=(p4[0], p4[1]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(p4[0], p4[1]), 
    point2=(p1[0], p1[1]))
mdb.models['Model-1'].Part(dimensionality=AXISYMMETRIC, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

#Sets
mdb.models['Model-1'].parts['Part-1'].Set(faces=
    mdb.models['Model-1'].parts['Part-1'].faces.findAt(((r1,0,0), 
    )), name='All')
mdb.models['Model-1'].parts['Part-1'].Set(name='P1', vertices=
    mdb.models['Model-1'].parts['Part-1'].vertices.findAt(((p1[0], p1[1],0), )))
mdb.models['Model-1'].parts['Part-1'].Set(name='P2', vertices=
    mdb.models['Model-1'].parts['Part-1'].vertices.findAt(((p2[0], p2[1],0), )))
mdb.models['Model-1'].parts['Part-1'].Set(name='P3', vertices=
    mdb.models['Model-1'].parts['Part-1'].vertices.findAt(((p3[0], p3[1],0), )))
mdb.models['Model-1'].parts['Part-1'].Set(name='P4', vertices=
    mdb.models['Model-1'].parts['Part-1'].vertices.findAt(((p4[0], p4[1],0), )))
    
    
    
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-1'].sets['All'], sectionName='Section-1', 
    thicknessAssignment=FROM_SECTION)
    

mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=0.0001)
mdb.models['Model-1'].parts['Part-1'].setMeshControls(algorithm=MEDIAL_AXIS, 
    elemShape=QUAD, regions=mdb.models['Model-1'].parts['Part-1'].faces.findAt(
    ((r1,0,0), )))
mdb.models['Model-1'].parts['Part-1'].generateMesh()
mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
    elemCode=CAX4I, elemLibrary=STANDARD), ElemType(elemCode=CAX3, 
    elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Part-1'].faces.findAt(((r1,0,0),)), ))
    
    
#Copy Part
mdb.models['Model-1'].Part(compressFeatureList=ON, mirrorPlane=XZPLANE, name=
    'Part-2', objectToCopy=mdb.models['Model-1'].parts['Part-1'])
    
#Sets
mdb.models['Model-1'].parts['Part-2'].Set(faces=
    mdb.models['Model-1'].parts['Part-2'].faces.findAt(((r1,0,0), 
    )), name='All')
mdb.models['Model-1'].parts['Part-2'].Set(name='P1', vertices=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt(((p1[0], -p1[1],0), )))
mdb.models['Model-1'].parts['Part-2'].Set(name='P2', vertices=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt(((p2[0], -p2[1],0), )))
mdb.models['Model-1'].parts['Part-2'].Set(name='P3', vertices=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt(((p3[0], -p3[1],0), )))
mdb.models['Model-1'].parts['Part-2'].Set(name='P4', vertices=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt(((p4[0], -p4[1],0), )))
    
    
    
mdb.models['Model-1'].parts['Part-2'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-2'].sets['All'], sectionName='Section-1', 
    thicknessAssignment=FROM_SECTION)
    

mdb.models['Model-1'].parts['Part-2'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=0.0001)
mdb.models['Model-1'].parts['Part-2'].setMeshControls(algorithm=MEDIAL_AXIS, 
    elemShape=QUAD, regions=mdb.models['Model-1'].parts['Part-2'].faces.findAt(
    ((r1,0,0), )))
mdb.models['Model-1'].parts['Part-2'].generateMesh()
mdb.models['Model-1'].parts['Part-2'].setElementType(elemTypes=(ElemType(
    elemCode=CAX4I, elemLibrary=STANDARD), ElemType(elemCode=CAX3, 
    elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Part-2'].faces.findAt(((r1,0,0),)), ))
    

    
    
#Assembly
mdb.models['Model-1'].rootAssembly.DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 
    0.0, -1.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
    part=mdb.models['Model-1'].parts['Part-2'])
    
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-1', ), 
    vector=(0.0, -1.6e-3+2*h, 0.0))
    

#Coupling
mdb.models['Model-1'].Coupling(controlPoint=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['P3'], 
    couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
    name='Constraint-1', surface=
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].sets['P3'], u1=ON, 
    u2=ON, ur3=ON)
    
#BC
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].sets['P1'], 
    u1=UNSET, u2=SET, ur3=UNSET)
    
#DynamicImplicit
mdb.models['Model-1'].ImplicitDynamicsStep(alpha=DEFAULT, application=
    TRANSIENT_FIDELITY, initialConditions=ON, initialInc=0.05, maxInc=0.05, maxNumInc=100000
    , name='Step-1', nlgeom=ON, nohaf=OFF, previous='Initial')
    
# mdb.models['Model-1'].StaticRiksStep(dof=2, initialArcInc=0.05, maxArcInc=0.05, 
    # maxNumInc=1000000, maximumDisplacement=-0.003, name='Step-1', nlgeom=ON, 
    # nodeOn=ON, previous='Initial', region=
    # mdb.models['Model-1'].rootAssembly.allInstances['Part-1-1'].sets['P1'])
    
#Gravity
mdb.models['Model-1'].Gravity(comp2=9.81, createStepName='Step-1', 
    distributionType=UNIFORM, field='', name='Gravity')
    
#Disp
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'BC-2', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['P1'], u1=
    UNSET, u2=-0.003, ur3=UNSET)

##Point Mass????
mdb.models['Model-1'].rootAssembly.engineeringFeatures.PointMassInertia(alpha=
    0.0, composite=0.0, mass=0.0113e3, name='Inertia-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['P1'])
    
    
#Output
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(frequency=1)
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(frequency=1)
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name=
    'H-Output-2',frequency= 1, rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.allInstances['Part-1-1'].sets['P1'], 
    sectionPoints=DEFAULT, variables=('U2', ))
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name=
    'H-Output-3',frequency= 1, rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.allInstances['Part-2-1'].sets['P1'], 
    sectionPoints=DEFAULT, variables=('RF2', ))

    
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
   



# ################################   
    
# mdb.Model(name='Model-1-Copy', objectToCopy=mdb.models['Model-1'])
# mdb.models['Model-1-Copy'].ImplicitDynamicsStep(initialInc=0.05, 
    # maintainAttributes=True, maxNumInc=10000, name='Step-1', nlgeom=ON, 
    # previous='Initial')
# mdb.models['Model-1-Copy'].boundaryConditions['BC-2'].setValues(u2=-0.003)
# mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    # explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    # memory=90, memoryUnits=PERCENTAGE, model='Model-1-Copy', modelPrint=OFF, 
    # multiprocessingMode=DEFAULT, name='Job-2', nodalOutputPrecision=SINGLE, 
    # numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    # ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)