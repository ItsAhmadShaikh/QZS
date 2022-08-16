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

## The following Scrpit is used to automate the 2D axisymetric simulations of disc spring stacks
##Use functions 'RunModel' and 'OutputModel' to automatically fun and create outputs
##Change the 'JobPath; variableif you wish to use the 'OutputModel' command
## Multiple Models are created through the 2D list h_all_models. Each list within the variable currosponts to the heights of each member of a single stack

# h_all_models = [[0.7e-3,0.7e-3,0.7e-3,0.7e-3,0.7e-3,0.7e-3], [0.67e-3,0.675e-3,0.68e-3,0.685e-3,0.69e-3,0.695e-3], [0.715e-3,0.725e-3,0.735e-3,0.745e-3,0.755e-3,0.765e-3]]
h_all_models = [np.ones(6)*0.7e-3,np.ones(6)*0.7e-3,np.ones(6)*0.7e-3];

friction_all_models = [0.05,0.3,0.7];

JobPath = 'D:/OneDrive - Indian Institute of Science/TRINA_QZS/Week4/Simulation/'



N_models = len(h_all_models)

for k in range(N_models):

    ModelName = 'Stack_'+str(k+1)
    JobName = ModelName

    mdb.Model(modelType=STANDARD_EXPLICIT, name=ModelName)
        
    t = 0.5e-3;
    # h = -0.67e-3;
    # h = -1*t*math.sqrt(2);
    # h_all = [0.67e-3, 0.675e-3, 0.68e-3, 0.685e-3, 0.69e-3, 0.695e-3];
    # h_all = [0.65e-3, 0.65e-3, 0.65e-3, 0.65e-3, 0.65e-3, 0.65e-3];
    #h_all = np.ones(2)*t*math.sqrt(2);
    # h_all = [0.67e-3, 0.675e-3];
    h_all = h_all_models[k];
    r1 = 11e-3;
    r2 = 17.25e-3;
    
    friction = friction_all_models[k];

    N_part = len(h_all);
    
    Dy1 = np.zeros(N_part);
    Dy2 = np.zeros(N_part);
    Dy3 = np.zeros(N_part);
    Dy4 = np.zeros(N_part);


    #Material
    mdb.models[ModelName].Material(name='Material-1')
    mdb.models[ModelName].materials['Material-1'].Elastic(table=((206000000000.0, 
        0.3), ))
    mdb.models[ModelName].materials['Material-1'].Density(table=((7900.0, ), ))
    mdb.models[ModelName].HomogeneousSolidSection(material='Material-1', name=
        'Section-1', thickness=None)




    for i in range(N_part):

        PartName = 'Part'+str(i+1)

        h = h_all[i]*math.pow(-1,i+1);

        #Calculations
        L1 = [r1-r2, -h, (r2-r1)*r1];
        L2 = [r1-r2, -h, h*h + (r2-r1)*r2];
        L3 = [h, r1-r2, (-h*r1)+0.5*t*math.sqrt(h*h*r1*r1+(r2-r1)*(r2-r1))];
        L4 = [h, r1-r2, (-h*r1)-0.5*t*math.sqrt(h*h*r1*r1+(r2-r1)*(r2-r1))];

        p1 = [ (L1[1]*L3[2] - L3[1]*L1[2]) / (L1[0]*L3[1]-L3[0]*L1[1]) , (L3[0]*L1[2]-L1[0]*L3[2])/(L1[0]*L3[1]-L3[0]*L1[1]) ]
        p2 = [ (L1[1]*L4[2] - L4[1]*L1[2]) / (L1[0]*L4[1]-L4[0]*L1[1]) , (L4[0]*L1[2]-L1[0]*L4[2])/(L1[0]*L4[1]-L4[0]*L1[1]) ]
        p3 = [ (L2[1]*L4[2] - L4[1]*L2[2]) / (L2[0]*L4[1]-L4[0]*L2[1]) , (L4[0]*L2[2]-L2[0]*L4[2])/(L2[0]*L4[1]-L4[0]*L2[1]) ]
        p4 = [ (L2[1]*L3[2] - L3[1]*L2[2]) / (L2[0]*L3[1]-L3[0]*L2[1]) , (L3[0]*L2[2]-L2[0]*L3[2])/(L2[0]*L3[1]-L3[0]*L2[1]) ]

        Dy1[i] = abs(p1[1])
        Dy2[i] = abs(p2[1])
        Dy3[i] = abs(p3[1])
        Dy4[i] = abs(p4[1])

        #Geometry
        mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=0.03)
        mdb.models[ModelName].sketches['__profile__'].sketchOptions.setValues(
            decimalPlaces=4, viewStyle=AXISYM)
        mdb.models[ModelName].sketches['__profile__'].ConstructionLine(point1=(0.0, 
            -0.015), point2=(0.0, 0.015))
        mdb.models[ModelName].sketches['__profile__'].geometry.findAt((0.0, 0.0))
        mdb.models[ModelName].sketches['__profile__'].FixedConstraint(entity=
            mdb.models[ModelName].sketches['__profile__'].geometry.findAt((0.0, 0.0), 
            ))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(p1[0], p1[1]), 
            point2=(p2[0], p2[1]))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(p2[0], p2[1]), 
            point2=(p3[0], p3[1]))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(p3[0], p3[1]), 
            point2=(p4[0], p4[1]))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(p4[0], p4[1]), 
            point2=(p1[0], p1[1]))
        mdb.models[ModelName].Part(dimensionality=AXISYMMETRIC, name=PartName, type=
            DEFORMABLE_BODY)
        mdb.models[ModelName].parts[PartName].BaseShell(sketch=
            mdb.models[ModelName].sketches['__profile__'])
        del mdb.models[ModelName].sketches['__profile__']

        #Sets
        mdb.models[ModelName].parts[PartName].Set(faces=
            mdb.models[ModelName].parts[PartName].faces.findAt(((r1,0,0), 
            )), name='All')
        mdb.models[ModelName].parts[PartName].Set(name='P1', vertices=
            mdb.models[ModelName].parts[PartName].vertices.findAt(((p1[0], p1[1],0), )))
        mdb.models[ModelName].parts[PartName].Set(name='P2', vertices=
            mdb.models[ModelName].parts[PartName].vertices.findAt(((p2[0], p2[1],0), )))
        mdb.models[ModelName].parts[PartName].Set(name='P3', vertices=
            mdb.models[ModelName].parts[PartName].vertices.findAt(((p3[0], p3[1],0), )))
        mdb.models[ModelName].parts[PartName].Set(name='P4', vertices=
            mdb.models[ModelName].parts[PartName].vertices.findAt(((p4[0], p4[1],0), )))
            
        #Surfaces
            
        x = (r1+r2)/2
        y = -1*(L3[0]*x + L3[2])/L3[1]
        mdb.models[ModelName].parts[PartName].Surface(name='Top', side1Edges=
            mdb.models[ModelName].parts[PartName].edges.findAt(((x, y,0), )))
        
        x = (r1+r2)/2
        y = -1*(L4[0]*x + L4[2])/L4[1]
        mdb.models[ModelName].parts[PartName].Surface(name='Bottom', side1Edges=
            mdb.models[ModelName].parts[PartName].edges.findAt(((x, y,0), )))
            
        #SectionAssignment
        mdb.models[ModelName].parts[PartName].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=
            mdb.models[ModelName].parts[PartName].sets['All'], sectionName='Section-1', 
            thicknessAssignment=FROM_SECTION)
        
        #Mesh
        mdb.models[ModelName].parts[PartName].seedPart(deviationFactor=0.1, 
            minSizeFactor=0.1, size=0.0001)
        mdb.models[ModelName].parts[PartName].setMeshControls(algorithm=MEDIAL_AXIS, 
            elemShape=QUAD, regions=mdb.models[ModelName].parts[PartName].faces.findAt(
            ((r1,0,0), )))
        mdb.models[ModelName].parts[PartName].generateMesh()
        mdb.models[ModelName].parts[PartName].setElementType(elemTypes=(ElemType(
            elemCode=CAX4I, elemLibrary=STANDARD), ElemType(elemCode=CAX3, 
            elemLibrary=STANDARD)), regions=(
            mdb.models[ModelName].parts[PartName].faces.findAt(((r1,0,0),)), ))
            
    r = [r2,r1];

    for i in range(2):    
        PartName = 'Block'+str(i+1);
        #Geometry
        mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=0.03)
        mdb.models[ModelName].sketches['__profile__'].sketchOptions.setValues(
            decimalPlaces=4, viewStyle=AXISYM)
        mdb.models[ModelName].sketches['__profile__'].ConstructionLine(point1=(0.0, 
            -0.015), point2=(0.0, 0.015))
        mdb.models[ModelName].sketches['__profile__'].geometry.findAt((0.0, 0.0))
        mdb.models[ModelName].sketches['__profile__'].FixedConstraint(entity=
            mdb.models[ModelName].sketches['__profile__'].geometry.findAt((0.0, 0.0), 
            ))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(r[i]-0.5e-3, Dy1[0]), 
            point2=(r[i]+0.5e-3, Dy1[0]))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(r[i]-0.5e-3, Dy1[0]+1.7e-3), 
            point2=(r[i]+0.5e-3, Dy1[0]+1.7e-3))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(r[i]-0.5e-3, Dy1[0]), 
            point2=(r[i]-0.5e-3, Dy1[0]+1.7e-3))
        mdb.models[ModelName].sketches['__profile__'].Line(point1=(r[i]+0.5e-3, Dy1[0]), 
            point2=(r[i]+0.5e-3, Dy1[0]+1.7e-3))
        mdb.models[ModelName].Part(dimensionality=AXISYMMETRIC, name=PartName, type=
            DEFORMABLE_BODY)
        mdb.models[ModelName].parts[PartName].BaseShell(sketch=
            mdb.models[ModelName].sketches['__profile__'])
        del mdb.models[ModelName].sketches['__profile__']
        
        #Sets
        mdb.models[ModelName].parts[PartName].Set(faces=
            mdb.models[ModelName].parts[PartName].faces.findAt(((r[i],0.85e-3,0), 
            )), name='All')
            
        #Surface
            
        mdb.models[ModelName].parts[PartName].Surface(name='Top', side1Edges=
            mdb.models[ModelName].parts[PartName].edges.findAt(((r[i], Dy1[0]+1.7e-3,0), )))
        mdb.models[ModelName].parts[PartName].Surface(name='Bottom', side1Edges=
            mdb.models[ModelName].parts[PartName].edges.findAt(((r[i], Dy1[0],0), )))
        
        #SectionAssignment
        mdb.models[ModelName].parts[PartName].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=
            mdb.models[ModelName].parts[PartName].sets['All'], sectionName='Section-1', 
            thicknessAssignment=FROM_SECTION)
            
        #Mesh
        mdb.models[ModelName].parts[PartName].seedPart(deviationFactor=0.1, 
            minSizeFactor=0.1, size=0.0001)
        mdb.models[ModelName].parts[PartName].setMeshControls(algorithm=MEDIAL_AXIS, 
            elemShape=QUAD, regions=mdb.models[ModelName].parts[PartName].faces.findAt(
            ((r1,0,0), )))
        mdb.models[ModelName].parts[PartName].generateMesh()
        mdb.models[ModelName].parts[PartName].setElementType(elemTypes=(ElemType(
            elemCode=CAX4I, elemLibrary=STANDARD), ElemType(elemCode=CAX3, 
            elemLibrary=STANDARD)), regions=(
            mdb.models[ModelName].parts[PartName].faces.findAt(((r[i],0.85e-3,0),)), ))
        
                

            
    #Assembly
    mdb.models[ModelName].rootAssembly.DatumCsysByThreePoints(coordSysType=
        CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 
        0.0, -1.0))
            
    PartName = 'Part'+str(1)
    mdb.models[ModelName].rootAssembly.Instance(dependent=ON, name=PartName +'-'+str(1), 
            part=mdb.models[ModelName].parts[PartName])
            
    adder = 0;
    for i in range(1,N_part):
        
        PartName = 'Part'+str(i+1)
        mdb.models[ModelName].rootAssembly.Instance(dependent=ON, name=PartName+'-'+str(1), 
            part=mdb.models[ModelName].parts[PartName])
            
        
        if i % 2 != 0:
        
            adder = adder + Dy4[i]+Dy3[i-1] + 1.7e-3;
            mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(1), ), 
                vector=(0.0, -1*(adder), 0.0))
                
        else:
        
            adder = adder + Dy2[i]+Dy1[i-1] + 1.7e-3;
            mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(1), ), 
                vector=(0.0, -1*(adder), 0.0))
                
                
                
    adder = 0;
    for i in range(N_part-1):
            
        if i % 2 != 0:
            PartName = 'Block'+str(2);
            mdb.models[ModelName].rootAssembly.Instance(dependent=ON, name=PartName+'-'+str(i+1), 
                part=mdb.models[ModelName].parts[PartName])
        
            adder = adder + Dy4[i]+Dy2[i] + 1.7e-3;
            mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(i+1), ), 
                vector=(0.0, -1*(adder), 0.0))
                
        else:
            PartName = 'Block'+str(1);
            mdb.models[ModelName].rootAssembly.Instance(dependent=ON, name=PartName+'-'+str(i+1), 
                part=mdb.models[ModelName].parts[PartName])
        
            adder = adder + Dy3[i] + Dy1[i] + 1.7e-3;
            mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(i+1), ), 
                vector=(0.0, -1*(adder), 0.0))

    # adder = 0;
    # for i in range(N_part):
    # # for i in [0,1,2,3]:
        
        # PartName = 'Part'+str(i+1)
        # mdb.models[ModelName].rootAssembly.Instance(dependent=ON, name=PartName+'-'+str(1), 
            # part=mdb.models[ModelName].parts[PartName])
            
        
        # if i % 2 != 0:
            # adder = adder + h_all[i];
            # mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(1), ), 
                # vector=(0.0, -1*(adder+(1.6e-3*i)), 0.0))
        # else:
            # mdb.models[ModelName].rootAssembly.translate(instanceList=(PartName+'-'+str(1), ), 
                # vector=(0.0, -1*(adder+(1.6e-3*i)), 0.0))   
            # adder = adder + h_all[i];
            
    #DynamicImplicit
    mdb.models[ModelName].ImplicitDynamicsStep(alpha=DEFAULT, application=
        QUASI_STATIC, initialConditions=ON, initialInc=0.05, maxInc=0.05, maxNumInc=100000
        , name='Step-1', nlgeom=ON, nohaf=OFF, previous='Initial')
                
    mdb.models[ModelName].ContactProperty('FrictionCoeff')
    mdb.models[ModelName].interactionProperties['FrictionCoeff'].TangentialBehavior(
        dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
        formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
        pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
        table=((friction, ), ), temperatureDependency=OFF)
        
        
    for i in range(N_part-1):
        PartName = 'Part'+str(i+1)+'-'+str(1)
        if i%2 == 0:
            BlockName = 'Block1'+'-'+str(i+1)
        else:
            BlockName = 'Block2'+'-'+str(i+1)
            
        mdb.models[ModelName].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
            clearanceRegion=None, createStepName='Step-1', datumAxis=None, 
            initialClearance=OMIT, interactionProperty='FrictionCoeff', master=
            mdb.models[ModelName].rootAssembly.instances[PartName].surfaces['Bottom'], 
            name='FrictionBlock'+str(i)+str(1), slave=
            mdb.models[ModelName].rootAssembly.instances[BlockName].surfaces['Top'], 
            sliding=FINITE, thickness=ON)   
            
        PartName = 'Part'+str(i+1+1)+'-'+str(1)

        mdb.models[ModelName].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
            clearanceRegion=None, createStepName='Step-1', datumAxis=None, 
            initialClearance=OMIT, interactionProperty='FrictionCoeff', master=
            mdb.models[ModelName].rootAssembly.instances[PartName].surfaces['Top'], 
            name='FrictionBlock'+str(i)+str(2), slave=
            mdb.models[ModelName].rootAssembly.instances[BlockName].surfaces['Bottom'], 
            sliding=FINITE, thickness=ON)            
    

    # for i in range (N_part):
        # PartName = 'Part'+str(i+1)
        
        # if i % 2 == 0:
            # try:
                # mdb.models[ModelName].Coupling(controlPoint=
                    # mdb.models[ModelName].rootAssembly.instances[PartName+'-'+str(1)].sets['P3'], 
                    # couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
                    # name='Constraint'+'-'+str(i+1), surface=
                    # mdb.models[ModelName].rootAssembly.instances['Part'+str(i+1+1)+'-'+str(1)].sets['P4'], u1=ON, 
                    # u2=ON, ur3=ON)
            # except:
                # print('End')
        # if i % 2 != 0:
            # try:
                # mdb.models[ModelName].Coupling(controlPoint=
                    # mdb.models[ModelName].rootAssembly.instances[PartName+'-'+str(1)].sets['P2'], 
                    # couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
                    # name='Constraint'+'-'+str(i+1), surface=
                    # mdb.models[ModelName].rootAssembly.instances['Part'+str(i+1+1)+'-'+str(1)].sets['P1'], u1=ON, 
                    # u2=ON, ur3=ON)
            # except:
                # print('End')
                
                
    #BC
    mdb.models[ModelName].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
        region=mdb.models[ModelName].rootAssembly.instances['Part'+str(N_part)+'-'+str(1)].sets['P2'], 
        u1=UNSET, u2=SET, ur3=UNSET)
        
    #Gravity
    mdb.models[ModelName].Gravity(comp2=-10, createStepName='Step-1', 
        distributionType=UNIFORM, field='', name='Gravity')
        
    #Disp
    mdb.models[ModelName].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=
        mdb.models[ModelName].rootAssembly.instances['Part'+str(1)+'-'+str(1)].sets['P1'], u1=
        UNSET, u2= -2*sum(h_all), ur3=UNSET)
        

    ##Point Mass
    mdb.models[ModelName].rootAssembly.engineeringFeatures.PointMassInertia(alpha=
        0.0, composite=0.0, mass=10.8, name='Inertia-1', region=
        mdb.models[ModelName].rootAssembly.instances['Part'+str(1)+'-'+str(1)].sets['P1'])
        
        
    #Output
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(frequency=1)
    mdb.models[ModelName].historyOutputRequests['H-Output-1'].setValues(frequency=1)
    mdb.models[ModelName].HistoryOutputRequest(createStepName='Step-1', name=
        'F-Output-2',frequency= 1, rebar=EXCLUDE, region=
        mdb.models[ModelName].rootAssembly.allInstances['Part'+str(1)+'-'+str(1)].sets['P1'], 
        sectionPoints=DEFAULT, variables=('U2', 'RF2'))
    mdb.models[ModelName].HistoryOutputRequest(createStepName='Step-1', name=
        'F-Output-3',frequency= 1, rebar=EXCLUDE, region=
        mdb.models[ModelName].rootAssembly.allInstances['Part'+str(N_part)+'-'+str(1)].sets['P2'], 
        sectionPoints=DEFAULT, variables=('RF2', ))
    mdb.models[ModelName].FieldOutputRequest(createStepName='Step-1', name=
        'F-Output-2', rebar=EXCLUDE, region=
        mdb.models[ModelName].rootAssembly.allInstances['Part'+str(1)+'-'+str(1)].sets['P1'], 
        sectionPoints=DEFAULT, variables=('U', 'RF'))
    mdb.models[ModelName].FieldOutputRequest(createStepName='Step-1', name=
        'F-Output-3', rebar=EXCLUDE, region=
        mdb.models[ModelName].rootAssembly.allInstances['Part'+str(N_part)+'-'+str(1)].sets['P2'], 
        sectionPoints=DEFAULT, variables=('RF', ))

        
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model=ModelName, modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=JobName, nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        
try:
    del mdb.models['Model-1']
except:
    print('No Model 1')
    

def RunModel():
    RunYN = getInput('Run? (Y/N)')

    if RunYN == 'Y' or RunYN == 'y':  

        for k in range(N_models):

            ModelName = 'Stack_'+str(k+1)   
            mdb.jobs[ModelName].submit(consistencyChecking=OFF)      

        for k in range(N_models):

            ModelName = 'Stack_'+str(k+1)
            mdb.jobs[ModelName].waitForCompletion() 
            
def RunModel2():
    RunYN = getInput('Run? (Y/N)')

    if RunYN == 'Y' or RunYN == 'y':  

        for k in range(N_models):

            ModelName = 'Stack_'+str(k+1)   
            mdb.jobs[ModelName].submit(consistencyChecking=OFF)      
                 

##The Following Section Automates Output Extraction
##Edit the JobPath Variable To Run

def OutputModel():
    global h_all_models
    global JobPath
    
            
    OutYN = getInput('extractData? (Y/N)')
        
    if OutYN == 'Y' or OutYN == 'y':  

        for k in range(N_models):
        
            ModelName = 'Stack_'+str(k+1) 
            
            N_part = len(h_all_models[k])
            
            o1 = session.openOdb(name = JobPath + str(ModelName)+ '.odb')
            odb = session.odbs[JobPath + str(ModelName)+ '.odb']
            
            session.viewports['Viewport: 1'].setValues(displayedObject=o1)
            session.viewports['Viewport: 1'].makeCurrent()
             
            SE = session.XYDataFromHistory(name='SE-'+str(k+1), odb=odb, 
                outputVariableName='Strain energy: ALLSE for Whole Model', steps=(
                'Step-1', ), __linkedVpName__='Viewport: 1')
              
            U = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=("PART1-1.P1", ))
            
            RFTop = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=("PART1-1.P1", ))
                
            RFBottom = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=("PART"+str(N_part)+"-1.P2", ))
                
            SEU = combine(-1*U[0], SE)
            tmpName = SEU.name
            session.xyDataObjects.changeKey(tmpName, 'SEU-'+str(k+1))
            
            # SEUDiff = differentiate ( SEU,  )
            # tmpName = SEUDiff.name
            # session.xyDataObjects.changeKey(tmpName, 'SEUDiff-'+str(k+1))
            
            RFTopU = combine(-1*U[0], -1*RFTop[0])
            tmpName = RFTopU.name
            session.xyDataObjects.changeKey(tmpName, 'RFTopU-'+str(k+1))
            
            RFBottomU = combine(-1*U[0], RFBottom[0])
            tmpName = RFBottomU.name
            session.xyDataObjects.changeKey(tmpName, 'RFBottomU-'+str(k+1))