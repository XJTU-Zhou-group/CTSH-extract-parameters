# -*- coding'': mbcs -*-
#Initialize
#From Abaqus:
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
#From Python:
import numpy as np
import math
import datetime
import xyPlot


#############################################################################################################################
# Created by Ning An
# 2019/12
# http://www.anning.me/create-virtual-nodes
# Created in Abaqus Version 2017
#----------------------------------------------------
# Function for creating virtual nodes
# mdb: model database
# NameModel: A string with the name of your model
# NameRef: A string with the name of a virtual node.
# Coord: A vector indicates the coordinates of the virtual node.
# Example: VirtualNodes(mdb, 'Model-1', 'Ref-0', [0.0, 0.0, 0.0])
#############################################################################################################################
def VirtualNodes(mdb, NameModel, NameRef, Coord):
    from part import THREE_D, DEFORMABLE_BODY
    #Create reference parts and assemble
    mdb.models[NameModel].Part(dimensionality=THREE_D, name=NameRef, type=
        DEFORMABLE_BODY)
    mdb.models[NameModel].parts[NameRef].ReferencePoint(point=(Coord[0], Coord[1], Coord[2]))
    mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef, 
        part=mdb.models[NameModel].parts[NameRef])
    mdb.models[NameModel].rootAssembly.Set(name=NameRef, referencePoints=(
        mdb.models[NameModel].rootAssembly.instances[NameRef].referencePoints[1],))


index = np.genfromtxt('index.csv', dtype=float, delimiter=',')
index = int(index)

index = index + 1
index = int(index)
with open('index.csv', 'w') as file_object:
    file_object.write(str(index))

E1_table = np.genfromtxt('E1.csv', dtype=float, delimiter=',')
E2_table = np.genfromtxt('E2.csv', dtype=float, delimiter=',')
V12_table = np.genfromtxt('V12.csv', dtype=float, delimiter=',')
G12_table = np.genfromtxt('G12.csv', dtype=float, delimiter=',')
G23_table = np.genfromtxt('G23.csv', dtype=float, delimiter=',')

E1=E1_table[index-1]
E2=E2_table[index-1]
V12=V12_table[index-1]
G12=G12_table[index-1]
G23=G23_table[index-1]

V23=E2/(2*G23)-1




mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(50.0, 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShellExtrude(depth=600.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
Plane_YZ_50 = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=50.0, 
    principalPlane=YZPLANE)
Plane_XY_0 = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XYPLANE)
Axis = mdb.models['Model-1'].parts['Part-1'].DatumAxisByTwoPlane(plane1=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_XY_0.id], plane2=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_YZ_50.id])
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=33.98, name='__profile__', 
    sheetSize=1359.55, transform=
    mdb.models['Model-1'].parts['Part-1'].MakeSketchTransform(
    sketchPlane=mdb.models['Model-1'].parts['Part-1'].datums[Plane_YZ_50.id], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=mdb.models['Model-1'].parts['Part-1'].datums[Axis.id], 
    sketchOrientation=RIGHT, origin=(50.0, 0.0, 300.0)))
mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])

mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(-210.0, 0.0)
    , direction=COUNTERCLOCKWISE, point1=(-210.0, 40.0), point2=(-210.0, 
    -40.0))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(210.0, 0.0), 
    direction=CLOCKWISE, point1=(210.0, 40.0), point2=(210.0, -40.0))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-210.0, 40.0), point2=(
    210.0, 40.0))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(210.0, -40.0), point2=(
    -210.0, -40.0))

mdb.models['Model-1'].parts['Part-1'].CutExtrude(flipExtrudeDirection=OFF, 
    sketch=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
    RIGHT, sketchPlane=mdb.models['Model-1'].parts['Part-1'].datums[2], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
    mdb.models['Model-1'].parts['Part-1'].datums[4])
del mdb.models['Model-1'].sketches['__profile__']


Plane_XY_40 = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=40.0, 
    principalPlane=XYPLANE)
Plane_XY_560 = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=560.0, 
    principalPlane=XYPLANE)
mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_XY_40.id], faces=
    mdb.models['Model-1'].parts['Part-1'].faces)
mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_XY_560.id], faces=
    mdb.models['Model-1'].parts['Part-1'].faces)

Plane_YZ = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=YZPLANE)
Plane_XZ = mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_YZ.id], faces=
    mdb.models['Model-1'].parts['Part-1'].faces)
mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[Plane_XZ.id], faces=
    mdb.models['Model-1'].parts['Part-1'].faces)

mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Density(table=((1.58076e-09, ), 
    ))
mdb.models['Model-1'].materials['Material-1'].Elastic(table=(
    (E1, E2, E2, V12, V12, V23, G12, G12, G23), ), type=
    ENGINEERING_CONSTANTS, temperatureDependency=OFF)


mdb.models['Model-1'].CompositeShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, layup=(SectionLayer(thickness=0.187, 
    orientAngle=45.0, material='Material-1', plyName='1'), SectionLayer(
    thickness=0.187, orientAngle=-45.0, material='Material-1', plyName='2')), 
    name='Section-1', poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=
    False, temperature=GRADIENT, thicknessModulus=None, thicknessType=UNIFORM, 
    useDensity=OFF)
#
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts['Part-1'].faces), 
    sectionName='Section-1', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(
    additionalRotationType=ROTATION_NONE, axis=AXIS_2, fieldName='', localCsys=
    None, orientationType=GLOBAL, region=Region(
    faces=mdb.models['Model-1'].parts['Part-1'].faces))

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].ExplicitDynamicsStep(improvedDtMethod=ON, name='Step-1', 
    previous='Initial',timePeriod=0.4)
mdb.models['Model-1'].FieldOutputRequest(createStepName='Step-1', name=
    'F-Output-1', timeInterval=0.005, variables=('S', 'LE', 'U', 'UR', 'V', 'RF', 'RM'), filter=ANTIALIASING)
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name=
    'H-Output-1', timeInterval=0.005, variables=('ALLAE', 'ALLIE', 'ALLKE', 'ALLSE', 'ALLVD', 'ALLWK', 'ETOTAL'), filter=ANTIALIASING)

mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
mdb.models['Model-1'].ContactExp(createStepName='Initial', name='Int-1')
mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')



VirtualNodes(mdb, 'Model-1', 'Ref-1', [0.0, 0.0, 0.0])
VirtualNodes(mdb, 'Model-1', 'Ref-2', [0.0, 0.0, 600.0])

mdb.models['Model-1'].rootAssembly.Surface(name='s_Surf-1', side2Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(
    ((-35.355339,-35.355339,0.),), 
    ((35.355339,-35.355339,0. ),),
    ((35.355339,35.355339,0. ),),
    ((-35.355339,35.355339,0. ),),))
mdb.models['Model-1'].Coupling(controlPoint=
    mdb.models['Model-1'].rootAssembly.sets['Ref-1'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1', 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['s_Surf-1'], u1=ON, u2=
    ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

mdb.models['Model-1'].rootAssembly.Surface(name='s_Surf-2', side2Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(
    ((-35.355339,-35.355339,600.),), 
    ((35.355339,-35.355339,600. ),),
    ((35.355339,35.355339,600. ),),
    ((-35.355339,35.355339,600. ),),))
mdb.models['Model-1'].Coupling(controlPoint=
    mdb.models['Model-1'].rootAssembly.sets['Ref-2'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['s_Surf-2'], u1=ON, u2=
    ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)


mdb.models['Model-1'].Gravity(comp3=-9870.0, createStepName='Step-1', 
    distributionType=UNIFORM, field='', name='Load-1')
mdb.models['Model-1'].SmoothStepAmplitude(data=((0.0, 0.0), (3.0,1.0)), name=
    'Amp-1', timeSpan=STEP)
mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Step-1'
    , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'BC-1', region=mdb.models['Model-1'].rootAssembly.sets['Ref-1'], u1=0.0, u2=0.0, 
    u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Step-1'
    , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'BC-2', region=mdb.models['Model-1'].rootAssembly.sets['Ref-2'], u1=0.0, u2=UNSET, 
    u3=UNSET, ur1=-1*pi, ur2=0.0, ur3=0.0)

mdb.models['Model-1'].rootAssembly.Surface(name='Surf-All', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces)
mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
    distributionType=VISCOUS, field='', magnitude=1.5e-08, name='Load-2', 
    refPoint=mdb.models['Model-1'].rootAssembly.sets['Ref-1'], region=
    mdb.models['Model-1'].rootAssembly.surfaces['Surf-All'])

mdb.models['Model-1'].rootAssembly.Set(faces=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces, name='Set-All')


mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=6.0)
mdb.models['Model-1'].parts['Part-1'].setMeshControls(elemShape=QUAD, algorithm=MEDIAL_AXIS, regions=
    mdb.models['Model-1'].parts['Part-1'].faces)
mdb.models['Model-1'].parts['Part-1'].generateMesh()

mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
    hourglassControl=ENHANCED), ElemType(elemCode=S3R, elemLibrary=EXPLICIT)), 
    regions=(mdb.models['Model-1'].parts['Part-1'].faces, ))


mdb.models['Model-1'].rootAssembly.regenerate()

jobName = "Job-1"
mdb.Job(explicitPrecision=DOUBLE_PLUS_PACK, model='Model-1', name=jobName, 
    nodalOutputPrecision=FULL, numCpus=5, numDomains=5)
mdb.saveAs(pathName=jobName)
mdb.jobs[jobName].submit()
mdb.jobs[jobName].waitForCompletion()


o2 = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
odb = session.odbs['Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'RM_ANTIALIASING', NODAL, ((COMPONENT, 'RM1'), )),('UR_ANTIALIASING', NODAL, ((COMPONENT, 'UR1'), 
    )), ), nodeSets=("REF-1", ))
session.xyDataObjects.changeKey(fromName='_RM_ANTIALIASING:RM1 PI: REF-1 N: 1', 
    toName='Moment')


Moment = session.xyDataObjects['Moment']
M_data = np.array(Moment)[:, 1:2]


with open('Moment.csv', 'a') as file_object:
    file_object.write(str(index) + ',')
    for i in range(np.size(M_data, 0)):
        if i != np.size(M_data, 0) - 1:
            file_object.write(str(M_data[i, 0]) + ',')
        else:
            file_object.write(str(M_data[i, 0]))
    file_object.write('\n')


