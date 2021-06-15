from dolfin import *
import sys
sys.path.append("../../../../FMECHS/")
import fenicsmechanics as fm
import numpy as np

Nr=100
Nt=100
Nz=1
Ly=5
Lx=2*Ly
lxa=1.0*Lx/20
Lz=1.00
meshname=BoxMesh(Point(0.0,0.0,0.0),Point(Lx,Ly,Lz),Nr,Nt,Nz)#
#Lzn=0.001
Lzn=0.005

x=meshname.coordinates()[:,0]
y=meshname.coordinates()[:,1]
z=meshname.coordinates()[:,2]
def transmesh(x,y,z,Lzn):#
    return [x,y,Lzn*z]
qn=transmesh(x,y,z,Lzn)
qnx=np.array(qn).transpose()
meshname.coordinates()[:]=qnx
tol=1E-14
##########

#xl =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
#xla =  CompiledSubDomain("(near(x[1], side) && on_boundary) && x[0]<lxa", side = 0.0, lxa=lxa)
#xlb =  CompiledSubDomain("(near(x[1], side) && on_boundary) && x[0]<lxa", side = Ly, lxa=lxa)
xlc =  CompiledSubDomain("(near(x[2], side) && on_boundary) && x[0]<lxa", side = 0.0, lxa=lxa)
xld =  CompiledSubDomain("(near(x[2], side) && on_boundary) && x[0]<lxa", side = Lzn, lxa=lxa)

#xr = CompiledSubDomain("near(x[0], side) && on_boundary", side = Lx)
#xra = CompiledSubDomain("(near(x[1], side) && on_boundary) && x[0]>lxb", side = 0.0, lxb=Lx-lxa)
#xrb = CompiledSubDomain("(near(x[1], side) && on_boundary) && x[0]>lxb", side = Ly, lxb=Lx-lxa)
xrc = CompiledSubDomain("(near(x[2], side) && on_boundary) && x[0]>lxb", side = 0.0, lxb=Lx-lxa)
xrd = CompiledSubDomain("(near(x[2], side) && on_boundary) && x[0]>lxb", side = Lzn, lxb=Lx-lxa)

boundary_parts = MeshFunction("size_t", meshname, meshname.topology().dim() - 1)
boundary_parts.set_all(0)

#xl.mark(boundary_parts, 1)
#xla.mark(boundary_parts,2)
#xlb.mark(boundary_parts,3)
xlc.mark(boundary_parts,1)
xld.mark(boundary_parts,2)

#xr.mark(boundary_parts, 4)
#xra.mark(boundary_parts,5)
#xrb.mark(boundary_parts,6)
xrc.mark(boundary_parts,3)
xrd.mark(boundary_parts,4)

bds=File("bmarksmarks.pvd")
bds << boundary_parts

###########################################################################
dlx=0.0

cr = Expression(("x0","kr*x[1]","z0"),x0 = 0.0, kr = 0.0, z0=0.0000,degree=1)
cra = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0000,degree=1)
crb = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0,degree=1)

cl = Expression(("x0","kl*x[1]","z0"),x0 = 0.0, kl = 0.0, z0=0.0,degree=1)
cla = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0,degree=1)
clb = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0,degree=1)
#cu = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0,degree=1)
#cd = Expression(("x0","y0","z0"),x0 = 0.0, y0 = 0.0, z0=0.0,degree=1)

material = {
    'type': 'elastic',
    'const_eqn': 'neo_hookean',
    'incompressible': True,
    'kappa': 10e9, # Pa
    'mu': 1.5e6 # Pa
}

#mesh_file, boundaries = fm.get_mesh_file_names("unit_domain", ret_facets=True,
#                                               refinements=[20, 20])
mesh = {
    'mesh_file': meshname,
    'boundaries': boundary_parts
}

formulation = {
    'element': 'p2-p1',
    'domain': 'lagrangian',
    'bcs': {
        'dirichlet': {
            #'displacement': [cr,cl,cu,cd],
            'displacement': [cl,cl,cr,cr],
            'regions': [1,2,3,4],
            }
    ##    'neumann': {
    #        'values': [[1e6, 0.]],
    #        'regions': [2],
    ##        'types': ['piola']
    #    }
    }
}

config = {
    'material': material,
    'mesh': mesh,
    'formulation': formulation
}


filen = "2ddisplacement.pvd"
problem = fm.SolidMechanicsProblem(config)
#c1.scale=0.01
solver = fm.SolidMechanicsSolver(problem, fname_disp=filen)
solver.set_parameters(linear_solver="mumps")
solver.set_parameters(newton_maxIters=500)
solver.set_parameters(newton_abstol=1e-7)
solver.set_parameters(newton_reltol=1e-6)
#solver.full_solve()
#print(mu)

d1=0.000
d2=0.000
#DX=0.01
#DY=0.01
ky=0.0
kx=0.0
for j in range(4000):#


    if j<50 or j> 8000:
        DX=0.01
    else:
        DX=0.001

    if j<5 or j>800:
        DY=0.001
    else:
        DY=1e-9

    d2+=DY
    print(d1)
    d1+=DX
    cr.x0=d1
    ky=d2/Ly
    cl.kl=ky
    cr.kr=ky

    solver.full_solve()
