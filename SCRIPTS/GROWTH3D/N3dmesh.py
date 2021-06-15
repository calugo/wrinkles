from dolfin import *
import numpy as np

Nx=40
Ny=20
Nz=41
Ly=2.0
Lx=1.0#*Ly
Lz=0.1
Ho=Lz/Nz
#Ho=0.5
DH=Ho
DHZ=Ho
print(Lx,Ly,Lz,Ho)

#mesh=BoxMesh(Point(-Lx,-Ly,0.0),Point(Lx,Ly,Lz),Nr,Nt,Nz)#
#x=mesh.coordinates()[:,0]
#y=mesh.coordinates()[:,1]
#z=mesh.coordinates()[:,2]
#def transmesh(x,y,z,n):#
#    return [x,y,z**(1/n)]
#qn=transmesh(x,y,z,5)
#qnx=np.array(qn).transpose()
#mesh.coordinates()[:]=qnx

###############################################################

mesh = UnitCubeMesh.create(Nx, Ny, Nz, CellType.Type.hexahedron)
# mesh size is smaller near x=y=0
mesh.coordinates()[:, 0] = mesh.coordinates()[:, 0]*Lx
mesh.coordinates()[:, 1] = mesh.coordinates()[:, 1]*Ly
mesh.coordinates()[:, 2] = mesh.coordinates()[:, 2]*Lz



################################################################
#tol=1E-14
#materials=MeshFunction("size_t",mesh,mesh.topology().dim())
## Mark boundary subdomians
#film=CompiledSubDomain("x[2]>=Rcut",Rcut=Lz-Ho) #,size=1.0)

#wall=CompiledSubDomain("(x[2]<=dxwz && x[0]<=Ax && x[0]>=-Ax && x[1]<=Ax && x[1]>=-Ax) ||\
#                        ((x[2]<=Rcut) && (\
#                        ((x[0]>=Ax-dxw)&&(x[0]<=Ax)) && ((x[1]<=Ax)&&(x[1]>=-Ax)) ||\
#                        ((x[0]>=-Ax)&&(x[0]<=-Ax+dxw)) && ((x[1]<=Ax)&&(x[1]>=-Ax)) ||\
#                        ((x[1]>=-Ax)&&(x[1]<=-Ax+dxw)) && ((x[0]<=Ax)&&(x[0]>=-Ax)) ||\
#                        ((x[1] >= (Ax-dxw)) && (x[1]<=Ax)) && ((x[0]<=Ax) && (x[0]>=-Ax))\
#                        ))",Rcut=Lz-Ho,dxw=DH,Ax=A,dxwz=DHZ)
#wall=CompiledSubDomain("(x[2]<=dxwz && x[0]<=AX && x[0]>=-AX && x[1]<=AY && x[1]>=-AY) ||\
#                        ((x[2]<=Rcut) && (\
#                        ((x[0]>=Ax-dxw)&&(x[0]<=Ax)) && ((x[1]<=Ay)&&(x[1]>=-Ay)) ||\
#                        ((x[0]>=-Ax)&&(x[0]<=-Ax+dxw)) && ((x[1]<=Ay)&&(x[1]>=-Ay)) ||\
#                        ((x[1]>=-Ay)&&(x[1]<=-Ay+dxw)) && ((x[0]<=Ax)&&(x[0]>=-Ax)) ||\
#                        ((x[1] >= (Ay-dxw)) && (x[1]<=Ay)) && ((x[0]<=Ax) && (x[0]>=-Ax))\
#                        ))",Rcut=Lz-Ho,dxw=DH,Ax=Ax,Ay=Ay,dxwz=DHZ,AX=Lx,AY=Ly)

#wall=CompiledSubDomain("(x[2]<=dxwz )",dxwz=DH)




#materials.set_all(0)
#film.mark(materials,1)
#wall.mark(materials,2)
#bmats=File("matx.pvd")
#bmats << materials

f = HDF5File(mesh.mpi_comm(), "meshq.hdf5", 'w')
f.write(mesh, "mesh")

#f = HDF5File(mesh.mpi_comm(), "meshq.hdf5", 'a')
#f.write(materials, "materials")
## Write out both the `Mesh` and the `MeshFunction`:
#File("mesh.xml") << mesh
##File("mf.xml") << materials
