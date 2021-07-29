



import numpy as np
import matplotlib.pyplot as plt


from mesh import Mesh


def main():

    NX         =   200
    NY         =   200
    xsize      =   1.
    ysize      =   1.
    t          =   0.
    tEnd       =   2.
    tOut       =   0.02
    dt         =   0.01

    plot       =   True

    mesh = Mesh(NX, NY, -1,1, -1,1)
 
    print(mesh.xlin)

    z = np.sin(mesh.x**2 + mesh.y**2) / (mesh.x**2 + mesh.y**2)
    #while t < tEnd:
    #    t = t+dt
    #    print(t)
    h = plt.contourf(mesh.x,mesh.y,z)
    plt.show()


if __name__== "__main__":
    main()
