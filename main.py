


import numpy as np
import matplotlib.pyplot as plt


from mesh import Mesh
import hydro

def plot_realtime(x, rho):
    axes = plt.gca()
    axes.set_ylim([0.,1.2])

    plt.plot(x, rho)
    plt.pause(0.0001)
    plt.clf()






def main():

    NX         =   100
    NY         =   100
    xsize      =   100.
    ysize      =   100.
    t          =   0.
    tEnd       =   2.
    tOut       =   0.02
    dt         =   0.01

    mesh = Mesh(NX, NY, -1,1, -1,1)
 
    # initial conditions
    rho = np.where(mesh.xlin < 0, 1., 0.)

    z = np.sin(mesh.xc**2 + mesh.yc**2) / (mesh.xc**2 + mesh.yc**2)

    # prep figure
    fig = plt.figure(figsize=(4,4), dpi=80)

    while t < tEnd:
        plot = False
        u = 0.1
        rho = hydro.boundary_conditions(rho)
        rho[1:-1] = hydro.advect_1d(rho, dt, u, mesh.dx)
        #if (t%tOut < 1.e-8):
        #    plot_realtime(mesh.xlin, rho)
        t = t+dt
        print(t)
    #print(z)
    #h = plt.contourf(mesh.x,mesh.y,z)
    plt.plot(mesh.xlin, rho)
    plt.show()

if __name__== "__main__":
    main()
