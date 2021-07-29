

import numpy as np


class Mesh():


    def __init__(self, nx, ny, x1, x2, y1, y2, ng=1):
        # cell centers
        self.nx = nx
        self.ny = ny
        self.ng = ng
        xsize = x2 - x1
        ysize = y2 - y1
        self.dx = xsize / nx
        self.dy = ysize / ny
        # linear centered coordinates
        self.xlin = np.linspace(x1+(0.5-ng)*self.dx, x2+(ng-0.5)*self.dx, nx+2*ng)
        self.ylin = np.linspace(y1+(0.5-ng)*self.dy, y2+(ng-0.5)*self.dy, ny+2*ng)
        # linear left interface coordinates
        self.xl = np.linspace(x1+(-ng)*self.dx, x2+(ng-1)*self.dx, nx+2*ng)
        self.yl = np.linspace(y1+(-ng)*self.dy, y2+(ng-1)*self.dy, ny+2*ng)
        # linear right interface coordinates
        self.xr = np.linspace(x1+(-ng+1)*self.dx, x2+(ng)*self.dx, nx+2*ng)
        self.yr = np.linspace(y1+(-ng+1)*self.dy, y2+(ng)*self.dy, ny+2*ng)

        self.xc, self.yc = np.meshgrid( self.xlin, self.ylin )

    def scratch_array(self):
        return np.zeros((self.nx+2*self.ng, self.ny+2*self.ng), dtype=np.float64)



if __name__ == "__main__":
    mesh = Mesh(20,20,-1,1,-1,1,ng=1)
    print(mesh.scratch_array())



