
import numpy as np



def get_gradient(f, dx, dy):
    """
    Calculate the gradients of a field
    f        is a matrix of the field
    dx       is the cell size
    f_dx     is a matrix of derivative of f in the x-direction
    f_dy     is a matrix of derivative of f in the y-direction
    """
    # directions for np.roll() 
    R = -1   # right
    L = 1    # left
    
    f_dx = ( np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (2*dx)
    f_dy = ( np.roll(f,R,axis=1) - np.roll(f,L,axis=1) ) / (2*dy)
    
    return f_dx, f_dy



def interface_flux(q, u, dt, dx):
    """
    """
    theta = np.sign(u)
    slope = slope(q, u)
    limiter = limiter(slope)
    flux = 0.5 * u * ( (1.0+theta) * q[:-1] ) + (1.0-theta) * q[1:] +
           0.5 * np.absolute(u) * (1.0 - np.abs(v*dt/dx))*limiter*(q[1:]-q[:-1])



def update(q, u, dt, dx, ng=1):
    pass
    


def slope(q, u):
    """compute slope at cell interfaces.
    """
    r = np.zeros_like(q, dtype=float)
    diff = q[1:] - q[0:-1]
    
    diff_l = diff[:-2]
    diff_r = diff[2:]
    diff_c = diff[1:-1]

    r[2:-1] = np.where(diff_c > 0.0, 0.0, np.where( u[2:-1] > 0.0, 
        diff_l / diff_c, diff_r / diff_c))
    #non_zero = (abs(diff_c) > 0.0)
    #r[2:-1] = np.where(u[2:-1] >= 0.0, diff_l / diff_c, diff_r / diff_c)
    #r[2:-1] = np.where(np.logical_not(non_zero), 0.0, r[2:-1])

    return r



def limiter(slope, name='donor-cell'):
    if name=='donor-cell':
        return 0.0
    elif name=='Lax-Wendroff':
        return 1.0
    else:
        return 0.0

def centered(f, dx, xs=1):
    return ( f[xs+1:] - f[:-xs-1] ) / (2*dx)

def upwind(f, dx, xs=1):
    return ( f[xs:-xs] - f[:-xs-1] ) / dx

def advect_1d(q, dt, vx, dx, xs=1):
    q_dx = centered(q, dx, xs)
    return q[xs:-xs] - dt * vx * q_dx 


def boundary_conditions(q, xs=1):
    q[0] = q[1]
    q[-1] = q[-2]
    return q


def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    return u0, u


if __name__ == "__main__":
    print('hydro test')
    q = np.array((1,5,4,2,1,9),dtype=np.float)
    u = np.array((0,0.1,0.5,1,0,0), dtype=np.float)
    r = slope(q,u)
    print(q)
    print(u)
    print(r)

#def interface_flux(q1, q2, q3, q4, vface, dt):
#    theta = sign(1.0,v_face)
# 
#    if(abs(q3-q2).gt.0.0) then
#         if(v_face.ge.0.0) then
#            r = (q2-q1)/(q3-q2)
#         else 
#            r = (q4-q3)/(q3-q2)
#         endif
#      else
#         r = 0.0
#      endif
#      !
#      limiter = phi(fl,r)
#      !
#      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
#             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*limiter*(q3-q2)
#             !
#
