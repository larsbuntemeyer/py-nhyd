!
! Hydro module
!
module Hydro
   !
   implicit none
   !
   public  :: Hydro_init, Hydro_solve
   !
   private :: hydro1D, fill_guardcells, minmod, &
              interface_flux, cfl_timestep
   !
contains
   !
   subroutine advect

   end subroutine advect
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_init
      !
      implicit none
      !
      write(*,*) '----- Hydro_init ------------------'
      write(*,*) 'Hydro is not doing anything yet,   '
      write(*,*) 'because i have to admit, i was lazy'
      write(*,*) '----- Hydro_init done -------------'
      !
   end subroutine Hydro_init
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_solve(dt)
      !
      use Grid, only: ndim,ib,ie
      use Database
      use RuntimeParameters, only: bc
      !
      implicit none
      !
      real, intent(inout) :: dt
      integer :: i,j,k
      real,dimension(nvar,nx+1) :: F_l,F_r
      !
      if(ndim>1) then
         write(*,*) '------------------------'
         write(*,*) 'ERROR in Hydro_solve:'
         write(*,*) 'only 1D possilbe!'
         write(*,*) '------------------------'
         stop
      endif
      !
      !call fill_guardcells
      !
      do k=kbg,keg
        do j=jbg,jeg
           call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),ibg,ieg,ieg,bc)
        enddo  
      enddo  
      !do j=jbg,jeg
      !  do i=ibg,ieg
      !     call fill_guardcells_1D(dens(i,j,:),pres(i,j,:),eint(i,j,:),u(i,j,:),v(i,j,:),w(i,j,:),kbg,keg,keg,2)
      !  enddo  
      !enddo  
      call hydro1D(dt)
      !do k=kb,ke
      !  do j=jb,je
      !    !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
      !    !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
      !    call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),ibg,ieg,ieg,2)
      !    call sw
      !  enddo  
      !enddo  
      !
      !
   end subroutine Hydro_solve
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine hydro1D(dt)
      !
      use Grid, only: ib,ie,kbg,keg,jbg,jeg,ibg,ieg,jb,je,kb,ke, &
                      xcCoord,ycCoord,zcCoord,       &
                      xlCoord,xrCoord,zlCoord,zrCoord, &
                      ylCoord,yrCoord,k2d,k3d
      use Database, only: dens,pres,u,v,w,eint,nvar,uf,vf,wf,ener
      use Eos, only: Eos_gamma
      !
      implicit none
      !
      real, intent(in) :: dt
      !
      ! the current fluxes:
      ! the flux is defined at the cell interfaces, so we
      ! for an x-grid with nx cells, we have nx+2*nguard+1
      ! cell interfaces 
      !
      real, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: xflux,yflux,zflux
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: q
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: dq
      !
      real :: rho,ei,ekin,vtot2,dp,dx,dy,dz,rhoinv
      real :: fx,fy,fz
      !
      integer :: i,j,k,ivar
      !
      ! Compute conserved variables
      !
      do k=kbg,keg
         do j=jbg,jeg
            do i=ibg,ieg
               !
               vtot2 = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               rho   = dens(i,j,k)
               ekin  = 0.5d0 * vtot2
               !
               q(1,i,j,k) = rho 
               q(2,i,j,k) = rho * u(i,j,k)
               q(3,i,j,k) = rho * v(i,j,k)
               q(4,i,j,k) = rho * w(i,j,k)
               q(5,i,j,k) = rho * (eint(i,j,k) + ekin)
               !
            enddo
         enddo
      enddo 
      !
      ! Interpolate velocities at cell interfaces
      !
      do k=kb,ke+k3d
         do j=jb,je+k2d
            do i=ib,ie+1
               !
               uf(i,j,k) = 0.5d0 * (u(i,j,k) + u(i-1,j,k))
               if(k2d==1) vf(i,j,k) = 0.5d0 * (v(i,j,k) + v(i,j-1,k))
               if(k3d==1) wf(i,j,k) = 0.5d0 * (w(i,j,k) + w(i,j,k-1))
               !
            enddo
         enddo
      enddo 
      !
      ! construct the fluxes at cell interfaces       
      !
      do k=kb,ke+k3d
         do j=jb,je+k2d
            do i=ib,ie+1
               !
               do ivar=1,nvar
                 xflux(ivar,i,j,k) = interface_flux(q(ivar,i-2,j,k),q(ivar,i-1,j,k), &
                                                q(ivar,i,j,k),  q(ivar,i+1,j,k), &
                                                uf(i,j,k),dt)
                 if(k2d==1) then
                 yflux(ivar,i,j,k) = interface_flux(q(ivar,i,j-2,k),q(ivar,i,j-1,k), &
                                                q(ivar,i,j,k),  q(ivar,i,j+1,k), &
                                                vf(i,j,k),dt)
                 endif
                 if(k3d==1) then
                 zflux(ivar,i,j,k) = interface_flux(q(ivar,i,j,k-2),q(ivar,i,j,k-1), &
                                                q(ivar,i,j,k),  q(ivar,i,j,k+1), &
                                                wf(i,j,k),dt)
                 endif
               enddo
               !
            enddo
         enddo
      enddo 
      !
      ! Now advect
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               do ivar=1,nvar
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(xrCoord(i)-xlCoord(i)) * (xflux(ivar,i,j,k) - xflux(ivar,i+1,j,k))
                  if(k2d==1) then
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(yrCoord(i)-ylCoord(i)) * (yflux(ivar,i,j,k) - yflux(ivar,i,j+1,k))
                  endif
                  if(k3d==1) then
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(zrCoord(i)-zlCoord(i)) * (zflux(ivar,i,j,k) - zflux(ivar,i,j,k+1))
                  endif
               enddo
            enddo
         enddo
      enddo 
      !
      ! Update simple variables
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               rho = q(1,i,j,k)
               dens(i,j,k) = rho 
               u(i,j,k)   = q(2,i,j,k) / rho
               v(i,j,k)   = q(3,i,j,k) / rho
               w(i,j,k)   = q(4,i,j,k) / rho
               vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               ekin        = 0.5d0 * vtot2
               ener(i,j,k) = q(5,i,j,k) / rho
               eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo 
      !
      call Eos_gamma 
      !
      ! Now add pressure force
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               rhoinv = 1.0/dens(i,j,k)
               dp = pres(i+1,j,k) - pres(i-1,j,k)
               dx = xcCoord(i+1)  - xcCoord(i-1)
               u(i,j,k)   = u(i,j,k) - dt * dp/dx * rhoinv
               !if(k2d==1) then 
               !  dp = pres(i,j+1,k) - pres(i,j-1,k)
               !  dy = ycCoord(i-1)  - ycCoord(i+1)
               !  v(i,j,k)   = v(i,j,k) - dt/dy * dp * rhoinv
               !endif 
               !if(k3d==1) then 
               !  dp = pres(i,j+1,k) - pres(i,j-1,k)
               !  dz = zcCoord(i-1)  - zcCoord(i+1)
               !  w(i,j,k)   = w(i,j,k) - dt/dz * dp * rhoinv
               !endif 
               ! 
               !vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               !ekin        = 0.5d0 * vtot2
               !ener(i,j,k) = q(5,i,j,k) / rho
               !eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo
      call Eos_gamma 
      !
   end subroutine hydro1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine fill_guardcells_1D(dens,pres,eint,u,v,w,ib,ie,n,bc)
     use Grid, only: nguard
     use RuntimeParameters, only: outflow, periodic
     implicit none
     real ,intent(inout), dimension(n) :: dens,pres,eint,u,v,w
     integer, intent(in) :: ib,ie,n 
     integer, intent(in) :: bc
     integer :: i,j,k 
     select case(bc)
       case(outflow) 
         do i=ib,ib+nguard-1
           dens(i) = dens(ib+nguard)
           pres(i) = pres(ib+nguard)
           eint(i) = eint(ib+nguard)
           u(i)    = u(ib+nguard)
           v(i)    = v(ib+nguard)
           w(i)    = w(ib+nguard)
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(ie-nguard)
           pres(i) = pres(ie-nguard)
           eint(i) = eint(ie-nguard)
           u(i)    = u(ie-nguard)
           v(i)    = v(ie-nguard)
           w(i)    = w(ie-nguard)
         enddo
       case(periodic) 
         do i=ib,ib+nguard-1
           dens(i) = dens(ie-nguard-i+1)
           pres(i) = pres(ie-nguard-i+1)
           eint(i) = eint(ie-nguard-i+1)
           u(i)    = u(ie-nguard-i+1)
           v(i)    = v(ie-nguard-i+1)
           w(i)    = w(ie-nguard-i+1)
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(i-n+nguard+1)
           pres(i) = pres(i-n+nguard+1)
           eint(i) = eint(i-n+nguard+1)
           u(i)    = u(i-n+nguard+1)
           v(i)    = v(i-n+nguard+1)
           w(i)    = w(i-n+nguard+1)
         enddo
       case default
     end select 
   end subroutine fill_guardcells_1D
   !
   subroutine fill_guardcells
      !
      use Grid
      use Database
      !
      implicit none
      !
      integer :: i,j,k,ivar
      !
      ! Fill guardcells x-direction
      !
      do i=ibg,nguard
        do j=jb,je
          do k=kb,ke
            !
            dens(i,j,k) = dens(ib,j,k)
            pres(i,j,k) = pres(ib,j,k)
            eint(i,j,k) = eint(ib,j,k)
            u(i,j,k)    = u(ib,j,k)
            !
          enddo
        enddo
      enddo 
      do i=ie+1,ieg
        do j=jb,je
          do k=kb,ke
            !
            dens(i,j,k) = dens(ie,j,k)
            pres(i,j,k) = pres(ie,j,k)
            eint(i,j,k) = eint(ie,j,k)
            u(i,j,k)    = u(ie,j,k)
            !
          enddo
        enddo
      enddo 
      !
      ! Fill guardcells z-direction
      !
      do k=kbg,nguard
        do j=jb,je
          do i=ib,ie
            !
            dens(i,j,k) = dens(i,j,ke)
            pres(i,j,k) = pres(i,j,ke)
            eint(i,j,k) = eint(i,j,ke)
            w(i,j,k)    = w(i,j,kb)
            !
          enddo
        enddo
      enddo 
      do k=ke+1,keg
        do j=jb,je
          do i=ib,ie
            !
            dens(i,j,k) = dens(i,j,ke)
            pres(i,j,k) = pres(i,j,ke)
            eint(i,j,k) = eint(i,j,ke)
            w(i,j,k)    = w(ie,j,k)
            !
          enddo
        enddo
      enddo 
      !
   end subroutine fill_guardcells
   !
   !----------------------------------------------------------------------------------------------
   !
   real function interface_flux(q1,q2,q3,q4,v_face,dt)
      !
      use Grid
      !
      implicit none
      !
      real :: q1,q2,q3,q4
      real :: v_face,dt
      real :: r,limiter,flux,theta
      !
      theta = sign(1.0,v_face)
      !
      if(abs(q3-q2).gt.0.0) then
         if(v_face.ge.0.0) then
            r = (q2-q1)/(q3-q2)
         else 
            r = (q4-q3)/(q3-q2)
         endif
      else
         r = 0.0
      endif
      !
      limiter = phi(fl,r)
      !
      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*limiter*(q3-q2)
             !
      interface_flux = flux
      !
   end function interface_flux
   !
   !----------------------------------------------------------------------------------------------
   !
   real function phi(fl,r)
   !
   implicit none
   !
   real,         intent(in) :: r
   character(*), intent(in) :: fl
   !
   reaL :: limiter
   !
   select case(fl)
     !
     case('donor-cell')
       limiter = 0.0
     case('Lax-Wendroff')
       limiter = 1.0
     case('Beam-Warming')
       limiter = r
     case('Fromm')
       limiter = 0.5*(1.0+r)
     case('minmod')
       limiter = minmod(1.0,r)
     case('superbee')
       limiter = max(0.0,min(1.0,2.0*r),min(2.0,r))
     case('hyperbee')
       limiter = 1.0!hyperbee(r)
     case('MC')
       limiter = max(0.0,min(0.5*(1.0+r),2.0,2.0*r))
     case('van Leer')
       limiter = (r+abs(r))/(1.0+abs(r))
     case('van Albada 1')
       limiter = (r*r+r)/(r*r+1.0)
     case('van Albada 2')
       limiter = (2.0*r)/(r*r+1.0)
     case default
       limiter = 0.0
   end select
   !
   phi = limiter
   !
   end function phi
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine compute_slope(a,l,n,nvar,r)
     !
     implicit none
     !
     integer, intent(in) :: n,nvar
     real, dimension(nvar,n), intent(in)  :: a,l
     real, dimension(nvar,n), intent(out) :: r
     !
     integer :: i,ivar
     !
     do ivar=1,nvar
       do i=2,n-1
         if(abs(a(ivar,i)) > 0.0) then
            if(l(ivar,i) >= 0.0) then
              r(ivar,i) = a(ivar,i-1)/a(ivar,i)
            else
              r(ivar,i) = a(ivar,i+1)/a(ivar,i)
            endif
         else
            r(ivar,i) = 0.0
         endif
       enddo
       r(ivar,1) = r(ivar,3)
       r(ivar,2) = r(ivar,3)
       r(ivar,n) = r(ivar,n-2)
       r(ivar,n-1) = r(ivar,n-2)
     enddo
     !
   end subroutine compute_slope
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine roe_average(q,qa,w,n)
      !
      ! This subroutine computes Roe's average qa
      ! of a quantity a using the weights w.
      !
      implicit none
      !
      integer,               intent(in)  :: n
      real,    dimension(n), intent(in)  :: q,w
      real,    dimension(n), intent(out) :: qa
      !
      integer :: i
      !
      do i=2,n-1
        qa(i) = (q(i-1)*w(i-1)+q(i)*w(i))/(w(i-1)+w(i))
      enddo
      !
      qa(1) = qa(2)
      qa(n) = qa(n-1)
      !
   end subroutine roe_average
   !
   !----------------------------------------------------------------------------------------------
   !
   !real function slope_limiter(slope,q1,q2,q3,q4,v)
   !   !
   !   implicit none
   !   !
   !   character(80)    :: slope
   !   real :: q1,q2,q3,q4
   !   real :: v
   !   !
   !   select case(slope)
   !      !
   !   end select
   !   !
   !   slope_limiter = 0.d0
   !   !
   !end function slope_limiter
   real function cfl_timestep(n,nvar,l,dx)
   !
   use RuntimeParameters
   
   implicit none
   !
   integer, intent(in) :: n,nvar
   real, intent(in), dimension(nvar,n) :: l
   real, intent(in)    :: dx
   !
   integer :: i,ivar
   real :: dt,lmax,lmin
   real, dimension(n-1) :: dti
   !
   do i=1,n-1
     lmax = maxval(l(:,i))
     if(lmax.lt.0.0) lmax=0.0
     lmin = minval(l(:,i+1))
     if(lmin.gt.0.0) lmin=0.0
     if(lmax-lmin > 0.0) then
       dti(i) = dx/(lmax-lmin)
     else
       dti(i) = dtmax
     endif
   enddo 
   !
   dt = cfl * minval(dti)
   !
   if(dt < dtmin) then
      write(*,*) 'WARNING: cfl timestep is less than minimum timestep'
      write(*,*) 'using dtmin'
      dt = dtmin
   endif
   !
   cfl_timestep = dt
   !
   end function cfl_timestep
   !
   !----------------------------------------------------------------------------------------------
   !
   real function minmod(a,b)
      !
      implicit none
      !
      real :: a,b,c
      !
      if(a*b.gt.0.d0) then
         if (abs(a).lt.abs(b)) then 
             c = a
         else
             c = b
         endif
      else
         c = 0
      endif
      !
      minmod = c 
      !
   end function minmod
   !
   !----------------------------------------------------------------------------------------------
   !
end module Hydro
