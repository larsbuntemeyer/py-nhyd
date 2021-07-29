!
! Hydro module
!
module Hydro
   !
   implicit none
   !
   public  :: Hydro_init, Hydro_solve
   !
   private :: hydro1D, sweep1D, fill_guardcells, minmod, &
              interface_flux, cfl_timestep
   !
contains
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
   subroutine sweep1D(dt,dx,dir,rho,u,v,w,eint,pres,n,nvar)
      !
      use RuntimeParameters, only: gamma,fl
      use Grid, only: ibg,ieg,nx,ib,ie
      !
      implicit none
      !
      real,    intent(in) :: dx
      real,    intent(in) :: dt
      integer, intent(in) :: dir,n,nvar
      real, dimension(n), intent(in)    :: pres 
      real, dimension(n), intent(inout) :: rho,u,v,w,eint
      !
      ! the current fluxes:
      ! the flux is defined at the cell interfaces, so we
      ! for an x-grid with nx cells, we have nx+2*nguard+1
      ! cell interfaces...
      !
      real, dimension(nvar,n)   :: q
      real, dimension(n)        :: enth
      ! Roe averages at cell interfaces
      real, dimension(n+1)      :: velx_a,vely_a,velz_a
      real, dimension(n+1)      :: enth_a,vel_a,cs_a
      ! average eigenvalues at cell interfaces
      real, dimension(nvar,n+1) :: dq,t,p,r,l,a,K1,K2,K3,K4,K5
      real, dimension(nvar,n+1) :: fl_left,fl_right,fl_diss,fl_roe
      !
      real :: rho2,rho2m1,ei,ekin,vtot2,dq5a,ener,dtdx,pres_a,rho_a
      real :: fx,fy,fz
      !
      integer :: i,j,k,ivar
      !
      ! State Vector
      ! contains conserved variables and enthalpy (including guard cells)
      !
      do i=ibg,ieg
         !
         vtot2 = u(i)**2+v(i)**2+w(i)**2
         ekin  = 0.5 * vtot2
         !
         q(1,i)  = rho(i)
         q(2,i)  = rho(i) * u(i)
         q(3,i)  = rho(i) * v(i)
         q(4,i)  = rho(i) * w(i)
         q(5,i)  = rho(i) * (eint(i) + ekin)
         enth(i) = eint(i) + ekin + pres(i)/rho(i) 
         !
      enddo
      !
      ! Compute jumps in conserved quantities at cell interfaces
      !
      do ivar=1,nvar
        do i=ib,ie+1
          !
          dq(ivar,i) = q(ivar,i) - q(ivar,i-1)
          !
        enddo
        !dq(ivar,1)   = 0.0!dq(ivar,2)
        !dq(ivar,n+1) = 0.0!dq(ivar,n)
      enddo
      !
      ! Compute Roe averages at cell interfaces
      !
      do i=ib,ie+1
        !
        rho2m1 = sqrt(rho(i-1))
        rho2   = sqrt(rho(i))
        pres_a = (pres(i-1)*rho2m1 + pres(i)*rho2) / (rho2m1+rho2)
        rho_a = (rho(i-1)*rho2m1 + rho(i)*rho2) / (rho2m1+rho2)
        velx_a(i) = (u(i-1)*rho2m1 + u(i)*rho2) / (rho2m1+rho2)
        vely_a(i) = (v(i-1)*rho2m1 + v(i)*rho2) / (rho2m1+rho2)
        velz_a(i) = (w(i-1)*rho2m1 + w(i)*rho2) / (rho2m1+rho2)
        enth_a(i) = (enth(i-1)*rho2m1 + enth(i)*rho2) / (rho2m1+rho2)
        vel_a(i) = sqrt(velx_a(i)**2+vely_a(i)**2+velz_a(i)**2)
        !cs_a(i) = sqrt((gamma-1.0)*(enth_a(i)-0.5*vel_a(i)**2))
        !cs_a(i) = (gamma-1.0)*sqrt(enth_a(i)-0.5*vel_a(i)**2)
        cs_a(i) = sqrt(gamma*pres_a/rho_a) 
        !
      enddo
      !
      do i=ib,ie+1
        !
        ! Eigenvalues 
        !
        l(1,i) = velx_a(i) - cs_a(i)
        l(2,i) = velx_a(i)
        l(3,i) = velx_a(i)
        l(4,i) = velx_a(i)
        l(5,i) = velx_a(i) + cs_a(i)
        !
        ! K-Vectors (Eigenvectors)
        !
        K1(1,i) = 1.0
        K1(2,i) = velx_a(i) - cs_a(i)
        K1(3,i) = vely_a(i)
        K1(4,i) = velz_a(i)
        K1(5,i) = enth_a(i) - velx_a(i)*cs_a(i)
        !
        K2(1,i) = 1.0
        K2(2,i) = velx_a(i)
        K2(3,i) = vely_a(i)
        K2(4,i) = velz_a(i)
        K2(5,i) = 0.5*vel_a(i)**2
        !
        K3(1,i) = 0.0
        K3(2,i) = 0.0 
        K3(3,i) = 1.0 
        K3(4,i) = 0.0 
        K3(5,i) = vely_a(i)
        !
        K4(1,i) = 0.0
        K4(2,i) = 0.0 
        K4(3,i) = 0.0 
        K4(4,i) = 1.0 
        K4(5,i) = velz_a(i)
        !
        K5(1,i) = 1.0
        K5(2,i) = velx_a(i) + cs_a(i)
        K5(3,i) = vely_a(i)
        K5(4,i) = velz_a(i)
        K5(5,i) = enth_a(i) + velx_a(i)*cs_a(i)
        !
        ! wavestrengths
        !
        a(3,i) = dq(3,i) - vely_a(i)*dq(1,i)
        a(4,i) = dq(4,i) - velz_a(i)*dq(1,i)
        dq5a   = dq(5,i) - (dq(3,i)-vely_a(i)*dq(1,i))*vely_a(i) - (dq(4,i)-velz_a(i)*dq(1,i))*velz_a(i)
        a(2,i) = (gamma-1.0)/(cs_a(i)**2) * (dq(1,i)*(enth_a(i)-velx_a(i)**2)+velx_a(i)*dq(2,i)-dq5a)
        a(1,i) = 1.0/(2.0*cs_a(i)) * (dq(1,i)*(velx_a(i)+cs_a(i))-dq(2,i)-cs_a(i)*a(2,i))
        a(5,i) = dq(1,i)-(a(1,i)+a(2,i))
        !
      enddo
      !
      ! timestep
      !
      !dt = cfl_timestep(nx+1,nvar,l(:,ib:ie+1),dx)
      !
      ! compute left and right flux vectors at cell interfaces
      !
      do i=ib,ie+1
        fl_left(1,i)  =  q(1,i-1)*u(i-1)
        fl_left(2,i)  = (q(2,i-1)*u(i-1)+pres(i-1))
        fl_left(3,i)  =  q(3,i-1)*u(i-1)
        fl_left(4,i)  =  q(4,i-1)*u(i-1)
        fl_left(5,i)  = (q(5,i-1)+pres(i-1))*u(i-1)
        fl_right(1,i) =  q(1,i)*u(i)
        fl_right(2,i) = (q(2,i)*u(i)+pres(i))
        fl_right(3,i) =  q(3,i)*u(i)
        fl_right(4,i) = q(4,i)*u(i)
        fl_right(5,i) = (q(5,i)+pres(i))*u(i)
      enddo
      !
      ! compute slopes for the flux limiter
      !
      call compute_slope(a,l,n,nvar,r)
      !
      ! compute flip flop function and flux limiter
      !
      do ivar=1,nvar
        do i=ib,ie+1
          ! flip flop
          t(ivar,i) = sign(1.0,l(ivar,i)) 
          ! flux limiter
          p(ivar,i) = phi(fl,r(ivar,i))
        enddo
      enddo
      !
      ! compute the dissipative flux
      !
      dtdx = dt/dx
      do ivar=1,nvar
        do i=ib,ie+1
          fl_diss(ivar,i) = a(1,i)*l(1,i)*K1(ivar,i) * (t(1,i)+p(1,i)*(l(1,i)*dtdx-t(1,i))) &
                          + a(2,i)*l(2,i)*K2(ivar,i) * (t(2,i)+p(2,i)*(l(2,i)*dtdx-t(2,i))) &
                          + a(3,i)*l(3,i)*K3(ivar,i) * (t(3,i)+p(3,i)*(l(3,i)*dtdx-t(3,i))) &
                          + a(4,i)*l(4,i)*K4(ivar,i) * (t(4,i)+p(4,i)*(l(4,i)*dtdx-t(4,i))) &
                          + a(5,i)*l(5,i)*K5(ivar,i) * (t(5,i)+p(5,i)*(l(5,i)*dtdx-t(5,i)))
        enddo
      enddo
      !
      ! compute Roe flux at cell interfaces 
      !
      do ivar=1,nvar
        do i=ib,ie+1
          fl_roe(ivar,i) = 0.5*(fl_left(ivar,i)+fl_right(ivar,i)) - 0.5*fl_diss(ivar,i)
        enddo
      enddo
      !
      ! advect in flux conserving form
      !
      do ivar=1,nvar
        do i=ib,ie+1
          q(ivar,i) = q(ivar,i) + dtdx*(fl_roe(ivar,i)-fl_roe(ivar,i+1))
        enddo
      enddo 
      !
      ! Update simple variables
      !
      do i=ib,ie
         rho(i) = q(1,i)
         u(i)   = q(2,i) / q(1,i)
         v(i)   = q(3,i) / q(1,i)
         w(i)   = q(4,i) / q(1,i)
         ener   = q(5,i) / q(1,i)
         vtot2  = u(i)**2+v(i)**2+w(i)**2
         ekin   = 0.5 * vtot2
         eint(i) = ener - ekin
      enddo
      !
      contains
        !
        ! flip flop function
        !
        real function ff(x)
        implicit none
        real, intent(in) :: x
        ff = sign(1.0,x)
        end function ff
        !
   end subroutine sweep1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_solve(dt)
      !
      use Grid, only: ndim,ib,ie
      use Database
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
      do k=kb,ke
        do j=jb,je
          !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
          !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
          call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),ibg,ieg,ieg,2)
          call sweep1D(0.5*dt,dx,1,dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k), &
                       pres(:,j,k),nx+2*nguard,nvar) 
        enddo  
      enddo  
      !
      do j=jb,je
        do i=ib,ie
          call fill_guardcells_1D(dens(i,j,:),pres(i,j,:),eint(i,j,:),w(i,j,:),kbg,keg,keg,2)
          call sweep1D(dt,dz,1,dens(i,j,:),w(i,j,:),u(i,j,:),v(i,j,:),eint(i,j,:), &
                       pres(i,j,:),nz+2*nguard,nvar) 
        enddo  
      enddo
      ! 
      do k=kb,ke
        do j=jb,je
          !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
          !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
          call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),ibg,ieg,ieg,2)
          call sweep1D(0.5*dt,dx,1,dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k), &
                       pres(:,j,k),nx+2*nguard,nvar) 
        enddo  
      enddo  
      !
   end subroutine Hydro_solve
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine compute_xflux(rho,u,v,w,eint,pres,n,nvar,ib,ie,Fl,Fr)
    implicit none
    integer, intent(in) :: n,nvar,ib,ie
    real, intent(in)   , dimension(n+1)    :: rho,u,v,w,eint,pres
    real, intent(inout), dimension(nvar,n) :: Fl,Fr
    integer :: i
    real    :: enth,ekin,vtot2
    !
    ! compute left and right flux vectors at cell interfaces
    !
    do i=ib,ie+1
     vtot2    = u(i-1)**2+v(i-1)**2+w(i-1)**2
     ekin     = 0.5 * vtot2
     Fl(1,i)  =    rho(i-1)*u(i-1)
     Fl(2,i)  =    rho(i-1)*u(i-1)*u(i-1) + pres(i-1)
     Fl(3,i)  =    rho(i-1)*u(i-1)*v(i-1)
     Fl(4,i)  =    rho(i-1)*u(i-1)*w(i-1)
     Fl(5,i)  =  ( rho(i-1)*(eint(i-1)+ekin) + pres(i-1) ) * u(i-1)
     vtot2    = u(i)**2+v(i)**2+w(i)**2
     ekin     = 0.5 * vtot2
     Fr(1,i)  =    rho(i)*u(i)
     Fr(2,i)  =    rho(i)*u(i)*u(i) + pres(i)
     Fr(3,i)  =    rho(i)*u(i)*v(i)
     Fr(4,i)  =    rho(i)*u(i)*w(i)
     Fr(5,i)  =  ( rho(i)*(eint(i)+ekin) + pres(i) ) * u(i)
    enddo
    ! 
   end subroutine compute_xflux
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine compute_zflux(rho,u,v,w,eint,pres,n,nvar,ib,ie,Fl,Fr)
    implicit none
    integer, intent(in) :: n,nvar,ib,ie
    real, intent(in)   , dimension(n+1)    :: rho,u,v,w,eint,pres
    real, intent(inout), dimension(nvar,n) :: Fl,Fr
    integer :: i
    real    :: enth,ekin,vtot2
    !
    ! compute left and right flux vectors at cell interfaces
    !
    do i=ib,ie+1
     vtot2    = u(i-1)**2+v(i-1)**2+w(i-1)**2
     ekin     = 0.5 * vtot2
     Fl(1,i)  =    rho(i-1)*w(i-1)
     Fl(2,i)  =    rho(i-1)*w(i-1)*u(i-1)
     Fl(3,i)  =    rho(i-1)*w(i-1)*v(i-1)
     Fl(4,i)  =    rho(i-1)*w(i-1)*w(i-1) + pres(i-1)
     Fl(5,i)  =  ( rho(i-1)*(eint(i-1)+ekin) + pres(i-1) ) * w(i-1)
     vtot2    = u(i)**2+v(i)**2+w(i)**2
     ekin     = 0.5 * vtot2
     Fr(1,i)  =    rho(i)*w(i)
     Fr(2,i)  =    rho(i)*w(i)*u(i)
     Fr(3,i)  =    rho(i)*w(i)*v(i)
     Fr(4,i)  =    rho(i)*w(i)*w(i) + pres(i)
     Fr(5,i)  =  ( rho(i)*(eint(i)+ekin) + pres(i) ) * w(i)
    enddo
    ! 
   end subroutine compute_zflux
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine hydro1D(dt)
      !
      use Grid
      use Database
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
      real, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: flux
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: q,dq
      !
      real :: rho,ei,ekin,vtot2
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
      ! Compute jumps in conserved variables
      !
      do k=kbg,keg
         do j=jbg,jeg
            do i=ibg,ieg
               do ivar=1,nvar
                  !
                  dq(ivar,i,j,k) = q(ivar,i,j,k)-q(ivar,i-1,j,k)
                  !
               enddo
            enddo
         enddo
      enddo 
      !
      ! Interpolate velocities at cell interfaces
      !
      do i=ib,ie+1
         do j=jb,je+k2d
            do k=kb,ke+k3d
               !
               !uf(i,j,k) = 0.5d0*(q(2,i,j,k)   / q(1,i,j,k) +  &
               !                    q(2,i-1,j,k) / q(1,i-1,j,k)) 
               !vf(i,j,k) = 0.5d0*(q(3,i,j,k)   / q(1,i,j,k) +  &
               !                    q(3,i,j-1,k) / q(1,i,j-1,k)) 
               !wf(i,j,k) = 0.5d0*(q(4,i,j,k)   / q(1,i,j,k) +  &
               !                    q(4,i,j,k-1) / q(1,i,j,k-1)) 
               uf(i,j,k) = 0.5d0 * (u(i,j,k) + u(i-1,j,k))
               !
            enddo
         enddo
      enddo 
      !
      ! construct the fluxes at cell interfaces       
      !
      do i=ib,ie+1
         do j=jb,je+k2d
            do k=kb,ke+k3d
               !
               flux(1,i,j,k) = interface_flux(q(1,i-2,j,k),q(1,i-1,j,k), &
                                              q(1,i,j,k),  q(1,i+1,j,k), &
                                              uf(i,j,k),dt)
             !  if(uf(i,j,k) < 0.d0) then
             !     flux(1,i,j,k) = uf(i,j,k) *  q(1,i,j,k) 
             !     flux(2,i,j,k) = uf(i,j,k) *  q(2,i,j,k) + pres(i,j,k)
             !     flux(3,i,j,k) = uf(i,j,k) *  q(3,i,j,k)
             !     flux(4,i,j,k) = uf(i,j,k) *  q(4,i,j,k)
             !     flux(5,i,j,k) = uf(i,j,k) * (q(5,i,j,k) + pres(i,j,k))
             !  else 
             !     flux(1,i,j,k) = uf(i,j,k) *  q(1,i-1,j,k) 
             !     flux(2,i,j,k) = uf(i,j,k) *  q(2,i-1,j,k) + pres(i,j,k)
             !     flux(3,i,j,k) = uf(i,j,k) *  q(3,i-1,j,k)
             !     flux(4,i,j,k) = uf(i,j,k) *  q(4,i-1,j,k)
             !     flux(5,i,j,k) = uf(i,j,k) * (q(5,i-1,j,k) + pres(i,j,k))
             !  endif
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
                       dt/(xrCoord(i)-xlCoord(i)) * (flux(ivar,i,j,k) - flux(ivar,i+1,j,k))
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
               !u(i,j,k)   = q(2,i,j,k) / rho
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
   end subroutine hydro1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine fill_guardcells_1D(dens,pres,eint,u,ib,ie,n,bc)
     use Grid, only: nguard
     use RuntimeParameters, only: outflow, periodic
     implicit none
     real ,intent(inout), dimension(n) :: dens,pres,eint,u
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
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(ie-nguard)
           pres(i) = pres(ie-nguard)
           eint(i) = eint(ie-nguard)
           u(i)    = u(ie-nguard)
         enddo
       case(periodic) 
         do i=ib,ib+nguard-1
           dens(i) = dens(ie-nguard-i+1)
           pres(i) = pres(ie-nguard-i+1)
           eint(i) = eint(ie-nguard-i+1)
           u(i)    = u(ie-nguard-i+1)
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(i-n+nguard+1)
           pres(i) = pres(i-n+nguard+1)
           eint(i) = eint(i-n+nguard+1)
           u(i)    = u(i-n+nguard+1)
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
      real :: r,phi,flux,theta
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
      select case(fl)
         !
         case('donor-cell')
            phi = 0.0
         case('Lax-Wendroff')
            phi = 1.0
         case('Beam-Warming')
            phi = r
         case('Fromm')
            phi = 0.5*(1.0+r)
         case('minmod')
            phi = minmod(1.0,r)
         case('superbee')
            phi = max(0.0,min(1.0,2.0*r),min(2.0,r))
         case default
            phi = 0.0
      end select
      !
      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*phi*(q3-q2)
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
