c     Implicite Euler matrix assembler using k_zz corrections
c     for vertical eddy diffusivity corrections.
c
c
c     ===============================================================
      subroutine asmblmat1_corrected_diff (maxmx, mbc, mx, cou 
     &                                    , x_0, x_low, L, moption
     &                                    , SLHS, lda,dx)
c     ===============================================================
c
c     this subroutine assembles the matrix that is required for solution
c     of the one dimensional diffusion solver using Implicit Euler for
c     time stepping with variable diffusion coefficient. the expression
c     for k_zz used for corrections is,
c   
c     k/k0 = z/z0 * ( (1+ z0/L)/ (1+4.7*z/L) )
c 
c     where z0, L, k0 are parameters defined by user. 
c
c     maxmx: maximum number of FV cells
c     mbc: number of ghost cells used for boundary conditions
c     mx: number of cells used for computations
c     moption: (currently a dummy variable) holds options for different
c             storage options, boundary conditions etc.
c     SLHS: the laft hand side matrix stored in LINPACK format. For
c             more info on the storage format look at:
c             www.netlib.org/linpack
c
c
c
c
c     ##################################################################
c     This is a bit sloppy, we should definitely rewrite this in a more
c     'readable' manner
c     #################################################################
c

      implicit double precision (a-h,o-z)

      dimension SLHS(lda,maxmx)
      integer moption
      doubleprecision L, x_0
      dimension moption(2)


      common /diffrobin/ W_dep, W_cor, u_set
      common /diffslv3/ z_0, z_low, alphax, alphay, alphaz, 
     &                  diffL, k_zz_option
      common /winddata/ windvel, winddir, windtime, readwindfromfile,
     &                  verticalcorrection, alphashear, maxwind      
     
c      write(*,*) 'in asmbl mat',  u_set
c      read(*,*)

      do 20 j = 3, (mx-2)
c         find position of faces of the cell
          x_back = (j-0.5)*dx -0.5*dx 
          x_front = (j-0.5)*dx +0.5*dx 
          
c         k values are going to be constant 
c         for x < x_low
          if (x_back .le. x_low) then
             x_back = x_low
          endif
          if (x_front .le. x_low) then 
             x_front = x_low
          endif

c         calculate k*dt/dx^2 on the back face
          cou_back = cou*(x_back/x_0)*( sqrt(1- 15*x_0/L))/
     &     ( sqrt(1 -15*x_back/L))
c         calculate k*dt/dx^2 on the front face
          cou_front = cou*(x_front/x_0)*( sqrt(1- 15*x_0/L))/
     &                ( sqrt(1- 15*(x_front)/L))
c
c
c         assemble
c
          SLHS(2,j) = -1*cou_back
          SLHS(3,j) = cou_back + cou_front + 1
          SLHS(4,j) = -1*cou_front
   20 continue

c     ===============================
c     boundary conditions on x- lower
c     ===============================
c
c     this is a bit sloppy
c      
c     calculate values of k*dt/dx^2 on faces of bottom 
c     cells
c      
c     for cell j = 1 we consider constant k on both faces since we 
c     want diffusivity to be bounded away from zero thus we just need 
c     ks for the top face of first cell and second cell

      x_back = (2-0.5)*dx - 0.5*dx
      x_front = (2 - 0.5)*dx + 0.5*dx

c     we want constat diffusion coefficients for 
c     x < x_low
      if (x_back .le. x_low) then
         x_back = x_low
      endif
          
      if (x_front .le. x_low) then 
         x_front = x_low
      endif

c     k*dt/dx^2 for very bottom face
      cou_low =  cou*(x_low/x_0)*( sqrt(1- 15*x_0/L))/
     &                ( sqrt(1 -4.7*(x_low)/L))
      
      ! also retrieve value of k to find correct
      ! deposition coefficient
      alphaz_low =  alphaz*(x_low/x_0)*( sqrt(1- 15*x_0/L))/
     &               ( sqrt(1- 15*x_low/L))
      gamma = -2*alphaz_low/(dx*(W_dep - u_set))
      beta = (gamma -1)/(1 + gamma)

c     calculate k*dt/dx^2 on the back face
      cou_back = cou*(x_back/x_0)*( sqrt(1- 15*x_0/L))/
     &               ( sqrt(1- 15*x_back/L))
c     calculate k*dt/dx^2 on the front face
      cou_front = cou*(x_front/x_0)*( sqrt(1- 15*x_0/L))/
     &                ( sqrt(1- 15*(x_front)/L))
      
      if (moption(1) .eq. 1) then
c     Dirichlet
         
         SLHS(3,1) = 1 + 2.0*cou_low + cou_back
         SLHS(4,1) = -1*cou_back
         SLHS(2,2) = -1*cou_back
         SLHS(3,2) = cou_back + cou_front +1
         SLHS(4,2) = -1.0*cou_front
      else if (moption(1) .eq. 2) then
c     Neumann
         SLHS(3,1) = 1 + cou_back
         SLHS(4,1) = -1.0*cou_back
         SLHS(2,2) = -1.0*cou_back
         SLHS(3,2) = cou_back + cou_front +1
         SLHS(4,2) = -1.0*cou_front
      else if  (moption(1) .eq. 3) then
c     Robin
c         write(*,*) 'coefficients in robinbc'
c         write(*,*) 'gamma: ',gamma  
c         write(*,*) 'beta: ', beta
c         write(*,*) 'alpha_low: ', alphaz_low
c         write(*,*) 
c         read(*,*)
         SLHS(3,1) =  1 + cou_back + ( 1.0- beta)*cou_low
         SLHS(4,1) = -1.0*cou_back
         SLHS(2,2) = -1.0*cou_back
         SLHS(3,2) = cou_back + cou_front +1
         SLHS(4,2) = -1.0*cou_front
      else   
         write(*,*) 'asmblmat1: Boundary condition on lower face 
     & unidentified'
         read(*,*)
      end if

c     ===============================
c     boundary conditions on x- upper
c     ===============================
c      
c     calculate values of k*dt/dx^2 on faces of top
c     cells
c      
c     for cell j = mx we consider constant k on both faces since we 
c     are assuming we are in the surface layer?
c
c     I am not sure if we even need this since we are not considering 
c     deposition in the upper face of the domain.
c
      

      x_back = (mx-1 - 0.5)*dx - 0.5*dx
      x_front = (mx-1 - 0.5)*dx + 0.5*dx

c     calculate k*dt/dx^2 on the back face
      cou_back = cou*(x_back/x_0)*( sqrt(1- 15*x_0/L))/
     &               ( sqrt(1- 15*x_back/L))
c     calculate k*dt/dx^2 on the front face
      cou_front = cou*(x_front/x_0)*( sqrt(1- 15*x_0/L))/
     &                ( sqrt(1- 15*(x_front)/L))  

      if (moption(2) .eq. 1) then
c     Dirichlet
         SLHS(2, mx-1) = -1.0*cou_back
         SLHS(3, mx-1) = cou_front + cou_back +1
         SLHS(2, mx-1) = -1.0*cou_front
         SLHS(2, mx) = -1.0*cou_front
         SLHS(3, mx) = 1 + 3.0*cou_front

      else if (moption(2) .eq. 2) then
c     Neumann
         SLHS(2, mx-1) = -1.0*cou_back
         SLHS(3, mx-1) = cou_back + cou_front  +1
         SLHS(4, mx-1) = -1.0*cou_front
         SLHS(2, mx) = -1.0*cou_front
         SLHS(3, mx) = 1 + cou_front

      else
         write(*,*) 'asmblmat1: Boundary condition unidentified'
         read(*,*)
      end if

      return
      end
