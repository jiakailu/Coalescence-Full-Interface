	PROGRAM JET Inclined Spines

C LAST MODIFICATION 08-Apr-2013 JK Lu Total Number of Element Max 105
C EL CODIGO ESTA BASADO EN AIRW.FOR Y LAS MODIFICACIONES REALIZADAS
C ESTAN SE¥ALADAS CON LA LEYENDA "jet"

C   *******************************************************************
C   *   This numerical code solves axisimetric Navier-Stokes          *
C   *   and continuity                                                *
C   *   equations for a transient 2-D free surface (optional) flow    *
C   *   problem. It is particularly applied to follow the time evol-  * 
C   *   ution of jet breakup caused by capillary forces               *
C   *                                                                 *
C   *      # The spatial discretization of equations is accomplished  *
C   *        by the Galerkin/finite element method and a suitable     *
C   *        parameterization of the free surface.                    *
C   *                                                                 *
C   *      # Time discretization is done by means of a second order   *
C   *        finite-difference integrator. Epsil is the local time    *
C   *        truncation error (it must be entered as a data).         *
C   *                                                                 *
C   *      # The non-linear resulting equation set is solved by means *
C   *        of full Newton method or one-step Newton method (The op- *
C   *        tion is taken through the parameter INEW).               *
C   *   ------------------------------------------------------------- *
C   *   Physico-chemical properties of the liquid correspond to       *
C   *                                                                 *
C   *            viscosity (rmu) = 10^-3 kg/m.s                       *
C   *            density (rho) = 4275 kg/m^3                          * 
C   *            air-liquid surface tension (sigma) = 0.02 N/m        * 
C   *   ------------------------------------------------------------- *
C   *   The initialization can be done by:                            *
C   *       a)Imposing a sinusoidal perturbation of amplitude epsilon * 
C   *         and wavenumber wno to an initally motionless liquid     *
C   *         film. In this case, time must be set equal to zero      * 
C   *         and epsilon and wno are required parameters.            *
C   *       b)Imposing the solution corresponding to a given value    *
C   *         of time (that is, the solution corresponding to the     *
C   *         final time of another "run"). In this case, a file      *   
C   *         with the corresponding solution, the initial value of   *
C   *         time and the two last time steps are required.          *
C   *  -------------------------------------------------------------- *
C   *  The code consists in a main part and in a number of special-   *
C   *  ized subroutines.                                              *
C   *  In the main part the following steps are accomplished:         *
C   *    1- Initialization parameters are introduced as data.         *
C   *    2- Dimensionless numbers are evaluated:                      *
C   *                RE, the Reynolds number.                         *
C   *                WE, the Capilar number.                          *
C   *                F,  the dimensionless initial radius             *
C   *    3- Boundary conditions are imposed.                          *
C   *    4- Time integration is performed.                            *
C   *    5- The global Jacobian matrix (EQ(i,j)) is built from the    *
C   *       elementary ones and then this matrix is modified (eqns.   *
C   *       corresponding to boundary conditons are not included in   *
C   *       the final system of eqns.)                                *
C   *    6- The error corresponding to Newton-Raphson method is eval- *
C   *       uated (when full Newton method is employed) in order to   *
C   *       verify if the desired convergence error has been achieved.*
C   *                                                                 *
C   *   The subroutines are:                                          *
C   *    INPUT: reads initial data, builds the connectivity matrix,    *
C   *           evaluates the total number of elements, equations...  * 
C   *    INPUT0 is equivalent to input but it is used when the        *
C   *           initial solution is the corresponding one to time=0.  *
C   *           nex is the number of elements in the x-direction.     *
C   *           ney is the number of elements in the y-direction.     *
C   *           nnx is the number of nodes in the x-direction.        *
C   *           nny is the number of nodes in the y-direction.        *
C   *           nhadd is the total number of nodes (including the     *
C   *           free surface nodes).                                  *
C   *           ne is the total number of elements.                   *
C   *           np is the total number of unknowns.                   *
C   *           maxit is the maximum number of iterations for each    *
C   *           time step when full Newton method is employed.        *
C   *           rmsmax is the required error for Newton method.       *        
C   *    MESH: builds the mesh.                                       *
C   *    DERJAC: evaluates the derivatives of the spatial variables   *
C   *            with respect to the free surface location (h(i)).    *
C   *    DERIV: evaluates the derivatives of the unknowns with        *
C   *             respect to zita (X) and eta (Y).                    *
C   *    ABFIND: evaluates the elementary jacobian matrix (A(i,j))    *
C   *            and the elementary residual vectors (RHS(i)). The    * 
C   *            global residual vector (R2(NP)) is also built here.  *
C   *            Total kinetic and surface energies and the volume    *
C   *            of the liquid layers are also evaluated. They are    *
C   *            written in for036
C   *    GAUSSM: solves the system of equations.                      *
C   *    FIDA: builds the file in which all data (number of elements, *
C   *          value of dimensionless numbers, converged solution for *
C   *          imposed final time, etc.). This file is the one to be  *
C   *          used when b) initialization is used.                   *
C   *    SHORPTR: builds the 'readable' final file, that is, the con- *
C   *             verged solution for the imposed final time and all  *
C   *             the required information: mesh characteristics,     *
C   *             value of dimensionless numbers, time history, etc.  *
C   *    TFUNCT: evaluates the biquadratic and bilinear basis         *
C   *            functions  that expand velocities, pressure and free *
C   *            surface coefficients for specified values of sita    *
C   *            (X) and eta (Y).                                     *
C   *******************************************************************
  
        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	INTEGER W1,W2,W3,W4
	DIMENSION JMOD(500),DSOLV(15999),NCOD(15999)
        DIMENSION vdel(41)
        DIMENSION V1(9),V2(9),V3(4)
	COMMON/RES/R2(15999)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON1/NELL
        COMMON/FSLAB/IELFS(1700), IFZONE(1700)
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
	COMMON/FRON3/A(25,25)
	COMMON/FRON4/EQ(15999,15999),NK(25)
	COMMON/ITER/ITE,MAXIT,RMS,RMSMAX,AMAXUP
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/MESCO/YPT(7011),XPT(7011)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2			!JK_new_mesh
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2		!JK_new_mesh
        COMMON/MESHWT/WY(41),WX(180)			
	COMMON/PARAM/RE,F,OME,WE
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/TEMP2/SOLV(15999),DSOL(15999)
	COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME  !moidified for continuation run: JK
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh 
        common/temf/solp(15999)
c        common/test/TESM(22,22),r2tes(22)
        common/init/epsilon,wno
	COMMON/FACT/IVUEL
	COMMON/NONNEWPARA/alp,beta,pn			!JK_NONNEWTONIAN MODEL PARAMETERS
        COMMON/ZMESH/wxz(15)
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1 !JKMoving
	COMMON/oup/NOoP(1700,12)   !JKMoving
        DATA V1,V2,V3/1,4,6,9,11,13,15,18,20,
     1             2,5,7,10,12,14,16,19,21,
     1             3,8,17,22/

	OPEN(29,file='input.txt',status='old')
	OPEN(69,file='output.sal',status='new')
        write(69,*) ' STEPT, TIMEF (final time)'
	read(29,*) STEPT
	read(29,*)TIMEF
        write(69,*)' epsil'
        read(29,*)epsil
	write(69,*) ' IFS (=0, THE F.S. IS FIXED)'
	read(29,*) IFS
        write(69,*)' inew=?, inew=1, 1 step Newton is applied'
        read(29,*)inew
        write(69,*)' inital value of time is...'
        read(29,*)time
	OPEN(10,file='ouxyi.sal',status='new')
	OPEN(4,file='salxyi.sal',status='new')
	OPEN(44,file='velocity.sal',status='new')
        OPEN(38,file='viscosity.sal',status='new')
        OPEN(444,file='pressure.sal',status='new')
        OPEN(555,file='pressure1.sal',status='new')
        OPEN(212,file='initialvelocity.dat',status='old')
	OPEN(39,file='velotension.sal',status='new')	
        write(39,*) 'TITLE="velotension"'
        write(39,*) 'Variable="X","Y","Kappa1","Kappa2","kappa","Vn"
     1,"Vt"' 
C       **************************************************************
C       IVUEL is an index used to evaluate the vector NO in Subroutine
C       Gaussm. If IVUEL=1, NO is evaluates, otherwise the vector is 
C       used as it is.
	IVUEL=1
C       Initial time step. This is the value of the time step for
C       the first four time steps.

	H0=126.32d-09						
	RHO=4275						
	RMU=1.0D-03
	SIGMA=0.02					
        WE=1.0d+00						! jet
	Re=0
        IWRIT=0
          ITRACK=1
	pn=1
	beta=0.001
        alp=10
c	Read(*,*)RE
        OME=0							! MODIFICADO
	write(69,*) ' RE ',RE,' WE ',WE,' F ',F

C       Input is the subroutine that reads initial data and
C       builds the connectivity matrix, the solution vector ...
c         write(69,*)' epsilon and wno are...'
         read(29,*)epsilon
c		 read(29,*)wno
		
         CALL INPUT0
c        end if
	if(time .ne. 0.d0)then
         CALL CONTR
c        DO I=1,NHADD
c        IF(MDF(I) .EQ. 2)THEN
c        SOL(NOPP(I))=SOL(NOPP(I))*0.2/dabs(SOL(NOPP(NHADD-NNY)))
c        SOL(NOPP(I)+1)=SOL(NOPP(I)+1)*0.2/dabs(SOL(NOPP(NHADD-NNY)))
c        END IF
c        IF(MDF(I) .EQ. 3)THEN
c        SOL(NOPP(I))=SOL(NOPP(I))*0.2/dabs(SOL(NOPP(NHADD-NNY)))
c        SOL(NOPP(I)+1)=SOL(NOPP(I)+1)*0.2/dabs(SOL(NOPP(NHADD-NNY)))
c        END IF
c        END DO
        end if 
	if(time .eq. 0.d0)then
cINITIAL CONDITION FOR THE WHOLE SHEET
        DO I=1,NHADD
        IF(MDF(I) .EQ. 2)THEN
c        READ(212,*)SOL(NOPP(I)),SOL(NOPP(I)+1)
        SOL(NOPP(I)+1)=0
        SOL(NOPP(I))=0 
        END IF
        IF(MDF(I) .EQ. 3)THEN
c        READ(212,*)SOL(NOPP(I)),SOL(NOPP(I)+1)
c        READ(313,*)SOL(NOPP(I)+2)
        SOL(NOPP(I)+1)=0
        SOL(NOPP(I))=0 
        END IF
        END DO
cINITIAL CONDITION FOR THE WHOLE SHEET
	Close(212)
c        Close(313)
	end if
c	if(time .eq. 0.d0)then
c        DO I=1,NHADD
c        IF(MDF(I) .EQ. 2)THEN
c        READ(212,*)SOL(NOPP(I)),SOL(NOPP(I)+1)
c        SOL(NOPP(I))=SOL(NOPP(I))*0.2/SOL(NOPP(NHADD-NNY-1))
c        SOL(NOPP(I)+1)=SOL(NOPP(I)+1)*0.2/SOL(NOPP(NHADD-NNY-1))
c        END IF
c        IF(MDF(I) .EQ. 3)THEN
c        READ(212,*)SOL(NOPP(I)),SOL(NOPP(I)+1)
c        SOL(NOPP(I))=SOL(NOPP(I))*0.2/SOL(NOPP(NHADD-NNY-1))
c        SOL(NOPP(I)+1)=SOL(NOPP(I)+1)*0.2/SOL(NOPP(NHADD-NNY-1))
c        END IF
c        END DO
c	Close(212)
c	end if

Ccccccccccccc   Fixed boundary conditions.ccccccccccccccccccc
C         a) b.c. at y=0, v=0.uy=0	
c	  DO I=NNXZ1*(NNYZ1+1)+(NNY+1)*(NNX1-2)+1,NHADD,NNY+1	
c	  DO I=(NNY+1)*(NNX1-1)+1,NHADD,NNY+1	
c	  SOL(NOPP(I))=0.		
c	  SOL(NOPP(I)+1)=0.
c	  NCOD(NOPP(I))=1	
c	  NCOD(NOPP(I)+1)=1
c	  END DO


c         b)b.c. at x=xw, v=0, u=inlet.		
c	  NNI=(NNX-1)*(NNY+1)
c	  kk=0
c	  DO I=NNI+1,NHADD-1
c	  SOL(NOPP(I))=-0.25
c	  SOL(NOPP(I))=-0.35*(1-(kk*(1./(2*ney)))**2)
c	  NCOD(NOPP(I))=1
c          SOL(NOPP(I)+1)=0.
c          NCOD(NOPP(I)+1)=1
c	  kk=kk+1
c	  END DO


C         c)b.c. at x=0, u=0, !Jet Impingement
          DO I=1,NNYZ1,1
c	  DO I=1,(NNY+1)*(NNX1-1)+1,NNY+1	!regular JetJet
c	  DO I=1,NNY,1	
c	  SOL(NOPP(I))=0.
c	  NCOD(NOPP(I))=1
	  SOL(NOPP(I)+1)=0. !CFI
	  NCOD(NOPP(I)+1)=1 !CFI
	  END DO  	!JK 02-11-2013 test sheet half

          DO I=NNXZ1*(NNYZ1+1)+1,NNXZ1*(NNYZ1+1)
     1         +(NNY+1)*(NNX1-2)+1,NNY+1
c	  SOL(NOPP(I))=0.
c	  NCOD(NOPP(I))=1
	  SOL(NOPP(I)+1)=0.  !CFI
	  NCOD(NOPP(I)+1)=1  !CFI
	  END DO  	!JK 02-11-2013 test sheet half

c	  d)b.c. at y=yw 
c          DO I=1,NNY
c          SOL(NOPP(I))=0.
c          NCOD(NOPP(I))=1
c          END DO


c	  e)b.c. fix nozzel h 		! jet on wall
c        sol(nopp(nhadd))=1.
c        ncod(nopp(nhadd))=1 !CFI
	if(ifs .eq. 0)then
	do i=1,nhadd
        if(MDF(i).eq.1)then
	ncod(nopp(i))=1
        end if
c	sol(nopp(i-1))=0.
c	ncod(nopp(i-1))=1
c	sol(nopp(i-1)+1)=0.
c	ncod(nopp(i-1)+1)=1.
	end do
        do i=nnxz1*(nnyz1+1)+nny+1, nhadd, nny+1
         ncod(nopp(i))=1
        end do
	end if


c         Do I=1,((NNX-1)/2)*(NNY+1)
c            If (mdf(i).ne.1) then
c             SOL(NOPP(I))=0.375
c             SOL(NOPP(I)+1)=0.
c            End If
c         End Do                         !JK01-31-2013 test for sheet formation

c         Do I=((NNX-1)/2+1)*(NNY+1)+1, NHADD-1
c          Do I=NNY+2,NHADD-1
c            If (mdf(i).ne.1) then
c             SOL(NOPP(I))=-30
c             SOL(NOPP(I)+1)=0.          !JK02-11-2013 test for sheet half 
c            End If
c         End Do
c	ncod(nopp(nny+1))=1

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



C       Time integration is performed here:
c        if(time .ne. 0.d0)then
c        write(69,*)' stept and stpetv are...'
c        read(*,*)stept,steptv
c        end if
	if(time .eq. 0.d0)then
        ntime=time/stept
        if(ntime .le. 5)ntime=time/stept
        if(ntime .gt. 5)ntime=5
	end if
        int=0
81      NTIME=NTIME+1
        invu=(ntime+1)/51
        rnvu=(ntime+1)/51.
        if(dabs(invu-rnvu) .le.  1.d-10)then
        write(4,201)ntime
201	format(//10x,' stept values, ntime =',i4,//)
        do j=1,5
        i1=(j-1)*10+1
        i2=i1+10
        write(4,202)(vdel(i),i=i1,i2)
        end do
202     format(5x,11e11.4)
        do i=1,50
        vdel(i)=0.
        end do
        int=0
        end if
        IF(NTIME .EQ. 1)THEN
        ct1=1.
        ct2=0.
        int=int+1
        vdel(int)=stept
        kh=0
        do i=1,nhadd
        k=nopp(i)
        if(mdf(i) .ne. 1)then
         solv(k)=sol(k)
         sol(k)=solv(k)+stept*dsol(k)
         solv(k+1)=sol(k+1)
         sol(k+1)=solv(k+1)+stept*dsol(k+1)
        else
         kh=kh+1
         solv(k)=sol(k)
         sol(k)=solv(k)+stept*dsol(k)
         h(kh)=sol(k)
        end if
        end do
         solv(NP+1)=sol(NP+1)   !JK-artificial moving front
c         sol(NP+1)=solv(NP+1)+stept*dsol(NP+1)  !JK-artifical moving front
          sol(NP+1)=-stept*solv(nopp(nnyz1))+solv(NP+1) !CFI
        END IF
        IF(NTIME .GE. 2 .AND .NTIME .LT. 5)THEN
        ct1=1.
        ct2=0.
        int=int+1
        vdel(int)=stept
        kh=0
        do i=1,nhadd
        k=nopp(i)
        if(mdf(i) .ne. 1)then
         dsolv(k)=dsol(k)
         dsol(k)=(sol(k)-solv(k))/stept
         temp=solv(k)
         solv(k)=sol(k)
         sol(k)=2.*solv(k)-temp
         dsolv(k+1)=dsol(k+1)
         dsol(k+1)=(sol(k+1)-solv(k+1))/stept
         temp=solv(k+1)
         solv(k+1)=sol(k+1)
         sol(k+1)=2.*solv(k+1)-temp
        else
         kh=kh+1
         dsolv(k)=dsol(k)
         dsol(k)=(sol(k)-solv(k))/stept
         temp=solv(k)
         solv(k)=sol(k)
         sol(k)=2.*solv(k)-temp
         h(kh)=sol(k)
        end if
        end do
         dsolv(NP+1)=dsol(NP+1)  !JK-artificial moving front
         dsol(NP+1)=(sol(NP+1)-solv(NP+1))/stept  !JK-artificial moving front
         temp=solv(NP+1)  !JK-artificial moving front
         solv(NP+1)=sol(NP+1)  !JK-artificial moving front
          sol(NP+1)=-stept*solv(nopp(nnyz1))+solv(NP+1) !CFI
c         sol(NP+1)=2.*solv(NP+1)-temp  !JK-artificial moving front
c          sol(NP+1)=stept*solv(nopp(nnyz1)+1)+solv(NP+1)
c          sol(NP+1)=-stept*solv(nopp(nnyz1)+1)+solv(NP+1)
        END IF
        IF(NTIME .EQ. 5)THEN
        ct1=2.
        ct2=-1.
        int=int+1
        vdel(int)=stept
        kh=0
        steptv=stept
        do i=1,nhadd
        k=nopp(i)
         if(mdf(i) .ne. 1)then
          dsolv(k)=dsol(k)
          dsol(k)=2.*(sol(k)-solv(k))/stept-dsolv(k)
          solv(k)=sol(k)
          sol(k)=solv(k)+stept/2.*(3.*dsol(k)-dsolv(k))
          solp(k)=sol(k)
          dsolv(k+1)=dsol(k+1)
          dsol(k+1)=2.*(sol(k+1)-solv(k+1))/stept-dsolv(k+1)
          solv(k+1)=sol(k+1)
          sol(k+1)=solv(k+1)+stept/2.*(3.*dsol(k+1)-dsolv(k+1))
          solp(k+1)=sol(k+1)
         else
          kh=kh+1
          dsolv(k)=dsol(k)
          dsol(k)=2.*(sol(k)-solv(k))/stept-dsolv(k)
          solv(k)=sol(k)
          sol(k)=solv(k)+stept/2.*(3.*dsol(k)-dsolv(k))
          solp(k)=sol(k)
          h(kh)=sol(k)
         end if
        end do
          dsolv(NP+1)=dsol(NP+1)!JK-artificial moving front
          dsol(NP+1)=2.*(sol(NP+1)-solv(NP+1))/stept-dsolv(NP+1)!JK-artificial moving front
          solv(NP+1)=sol(NP+1)!JK-artificial moving front
          sol(NP+1)=-stept*solv(nopp(nnyz1)+1)+solv(NP+1)
c          sol(NP+1)=-1.5*stept*solv(nopp(nnyz1)+1)+solv(NP+1)
c          sol(NP+1)=solv(NP+1)+stept/2.*(3.*dsol(NP+1)-dsolv(NP+1))!JK-artificial moving front
          solp(NP+1)=sol(NP+1)!JK-artificial moving front
        END IF
        IF(NTIME .GE. 6)THEN
        ct1=2.
        ct2=-1.
        int=int+1
        kh=0
        sumu=0.
        sumv=0.
        sumh=0.
        umax=0.
        vmax=0.
        hmax=0.
        do i=1,nhadd
        k=nopp(i)
        if(mdf(i) .ne. 1)then
        sumu=sumu+(sol(k)-solp(k))**2
        if(dabs(sol(k)) .ge. umax)umax=dabs(sol(k))
        sumv=sumv+(sol(k+1)-solp(k+1))**2
        if(dabs(sol(k+1)) .ge. vmax)vmax=dabs(sol(k+1))
        else
        sumh=sumh+(sol(k)-solp(k))**2
        if(dabs(sol(k)) .ge. hmax)hmax=dabs(sol(k))
        end if
        end do
        rnorv=dsqrt((sumu/umax**2+sumv/vmax**2)/(nnx*nny))
        rnorh=dsqrt((sumh/hmax**2)/nnx)
        temp=stept
        stept1=stept*(3.*epsil*(1.+steptv/stept)/rnorv)**(1./3.)
        stept2=stept*(3.*epsil*(1.+steptv/stept)/rnorh)**(1./3.)
        write(69,*)stept1,stept2
        stept=stept1
        if(stept2 .lt. stept)stept=stept2
        if(stept .gt. (steptv*2.))stept=2.*steptv
        steptvv=steptv
        steptv=temp
        if((stept/steptv) .lt. 1. .and. (stept/steptv) .ge. 0.8)
     1   stept=steptv
        ioa=0

c        if((stept/steptv) .lt. 0.8)then
c         stept=0.8*steptv

c         time=time-steptv
c         steptv=steptvv
c         int=int-1
c         ioa=1
c        end if 
	if(stept.ge.0.001)then
	stept=0.001
	end if				!JK max timestep
       
        vdel(int)=stept
         coe1=solv(nopp(nnyz1))   !coe
        do i=1,nhadd
        k=nopp(i)
         if(mdf(i) .ne. 1)then
         if(ioa .ne. 1)then
         dsolv(k)=dsol(k)
         dsol(k)=2./steptv*(sol(k)-solv(k))-dsolv(k)
         solv(k)=sol(k)
         end if
         sol(k)=solv(k)+stept/2.*((2.+stept/steptv)*dsol(k)-
     1           stept/steptv*dsolv(k))
         solp(k)=sol(k)
         if(ioa .ne. 1)then
         dsolv(k+1)=dsol(k+1)
         dsol(k+1)=2./steptv*(sol(k+1)-solv(k+1))-dsolv(k+1)
         solv(k+1)=sol(k+1)
         end if
         sol(k+1)=solv(k+1)+stept/2.*((2.+stept/steptv)*dsol(k+1)-
     1              stept/steptv*dsolv(k+1))
         solp(k+1)=sol(k+1)
         else
          kh=kh+1
          if(ioa .ne. 1)then
          dsolv(k)=dsol(k)
          dsol(k)=2./steptv*(sol(k)-solv(k))-dsolv(k)
          solv(k)=sol(k)
          end if
          sol(k)=solv(k)+stept/2.*((2.+steptv/stept)*dsol(k)-
     1              steptv/stept*dsolv(k))
          solp(k)=sol(k)
          h(kh)=sol(k)
         end if
        end do
cccccccccccccccc  CFI  ccccccccccccccccccccccccccc
          if(ioa .ne. 1)then
          dsolv(NP+1)=dsol(NP+1)!JK-artificial moving front
          dsol(NP+1)=2./steptv*(sol(NP+1)-solv(NP+1))-dsolv(NP+1)!JK-artificial moving front
          solv(NP+1)=sol(NP+1)!JK-artificial moving front
          end if
          
          coe=solv(nopp(nnyz1))/coe1
         if (solv(nopp(nnyz1)).GT.0) then
          if(solv(nopp(nnyz1+1)).GE.0.8*0.2.AND.solv(nopp(nnyz1
     1     +1)).LE.1.2*0.2) then
          coe=coe
          else if(solv(nopp(nnyz1+1)).LT.0.8*0.2) then
          coe=1.2*coe
          else if(solv(nopp(nnyz1+1)).GT.1.2*0.2) then
          coe=0.8*coe
          end if
         else 
          if(solv(nopp(nnyz1+1)).GE.0.8*0.2.AND.solv(nopp(nnyz1
     1     +1)).LE.1.2*0.2) then
          coe=coe
          else if(solv(nopp(nnyz1+1)).LT.0.8*0.2) then
          coe=0.8*coe
          else if(solv(nopp(nnyz1+1)).GT.1.2*0.2) then
          coe=1.2*coe
          end if
         end if
          sol(NP+1)=-coe*stept*solv(nopp(nnyz1))+solv(NP+1)
c          sol(NP+1)=solv(NP+1)+stept/2.*((2.+steptv/stept)*dsol(NP+1)-!JK-artificial moving front
c     1              steptv/stept*dsolv(NP+1))!JK-artificial moving front
          solp(NP+1)=sol(NP+1)!JK-artificial moving front
        END IF
ccccccccccccccccc CFI ccccccccccccccccccccccccccccccc
	TIME=TIME+STEPT

cccccccccccc  Defining the Inlet Velocity BC ccccccccccccccccc
c          RE= 2000+100*TIME
c          IF (RE.GT.20000)then
c          RE=108000
c          END IF
c          RE=108000
          IF(ITRACK.EQ.1)then
          TIMEI=TIME
          END IF
          AAA=0.2
          ITRACK=1+ITRACK
c          AAA=0.2+0.05*(TIME-TIMEI)
          if ((time-timei).ge.5) then
          AAA=0.2+0.01*dsin((2*3.14/10)*(TIME-timei-5))
c         AAA=0.4
          end if
c          write(*,*)(time-timei)
c          if (zzz.ge.1) then
          zzz=1
c          end if
c          NNI=(NNX-1)*(NNY+1)
          NNI=NHADD-NNY-1
          kk=1
c	  AAA=0.35
          DO I=NNI+1,NHADD-1
c           Do I=1, NNY, 1
          SOL(NOPP(I))=0.
c          SOL(NOPP(I))=-AAA*(1-zzz*(wy(kk))**2)
          NCOD(NOPP(I))=1
c          SOL(NOPP(I)+1)=0.
c          NCOD(NOPP(I)+1)=1
          kk=kk+1
          END DO                !JK_OSCILATING_INLET

c          DO I= 1, NNY
c          SOL(NOPP(I))=AAA*(1-zzz*(wy(kk))**2)
c          NCOD(NOPP(I))=1
c          SOL(NOPP(I)+1)=0.
c          NCOD(NOPP(I)+1)=1
c          kk=kk+1
c          END DO   !JK01-31-2013 test for sheet formation 


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        write(36,98)ntime,time,stept,steptv
	ITE=0
10      ITE=ITE+1
	write(69,*) ITE, Re
	CALL MESH
C       The global Jacobian matrix is built here:  
	  NELL=0
	  DO J=1,NP
	     R2(J)=0.!Nullify the Global Residual Vector
             DO I=1,NP
		EQ(I,J)=0.!Nullify the Global Jacobian matrix
	     END DO
	  END DO

18        NELL=NELL+1
	  call abfind
          if(IFZONE(NELL).EQ.2)THEN
            DO i=1,25 
             A(i,23)=0.d0
             A(i,24)=0.d0
             A(i,25)=0.d0
            END DO
           end if
	  N=NELL
	  DO I=1,25
	  NK(I)=0
	  END DO
	  KC=0
	  JEND=NBN(N)
	    DO J=1,JEND
	    M=NOP(N,J)
	    K=NOPP(M)
	    IDF=MDF(M)
		DO I=1,IDF
		KC=KC+1
		II=K+I-1
		NK(KC)=II
		END DO
	    END DO
	  IF(JEND .NE. 12)THEN
	   DO NO=1,3
	   NK(22+NO)=W4(NO)
	   END DO
	  END IF
	  LEND=25
	    DO L=1,LEND
	    LL=NK(L)
	    KEND=NCN(N)
		DO K=1,KEND
		KK=NK(K)
		EQ(KK,LL)=EQ(KK,LL)+A(K,L)
		END DO
	    END DO

	  IF(NELL .LT. NE)GO TO 18
	  write(69,561) NP
c	Do i=1,NP 
c	EQ(NP,i)=0				! forcing last two h equal
c	End Do
c        EQ(NP,NP)=1
c	EQ(NP,nopp(nhadd-nny-1))=-1
c	R2(NP)=0
C         **************************************************
	  IR=0
	  DO K=1,NP
	  IF(NCOD(K) .EQ. 1)THEN
	  IR=IR+1
	  JMOD(IR)=K
	  END IF
	  END DO
	  write(69,*) IR
	  NCC=IR
C         *******************************************************
C         Rows corresponding to boundary conditions are modified:
	  DO IRR=1,IR
	  K=JMOD(IRR)
	  R2(K)=0.
	    DO L=1,NP
	    EQ(K,L)=0.
	    IF(L .EQ. K)EQ(K,L)=1.
	    END DO
	  END DO
C       **************************************************************
C       Rows and columns corresponding to b.c. are put down and wright
C       in the Jacobian matrix in order to exclude them when solving
C       the equation system. 
	  K=0
	  KS=0
	  DO I=1,NP
	  IF(NCOD(I) .EQ. 1)THEN
	  K=K+1
	  ELSE
	  KS=I-K
	  R2(KS)=R2(I)
	     DO J=1,NP
	     EQ(KS,J)=EQ(I,J)
	     END DO
	  END IF
	  END DO
	  K=0
	  KS=0
	  DO I=1,NP
	  IF(NCOD(I) .EQ. 1)THEN
	  K=K+1
	  ELSE
	  KS=I-K
	     DO J=1,NP
	     EQ(J,KS)=EQ(J,I)
	     END DO
	  END IF
	  END DO
C       NCC= number of boundary conditions.
	NPM=NP-NCC
C       NPM is the number of equations to be solved by GAUSSM.
	write(69,524) NP,NPM
	CALL GAUSSM(EQ,R2,NPM)
	KS=NPM+1
        if(ncod(1).NE.1)then    !JK FIX
        temp=R2(1)              !JK FIX
        end if                  !JKFIX
	DO I=NP,1,-1
	IF(NCOD(I) .NE. 1)THEN
	KS=KS-1
	R2(I)=R2(KS)
	R2(KS)=0.
	END IF
	END DO
        if(ncod(1).ne.1)then !JK Fix
        R2(1)=temp  !JK FIX
        end if !JK FIX
	RMS=0.
	AMAXUP=0.
	DO J=1,NP
	IF(NCOD(J) .EQ. 1)THEN 
	R2(J)=0.
	ELSE
 	  SOL(J)=SOL(J)+R2(J)
	  R22=R2(J)**2
	  IF(R22 .GT. AMAXUP)THEN
	  IRE=J
	  AMAXUP=R22
	  END IF
	RMS=RMS+R22
	END IF
	END DO
300     format(5x,' ite:',i4,' i:',i4,
     1	' r2(i):',e15.5,' sol(i):',e15.5)
C       Errors evaluation:
        IF(inew .ne. 1)then
           AMAXUP=DSQRT(AMAXUP)
	   RMS=DSQRT(RMS)
	   if(ITE .EQ. 1)RNORMA=RMS
	   write(69,13) AMAXUP,RMS,IRE
        Else 
           rms=0.d0
        End if
 
C       H(I) is here updated.
	K=0
	DO I=1,NHADD
        IF(MDF(I).EQ.1)THEN
	KN=NOPP(I)
	K=K+1
	H(K)=SOL(KN)
	write(69,91) I,K,KN,SOL(NOPP(I))
	IF(H(K) .LT. 0.)STOP' H < 0'
        END IF
	END DO
	IF(ITE .GE. MAXIT)THEN
	write(69,*) ITE,MAXIT
	write(69,100)
c        call fida
c        call shortpr
	END IF
	IF(RMS .GT. RMSMAX .AND. ITE .LT. MAXIT)THEN
	write(69,11) ITE,RMS
	   GO TO 10
	END IF
	IF(RMS .LE. RMSMAX)THEN
c         write(36,*)h(1),h(nnx)
          if(time .ge. timef)then
            write(4,201)ntime
            do j=1,5
               i1=(j-1)*10+1
               i2=i1+10
               write(4,202)(vdel(i),i=i1,i2)
            end do
c            call fida
	    if(NTIME.GE.6)THEN
	    call CONTW
	    end if
            call shortpr

	  else
            ivuel=1
	    if(NTIME.GE.6)THEN
	    call CONTW
	    end if
            IWRIT=IWRIT+1
            IF(IWRIT.EQ.2.OR.IWRIT.EQ.1)then
            call shortpr
            IWRIT=1
            END IF
            call calvelot
            ivuel=1  !what is ivuel ??
	    go to 81
	  end if
	END IF
	 Close(29)
	 Close(69)
	 Close(44)
	 Close(444)
	 Close(555)
         Close(39)
98      FORMAT(5X,' NTIME=',I4,' TIME=',E20.10,' STEPT=',E20.10,
     1            ' STEPTV=',E20.10)
99      FORMAT(10X,' MAIN ITE:',I4)
561     FORMAT(10X,' The evaluation of EQ(i,j) is finished, NP =',I5)
524     FORMAT(10X,' NP:',I5,' NPM:',I5)
13      FORMAT(10X,' AMAXUP =',E20.10,' RMS =',E20.10,' IRE =',I5)
91      FORMAT(10X,' I:',I4,' K:',I4,' KN:',I4,' H:',E20.10)
11      FORMAT(10X,' ITE =',I5,' ERROR =',E20.10)
100     FORMAT(///,'  NO CONVERGENCE')
	STOP
	END


C 	*****************************************

	SUBROUTINE CONTR

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON1/NELL
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/TEMP2/SOLV(15999),DSOL(15999)
        common/temf/solp(15999)
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
        COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME
	
	OPEN(529,file='continuation.sal',status='OLD')
	read(529,*)STEPT
	read(529,*)STEPTV
	read(529,*)TIME
	read(529,*)NTIME
	DO i=1,np+1
           read(529,*)sol(i)
        END DO
        kh=1
        DO i=1,nhadd
           if(mdf(i).eq.1)then
           k=nopp(i)
           h(kh)=sol(k)
           kh=kh+1
           end if
        END DO
        DO i=1,np+1
           read(529,*)solv(i)
        END DO
        DO i=1,np+1
           read(529,*)solp(i)
        END DO
        DO i=1,np+1
           read(529,*)dsol(i)
        END DO

	Close(529) 
	
	RETURN
	END
	

	SUBROUTINE CONTW
	

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON1/NELL
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/TEMP2/SOLV(15999),DSOL(15999)
        common/temf/solp(15999)
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
        COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME
	
	OPEN(529,file='continuation.sal',status='REPLACE')
	write(529,*) STEPT
	write(529,*) STEPTV
	write(529,*) TIME
	write(529,*) NTIME
	DO i=1,NP+1
	   write(529,*)SOL(i)
	END DO 
	DO i=1,NP+1
	   write(529,*)SOLV(i)
	END DO 
	DO i=1,NP+1
	   write(529,*)SOLP(i)
	END DO 
	DO i=1,NP+1
	   write(529,*)DSOL(i)
	END DO
	Close(529) 
	RETURN
	END

	SUBROUTINE INPUT0

C Last modification 15-may-98
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/FSLAB/IELFS(1700), IFZONE(1700)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
	COMMON/MESPA/XW,YW
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/ITER/ITE,MAXIT,RMS,RMSMAX,AMAXUP
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2
	COMMON/MESHWT/WY(41),WX(180)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON1/NELL
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
        COMMON/INIT/epsilon,wno
	INTEGER W1,W2,W3,W4
	COMMON/PARAM/RE,F,OME,WE	! dimensionless parameter OME??
	COMMON/NMESH/DIR(180,2)				!Added: Jun-2010 By: JK Store the sin and cos of the angle of each spine
        COMMON/NONNEWPARA/alp,beta,pn 
        COMMON/ZMESH/wxz(15) !JKMoving
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1 !JKMoving
	COMMON/oup/NOoP(1700,12)   !JKMoving
C 	*****************************************
	OPEN(129,file='xcoor.dat',status='OLD') 
	OPEN(229,file='ycoor.dat',status='OLD') 
	OPEN(629,file='xcoorz.dat',status='OLD') !JKMoving 
	OPEN(329,file='direction.dat',status='OLD')
	OPEN(429,file='initialh.dat',status='OLD')  
	Write(*,*)'number of x1 elements'	!JK_MESh
	Read(29,*)nex1
	Write(*,*)'number of x2 elements'
	Read(29,*)nex2
	nex=nex1+nex2
	Write(*,*)'number of y elements'
        Read(29,*)ney
        Write(*,*)'xw'
	Read(29,*)xw
	Write(*,*)'yw'
	Read(29,*)yw
	Write(*,*)'RMSMAX'
        Read(29,*)rmsmax
	Write(*,*)'Maximum of newton iteration'
        Read(29,*)maxit
	Write(*,*)'Npatternx'
	Read(29,*)Npatternx
	Write(*,*)'Npatterny'
	Read(29,*)Npatterny
        Read(29,*)NEXZ1!JKMoving
        NEYZ1=NEY-NEXZ1!JKMoving
        NEZ1=NEYZ1*NEXZ1!JKMoving
        NNXZ1=2*NEXZ1+1!JKMoving
        NNYZ1=2*NEYZ1+1!JKMoving
	nnx1=2*nex1+1				!JK_mesh
	nnx2=2*nex2				!JK_mesh
	NNX=2*NEX+1
	NNY=2*NEY+1
c*************************
	NE=NEX*NEY+NEXZ1*NEYZ1 !JKMoving
	NH=NNX*NNY+NNXZ1*NNYZ1-NNY !JKMoving
	NHADD=NH+NNX+NNXZ1-1 !JKMoving
	NP=2*NH+NNX+(NEX+1)*(NEY+1)
     1    +(NEYZ1+1)*(NEXZ1+1)-(NEY+1)+(NNXZ1-1)  !JKMoving
c**************************
c	NE=NEX*NEY
c	NH=NNX*NNY
c	NHADD=NH+NNX
c	NP=2*NH+NNX+(NEX+1)*(NEY+1)
c**************************
c        WRITE(18,3)NEX,NEY,NNX,NNY
c        WRITE(18,4)NH,NHADD,NP
3       FORMAT(10X,' NEX =',I4,' NEY =',I4,' NNX =',I5,' NNY =',I5)
4       FORMAT(10X,' NH =',I5,' NHADD =',I5,' NP =',I5)
        delx1=1./(nnx1-1)
	delx2=1./nnx2				!JK_mesh
	dely=1./(nny-1)
	IF(Npatternx.eq.1)then
         DO I=1,NNX1
          wx(i)=(i-1)*delx1
         END DO
          wx(1)=0
	 DO I=1,NNX2
	  j=i+nnx1
	  wx(j)=i*delx2*xw
	 END DO
	Else
         DO I=1,NNX
          Read(129,*)wx(i)
c          wx(i)=wx(i)+0.3
         END DO
	END IF
        IF(Npatterny.eq.1)then
         DO I=1,NNY
          wy(i)=(i-1)*dely
         END DO
        Else
         DO I=1,NNY-NNXZ1+1
          Read(229,*)wy(i)
         END DO
        END IF

c*************************************!JKMoving
        DO I=1,NNXZ1
          Read(629,*)wxz(i) 
        END DO
c*************************************!JKMoving


	DO I=1,NNX+NNXZ1-1 !JKMoving
	 Read(329,*)dir(i,1),dir(i,2)
c         dir(i,1)=0
c         dir(i,2)=1 !flow in the tube
	END DO
        DO I=1,NE
	  IELFS(I)=0
	  NBN(I)=0
	  NCN(I)=0
	END DO
c*************************************!JKMoving
        DO I=NEYZ1,NEZ1,NEYZ1
          IELFS(I)=1
        END DO
	DO I=NEY+NEZ1,NE,NEY
	  IELFS(I)=1
	END DO
c*************************************!JKMoving

c*************************************!JKMoving
        DO I=1,NEZ1
          IFZONE(I)=1
        END DO

        DO I=1,NEX
           DO J=1,NEY
             IF(J.LE.NEXZ1)THEN
               IFZONE((I-1)*NEY+J+NEZ1)=2
             ELSE 
               IFZONE((I-1)*NEY+J+NEZ1)=3
             END IF
           END DO
        END DO
c*************************************!JKMoving

	DO 13 I=1,NEXZ1
	  DO 13 J=1,NEYZ1
	    NEL=NEYZ1*(I-1)+J
	    DO K=1,3 
	      NOOP(NEL,K)=(NNYZ1)*(2*I+K-3)+2*J-1 
	      NOOP(NEL,K+3)=NOoP(NEL,K)+1
	      NOoP(NEL,K+6)=NOoP(NEL,K)+2
	    END DO
13      CONTINUE

	DO  I=1,1
	  DO  J=1,NEXZ1
	    NEL=NEY*(I-1)+J+NEZ1
	      DO K=1,1 
	        NOoP(NEL,K)=NOoP((J-1)*NEYZ1+1,1)
	        NOoP(NEL,K+3)=NOoP((J-1)*NEYZ1+1,2)
	        NOoP(NEL,K+6)=NOoP((J-1)*NEYZ1+1,3)
	      END DO
	      DO K=2,3 
	        NOoP(NEL,K)=(NNY)*(2*I+K-3-1)+2*J-1+NOoP(NEZ1,9) 
	        NOoP(NEL,K+3)=NOoP(NEL,K)+1
	        NOoP(NEL,K+6)=NOoP(NEL,K)+2
	      END DO
           END DO
	  DO  J=NEXZ1+1,NEY
	    NEL=NEY*(I-1)+J+NEZ1
	      DO K=1,1 
	        NOoP(NEL,K)=NOoP(J+NEYZ1*(NEXZ1-1)-NEXZ1,3)
	        NOoP(NEL,K+3)=NOoP(J+NEYZ1*(NEXZ1-1)-NEXZ1,6)
	        NOoP(NEL,K+6)=NOoP(J+NEYZ1*(NEXZ1-1)-NEXZ1,9)
	      END DO
	      DO K=2,3 
	        NOoP(NEL,K)=(NNY)*(2*I+K-3-1)+2*J-1+NOoP(NEZ1,9) 
	        NOoP(NEL,K+3)=NOoP(NEL,K)+1
	        NOoP(NEL,K+6)=NOoP(NEL,K)+2
	      END DO
           END DO
        END DO

	DO  I=2,NEX
	  DO  J=1,NEY
	    NEL=NEY*(I-1)+J+NEZ1
	      DO K=1,3 
	        NOoP(NEL,K)=(NNY)*(2*I+K-3-1)+2*J-1+NOoP(NEZ1,9) 
	        NOoP(NEL,K+3)=NOoP(NEL,K)+1
	        NOoP(NEL,K+6)=NOoP(NEL,K)+2
	      END DO
           END DO
        END DO

C 	Connectivity matrix

	DO 12 I=1,NEXZ1
	  DO 12 J=1,NEYZ1
 	    NEL=NEYZ1*(I-1)+J
 	    DO K=1,3 
	      NOP(NEL,K)=(NNYZ1+1)*(2*I+K-3)+2*J-1 
	      NOP(NEL,K+3)=NOP(NEL,K)+1
	      NOP(NEL,K+6)=NOP(NEL,K)+2
 	    END DO
	  IF(IELFS(NEL) .EQ. 1)THEN
 	    DO K=1,3
	      NOP(NEL,K+9)=NOP(NEL,6+K)+1
 	    END DO 
	  END IF
12      CONTINUE

	DO  I=1,1
	  DO  J=1,NEXZ1
 	    NEL=NEY*(I-1)+J+NEZ1
 	      DO K=1,1 
	        NOP(NEL,K)=NOP((J-1)*NEYZ1+1,1)
	        NOP(NEL,K+3)=NOP((J-1)*NEYZ1+1,2)
	        NOP(NEL,K+6)=NOP((J-1)*NEYZ1+1,3)
	      END DO
	      DO K=2,3 
	        NOP(NEL,K)=(NNY+1)*(2*I+K-3-1)+2*J-1+NOP(NEZ1,12) 
	        NOP(NEL,K+3)=NOP(NEL,K)+1
	        NOP(NEL,K+6)=NOP(NEL,K)+2
	      END DO
           END DO
	  DO  J=NEXZ1+1,NEY
 	    NEL=NEY*(I-1)+J+NEZ1
	      DO K=1,1 
	        NOP(NEL,K)=NOP(J+NEYZ1*(NEXZ1-1)-NEXZ1,3)
	        NOP(NEL,K+3)=NOP(J+NEYZ1*(NEXZ1-1)-NEXZ1,6)
	        NOP(NEL,K+6)=NOP(J+NEYZ1*(NEXZ1-1)-NEXZ1,9)
	      END DO
	      DO K=2,3 
	        NOP(NEL,K)=(NNY+1)*(2*I+K-3-1)+2*J-1+NOP(NEZ1,12) 
	        NOP(NEL,K+3)=NOP(NEL,K)+1
	        NOP(NEL,K+6)=NOP(NEL,K)+2
	      END DO
	     IF(IELFS(NEL) .EQ. 1)THEN
	       DO K=1,3
	         NOP(NEL,K+9)=NOP(NEL,6+K)+1
	       END DO 
	     END IF
           END DO
        END DO

	DO  I=2,NEX
	  DO  J=1,NEY
 	    NEL=NEY*(I-1)+J+NEZ1
	      DO K=1,3 
	        NOP(NEL,K)=(NNY+1)*(2*I+K-3-1)+2*J-1+NOP(NEZ1,12) 
	        NOP(NEL,K+3)=NOP(NEL,K)+1
	        NOP(NEL,K+6)=NOP(NEL,K)+2
	      END DO
           END DO
	   IF(IELFS(NEL) .EQ. 1)THEN
	      DO K=1,3
	        NOP(NEL,K+9)=NOP(NEL,6+K)+1
	      END DO 
	   END IF
        END DO
c        WRITE(18,100)
c        DO I=1,NE
c          WRITE(18,102)I
c          WRITE(18,101)(NOP(I,J),J=1,12)
c        END DO
100     FORMAT(//,10X,' CONNECTIVITY MATRIX',//)
101     FORMAT(10X,12I4)
102     FORMAT(10X,' NEL:',I4)

C 	Degrees of freedom

C       NBN(NEL)=I, I is the number of nodes of element nel
C       NCN(NEL)=J, J is the number of equations of element nel.
	DO NEL =1,NE 
	  NBN(nel)=9
	  NCN(NEL)=22
	  DO NODO=1,9
	    K=NOP(NEL,NODO)
	    MDF(K)=3
	  END DO
	  DO NODO=2,8,2
	    K=NOP(NEL,NODO)
	    MDF(K)=2
	  END DO
	  MDF(NOP(NEL,5))=2
	  IF(IELFS(NEL) .EQ. 1)THEN
	    DO NODO=10,12
	      K=NOP(NEL,NODO)
	      MDF(K)=1
	    END DO
	    NCN(NEL)=25
	    NBN(NEL)=12
	  END IF
	END DO

C 	Nopp(nhadd)

	NOPP(1)=1
	J=0
	DO I=2,NHADD
	  J=J+MDF(I-1)
	  NOPP(I)=J+1
	END DO

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c
C 	Perturbacion senoidal inicial        

C	pis2=3.14159/2.d0			! modificado

        do i=1,nnxz1
	  read(429,*)h(i)
          write(69,*)i,h(i)
        end do  
        do i=nnxz1+1,nnx+nnxz1-1
c          ax=wno*(i-1)*delx*xw
c          h(i)=F*(1.+epsilon*dcos(ax))		! modificacion
	  read(429,*)h(i)
c           h(i)=1
c            h(i)=0.5/wx(i+1-nnxz1)
          write(69,*)i,h(i)
        end do
        kh=1
        do i=1,nhadd
          if(MDF(i).eq.1)then
          k=nopp(i)
          sol(k)=h(kh)
          kh=kh+1
          end if
        end do
        Sol(NP+1)=WX(NNX1)-WX(1)!ADDADDADDADD
	Close(129)
	Close(229)
	Close(329)
	Close(429)
	Close(629) !JKM
	RETURN
	END


C	**********************************************

	SUBROUTINE MESH

C 	Last modification 09-may-98

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
	COMMON/MESCO/YPT(7011),XPT(7011)
	COMMON/MESPA/XW,YW
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2
	COMMON/MESHWT/WY(41),WX(180)
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh 
        COMMON/NONNEWPARA/alp,beta,pn
        COMMON/ZMESH/wxz(15) !JKMoving
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1 !JKMoving
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/oup/NOoP(1700,12)
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
C       Base points      
        DO I=1,NNXZ1
c          BCP(I,1)=wxz(i)
c          BCP(I,2)=-sol(NP+1)+wx(NNX1)
          BCP(I,2)=wxz(i) !CFI
          BCP(I,1)=-sol(NP+1)+wx(NNX1) !CFI
        END DO   
 
        DO I=1,NNX1-1
          j=I+NNXZ1
c	  BCP(j,2)=-(wx(I+1)-wx(NNX1))*
c     1    sol(NP+1)/(wx(1)-wx(NNX1))
c     2    +wx(NNX1)
c	  BCP(j,1)=0.0d0				! JK_mesh
	  BCP(j,1)=-(-wx(I+1)+wx(NNX1))*
     1    sol(NP+1)/(-wx(1)+wx(NNX1))
     2    +wx(NNX1)
	  BCP(j,2)=0.0d0				! CFI 
        END DO

	DO I=1,NNX2
	  j=I+NNX1
	  BCP(j,1)=wx(j)
	  BCP(j,2)=0.0d0				! JK_mesh
	END DO

        DO I=1,NNXZ1
          DO J=1,NNYZ1
             NODO=(NNYZ1+1)*(I-1)+J
	     XPT(NODO)=BCP(I,1) + H(I)*WY(J)*DIR(I,1)          ! JK_new_mesh
	     YPT(NODO)=BCP(I,2) + H(I)*WY(J)*DIR(I,2)	! JK_new_mesh
          END DO
        END DO
C       Mesh coordinates
	DO I=2,NNX
	  DO J=1,NNXZ1
	    NODO=(NNY+1)*(I-2)+J+NNXZ1*(NNYZ1+1)
c	    XPT(NODO)=BCP(I+NNXZ1-1,1) 
c     1       + wxz(NNXZ1)*(wxz(j)/wxz(NNXZ1))         ! JK_new_mesh
c	    YPT(NODO)=BCP(I+NNXZ1-1,2)! JK_new_mesh
	    YPT(NODO)=BCP(I+NNXZ1-1,2) 
     1       + wxz(NNXZ1)*(wxz(j)/wxz(NNXZ1))         ! CFI
	    XPT(NODO)=BCP(I+NNXZ1-1,1)! CFI
	  END DO
	  DO J=NNXZ1+1,NNY
	    NODO=(NNY+1)*(I-2)+J+NNXZ1*(NNYZ1+1)
c	    XPT(NODO)=BCP(I+NNXZ1-1,1)
c     1     + wxz(NNXZ1)
c     2     + H(I+NNXZ1-1)*WY(J-NNXZ1+1)*DIR(I+NNXZ1-1,1)          ! JK_new_mesh
c	    YPT(NODO)=BCP(I+NNXZ1-1,2) 
c     1     + H(I+NNXZ1-1)*WY(J-NNXZ1+1)*DIR(I+NNXZ1-1,2)	! JK_new_mesh
	    XPT(NODO)=BCP(I+NNXZ1-1,1)   !CFI
     1     + H(I+NNXZ1-1)*WY(J-NNXZ1+1)*DIR(I+NNXZ1-1,1)          !CFI
	    YPT(NODO)=BCP(I+NNXZ1-1,2) !CFI
     1     + wxz(NNXZ1)   !CFI
     2     + H(I+NNXZ1-1)*WY(J-NNXZ1+1)*DIR(I+NNXZ1-1,2)	! CFI
	  END DO
	END DO        
85      FORMAT(10X,' NODE:',I4,' XPT(I):',E20.10,' YPT(I):',E20.10)
86      FORMAT(//,10X,' MESH COORDINATES',//)
	
	RETURN
	END

c       ****************************************

	SUBROUTINE DERIV

C	Last modification 09-may-98

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/TFUN1/PHI(9),PSI(4),PHIC(9),PHIE(9),PHICC(9)
	COMMON/DER1/AXPT(9),AYPT(9),WYL(3),DIRL(3,2),WXL(3)
	COMMON/DER2/U,UC,UE,V,VC,VE,P,rs                !rs: second invariant of the 
	COMMON/DER3/XE,XC,YC,YE,DET,Y			! modificado
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/TEMP1/UV,VV,UPV,VPV,DYT1,DYT2,DXT1,DXT2,DYH(3)
	COMMON/TEMP2/SOLV(15999),DSOL(15999)    !*V: previous solution
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1 !JKMoving
        COMMON/ZMESH/wxz(15) 
	COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME
	COMMON/FSLAB/IELFS(1700),IFZONE(1700)
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
	COMMON/FRON1/NELL
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	INTEGER W1,W2,W3,W4

	XC=0.
	YC=0.
	YE=0.
	XE=0.
	Y =0.						! modificado
	DO I=1,9
	  XC=XC+AXPT(I)*PHIC(I)
	  YC=YC+AYPT(I)*PHIC(I)
	  YE=YE+AYPT(I)*PHIE(I)
	  XE=XE+AXPT(I)*PHIE(I)        !JK_Det_test
	  Y =Y +AYPT(I)*PHI(I) 				! modificado
	END DO
	DET = -XC*YE+XE*YC					!JK_Det_test
	U=0.
	V=0.
	P=0.
	UC=0.
	UE=0.
	VC=0.
	VE=0.
	UV=0.
	VV=0.
	UPV=0.
	VPV=0.
         rs=0.
 	DO I=1,9
	  U=U+SOL(W1(I))*PHI(I)
	  V=V+SOL(W2(I))*PHI(I)
	  UC=UC+SOL(W1(I))*PHIC(I)
	  UE=UE+SOL(W1(I))*PHIE(I)
	  VC=VC+SOL(W2(I))*PHIC(I)
	  VE=VE+SOL(W2(I))*PHIE(I)
	  UV=UV+SOLV(W1(I))*PHI(I)
	  VV=VV+SOLV(W2(I))*PHI(I)
	  UPV=UPV+DSOL(W1(I))*PHI(I)
	  VPV=VPV+DSOL(W2(I))*PHI(I) 
c          rs=rs+SLvis(W5(i))*PHI(I)
	END DO
	DO I=1,4
	  P=P+SOL(W3(I))*PSI(I)
	END DO
	DYT1=0.
	DYT2=0.
	DXT1=0.   !JK
	DXT2=0.   !JK
	DO K=1,3
          dyh(k)=0.
          DO L=1,3
            dyt1=dyt1+phi(3*(l-1)+k)*dirl(k,2)*wyl(l)*(sol(w4(k))
     & 	    -solv(w4(k))) 				! jet
            dyt2=dyt2+phi(3*(l-1)+k)*wyl(l)*dirl(k,2)*dsol(w4(k))	! jet
            dyh(k)=dyh(k)+wyl(l)*dirl(k,2)*phi(3*(l-1)+k)       	! jet
	    dxt1=dxt1+phi(3*(l-1)+k)*dirl(k,1)*wyl(l)*(sol(w4(k))
     &      -solv(w4(k))) 				! JK_X_Mesh_Movement
	    dxt2=dxt2+phi(3*(l-1)+k)*dirl(k,1)*wyl(l)*dsol(w4(k))	! JK_X_Mesh_Movement
	  END DO					
	END DO	
        DO K=1,3
          DO L=1,3
            dxt1=dxt1-phi(3*(l-1)+k)*wxl(k)*(sol(NP+1)
     & 	    -solv(NP+1)) 				
            dxt2=dxt2-phi(3*(l-1)+k)*wxl(k)*dsol(NP+1)	
c            dyt1=dyt1-phi(3*(l-1)+k)*wxl(k)*(sol(NP+1)
c     & 	    -solv(NP+1)) 				
c            dyt2=dyt2-phi(3*(l-1)+k)*wxl(k)*dsol(NP+1)	
	  END DO					
	END DO
	RETURN
	END

C       ******************************************************

	SUBROUTINE DERJAC

c	Last modification 09-may-98

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/TFUN1/PHI(9),PSI(4),PHIC(9),PHIE(9),PHICC(9)
	COMMON/DERJA/DYEDH(3),DYCDH(3),DH(3),DHE(3),DHC(3)
        COMMON/DER1/AXPT(9),AYPT(9),WYL(3),DIRL(3,2),WXL(3)
	COMMON/DER3/XE,XC,YC,YE,DET,Y			! modificado
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 
	DO I=1,3
	  DH(I)=0.
	  DHC(I)=0.
	  DHE(I)=0.
	END DO
	DO K=1,3
	  DO L=1,3
	    dh(k)=dh(k)+phi(3*(l-1)+k)*wyl(l)           ! JK 
	    dhc(k)=dhc(k)+phic(3*(l-1)+k)*wyl(l)        ! JK
	    dhe(k)=dhe(k)+phie(3*(l-1)+k)*wyl(l)        ! JK
	  END DO
	END DO
	RETURN
	END

C       *****************************************

	SUBROUTINE TFUNCT(X,Y)

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/TFUN1/PHI(9),PSI(4),PHIC(9),PHIE(9),PHICC(9)
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 
	RL1(C)=2.*C**2-3.*C+1.
	RL2(C)=4.*C-4.*C**2
	RL3(C)=2.*C**2-C
	DL1(C)=4.*C-3.
	DL2(C)=4.-8.*C
	DL3(C)=4.*C-1.
        D2L1=4.                    !01-17-2013 JK , for calculation curvature for output
        D2L2=-8.                   !01-17-2013 JK , for calculation curvature for output
        D2L3=4.                    !01-17-2013 JK , for calculation curvature for output
	PHI(1)=RL1(X)*RL1(Y)
	PHI(2)=RL2(X)*RL1(Y)
	PHI(3)=RL3(X)*RL1(Y)
	PHI(4)=RL1(X)*RL2(Y)
	PHI(5)=RL2(X)*RL2(Y)
	PHI(6)=RL3(X)*RL2(Y)
	PHI(7)=RL1(X)*RL3(Y)
	PHI(8)=RL2(X)*RL3(Y)
	PHI(9)=RL3(X)*RL3(Y)
	PHIC(1)=RL1(Y)*DL1(X)
	PHIC(2)=RL1(Y)*DL2(X)
	PHIC(3)=RL1(Y)*DL3(X)
	PHIC(4)=RL2(Y)*DL1(X)
	PHIC(5)=RL2(Y)*DL2(X)
	PHIC(6)=RL2(Y)*DL3(X)
	PHIC(7)=RL3(Y)*DL1(X)
	PHIC(8)=RL3(Y)*DL2(X)
	PHIC(9)=RL3(Y)*DL3(X)
	PHIE(1)=RL1(X)*DL1(Y)
	PHIE(2)=RL2(X)*DL1(Y)
	PHIE(3)=RL3(X)*DL1(Y)
	PHIE(4)=RL1(X)*DL2(Y)
	PHIE(5)=RL2(X)*DL2(Y)
	PHIE(6)=RL3(X)*DL2(Y)
	PHIE(7)=RL1(X)*DL3(Y)
	PHIE(8)=RL2(X)*DL3(Y)
	PHIE(9)=RL3(X)*DL3(Y)
	PSI(1)=(1.-Y)*(1.-X)
	PSI(2)=X*(1.-Y)
	PSI(3)=(1.-X)*Y
	PSI(4)=X*Y
        PHICC(1)=RL1(Y)*D2L1          !01-17-2013 JK , for calculation curvature for output
        PHICC(2)=RL1(Y)*D2L2          !01-17-2013 JK , for calculation curvature for output
        PHICC(3)=RL1(Y)*D2L3          !01-17-2013 JK , for calculation curvature for output
        PHICC(4)=RL2(Y)*D2L1          !01-17-2013 JK , for calculation curvature for output
        PHICC(5)=RL2(Y)*D2L2          !01-17-2013 JK , for calculation curvature for output
        PHICC(6)=RL2(Y)*D2L3          !01-17-2013 JK , for calculation curvature for output
        PHICC(7)=RL3(Y)*D2L1          !01-17-2013 JK , for calculation curvature for output
        PHICC(8)=RL3(Y)*D2L2          !01-17-2013 JK , for calculation curvature for output
        PHICC(9)=RL3(Y)*D2L3          !01-17-2013 JK , for calculation curvature for output
	RETURN
	END

C       *****************************************

	SUBROUTINE GAUSSM (BB,BS,NP)

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION BB(15999,15999),BS(15999),NO(15999)
	DIMENSION NPIVOT(15999),NC(15999)
	COMMON/FACT/IVUEL
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 

        KV=0
	IF(IVUEL .EQ. 1)THEN
	DO 1 J=1,NP
	   DO I=NP,1,-1
	   K=NP-I+1
	   IF(dabs(BB(I,J)) .gt. 1.d-99)THEN
	    NO(J)=K
	    GO TO 1
	   END IF
	   END DO
1       CONTINUE
	DO 10 I=1,NP
	   DO J=NP,1,-1
	   IF(dabs(BB(I,J)) .gt. 1.d-99)THEN
	   NC(I)=J
	   GO TO 10
	   END IF
	   END DO
10      CONTINUE
	IVUEL=2
	END IF
	DO K=1,NP-1
	NPIVOT(K)=K
	AMAX=DABS(BB(K,K))
	IF(K .GT. 1)THEN
	 IFIV=IFI
	 IFI=NP-NO(K)+1 
	 IF (IFI .LT. IFIV)IFI=IFIV
	ELSE
	 IFI=NP-NO(K)+1
	END IF
	  DO I=K+1,IFI 
	  ABSA=DABS(BB(I,K))
	  IF(ABSA .GT. AMAX)THEN
	   NPIVOT(K)=I
	   AMAX=ABSA
	  END IF
	  END DO
	IF(NPIVOT(K) .NE. K)THEN
	 KV=KV+1
	 I=NPIVOT(K)
	 TEM=BS(K)
	 BS(K)=BS(I)
	 BS(I)=TEM
	 KN=NC(K)
	 NC(K)=NC(I)
	 NC(I)=KN
	   DO J=K,NP
	   TEM=BB(K,J)
	   BB(K,J)=BB(I,J)
	   BB(I,J)=TEM
	   END DO
	END IF

	IFC=NC(K)
	  DO I=K+1,IFI
	  IF(NC(I) .GT. IFC)THEN
	   IFC=NC(I)
	  ELSE
	   NC(I)=IFC
	  END IF
          IF(dabs(BB(I,K)).GT.1.d-99)then
	  RMUL=BB(I,K)/BB(K,K)
	  BB(I,K)=RMUL
	  BS(I)=BS(I)-RMUL*BS(K)
	    DO J=K+1,IFC
            if(dabs(BB(K,J)).GT.1.d-99)then
	    BB(I,J)=BB(I,J)-RMUL*BB(K,J)
            end if
	    END DO
          END IF
	  END DO
	END DO
	DO I=NP,1,-1
	SUM=0.
	I1=I+1
	IF(I1 .LE. NP)THEN
	 DO J=I1,NP
	 SUM=SUM+BB(I,J)*BS(J)
	 END DO
	END IF
	BS(I)=(BS(I)-SUM)/BB(I,I)
	END DO
c       Determinant of BB(i,j)
c        nkv=kv/2
c        rkv=kv/2.
c        sign=1.
c        if(dabs(rkv-nkv) .gt. 1.d-10)sign=-1.
c        det=1.
c        do i=1,np
c        det=det*bb(i,i)/dabs(bb(i,i))
c        end do
c        det=det*sign
c        write(69,11)det
c11	format(5x,' det bb(i,i)=',e12.5)
	RETURN 
	END

C       ****************************************

	SUBROUTINE FIDA

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/ITER/ITE,MAXIT,RMS,RMSMAX,AMAXUP
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2
	COMMON/PARAM/RE,F,OME,WE
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/MESPA/XW,YW
	COMMON/MESHWT/WY(41),WX(180)
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
        common/temp2/solv(12265),dsol(12265)
        common/temf/solp(15999)
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 
	INTEGER W1,W2,W3,W4
C
	OPEN(10,file='ouxyi.sal',status='OLD')
	WRITE(10,1001)NEY,NEX
	WRITE(10,1002)RE,WE,F,OME
	WRITE(10,1006)XW
	WRITE(10,1004)RMSMAX,MAXIT
	DO I=1,NNX
	WRITE(10,1006)WX(I)
	END DO
	DO I=1,NNY
	WRITE(10,1006)WY(I)
	END DO
        do i=1,np
        write(10,1006)sol(i)
        end do
        do i=1,np
        write(10,1006)solv(i)
        end do
        do i=1,np
        write(10,1006)solp(i)
        end do
        do i=1,np
        write(10,1006)dsol(i)
        end do
1001    FORMAT(10X,2I4)
1002    FORMAT(10X,4E20.10)
1004    FORMAT(10X,E20.10,I4)
1006    FORMAT(10X,E20.10)
	RETURN
	END
C       **********************************
	SUBROUTINE SHORTPR
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/PARAM/RE,F,OME,WE
	COMMON/MESPA/XW,YW
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2
	COMMON/ITER/ITE,MAXIT,RMS,RMSMAX,AMAXUP 
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
	COMMON/MESCO/YPT(7011),XPT(7011)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
        COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME
	COMMON/NMESH/DIR(180,2)				!JK_new_mesh
        COMMON/NONNEWPARA/alp,beta,pn 
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1
        COMMON/ZMESH/wxz(15) 
	COMMON/oup/NOoP(1700,12)
	INTEGER W1,W2,W3,W4
	OPEN(4,file='salxyi.sal',status='OLD')
	OPEN(44,file='velocity.sal',status='OLD')
        OPEN(444,file='pressure.sal',status='OLD')
        OPEN(555,file='pressure1.sal',status='OLD')
c	if(TIME.ge.TIMEF)then
c	WRITE(4,1001)RE,WE,F,OME
c	WRITE(4,1002)
c	WRITE(4,1003)NEY,NEX
c	WRITE(4,1004)XW
c	WRITE(4,1005)TIME,RMS,ITE
c	DO I=1,NHADD
c	IF(MDF(I) .EQ. 2)THEN
c	UW=SOL(NOPP(I))
c	VW=SOL(NOPP(I)+1)
c	WRITE(4,1006)I,XPT(I),YPT(I),UW,VW
c	END IF
c	IF(MDF(I) .EQ. 3)THEN
c	UW=SOL(NOPP(I))
c	VW=SOL(NOPP(I)+1)
c	PW=SOL(NOPP(I)+2)
c	WRITE(4,1007)I,XPT(I),YPT(I),UW,VW,PW
c	END IF
c	END DO
c	end if
	
	IF(TIME.eq.STEPT)THEN
c	WRITE(44,1009)TIME,2*NEY+1,2*NEX+1
	WRITE(44,1009)TIME,NOop(NE,9),4*NE
	ELSE 
c	WRITE(44,1011)TIME,2*NEY+1,2*NEX+1
	WRITE(44,1011)TIME,Noop(NE,9),4*NE
	END IF
        DO I=1,NHADD
          if(MDF(i).ne.1)then
             UW=SOL(NOPP(I))
             VW=SOL(NOPP(I)+1)
            write(44,1006) XPT(i),YPT(i),UW,VW
          end if
        END DO     
        DO I=1,NEXZ1
           DO J=1,NEYZ1
             k=J+(i-1)*NEYZ1
             write(44,*) NoOP(k,1), NOoP(k,2), NoOP(k,5), NOoP(k,4)  
             write(44,*) NoOP(k,4), NOoP(k,5), NOoP(k,8), NOoP(k,7) 
           END Do 
           DO J=1,NEYZ1
             k=J+(i-1)*NEYZ1
             write(44,*) NoOP(k,2), NOoP(k,3), NoOP(k,6), NOoP(k,5)  
             write(44,*) NoOP(k,5), NOoP(k,6), NOoP(k,9), NOoP(k,8) 
           END Do 
        END DO 
        DO I=1,NEX
           DO J=1,NEY
             k=J+(i-1)*NEY+NEZ1
             write(44,*) NoOP(k,1), NOoP(k,2), NoOP(k,5), NOoP(k,4)  
             write(44,*) NoOP(k,4), NOoP(k,5), NOoP(k,8), NOoP(k,7) 
           END Do 
           DO J=1,NEY
             k=J+(i-1)*NEY+NEZ1
             write(44,*) NoOP(k,2), NOoP(k,3), NoOP(k,6), NOoP(k,5)  
             write(44,*) NoOP(k,5), NOoP(k,6), NOoP(k,9), NOoP(k,8) 
           END Do 
        END DO 
c	IF(TIME.eq.STEPT)THEN
c	WRITE(44,1009)TIME,2*NEY+1,2*NEX+1
c	ELSE 
c	WRITE(44,1011)TIME,2*NEY+1,2*NEX+1
c	END IF
c        DO I=1,NHADD
c        IF(MDF(I) .EQ. 2)THEN
c        UW=SOL(NOPP(I))
c        VW=SOL(NOPP(I)+1)
c        WRITE(44,1006)XPT(I),YPT(I),UW,VW
c        END IF
c        IF(MDF(I) .EQ. 3)THEN
c        UW=SOL(NOPP(I))
c        VW=SOL(NOPP(I)+1)
c        WRITE(44,1006)XPT(I),YPT(I),UW,VW
c        END IF
c        END DO

        IF(TIME.eq.STEPT)THEN
	WRITE(444,1010)TIME,NEYZ1+1,NEXZ1+1
	ELSE
	WRITE(444,1012)TIME,NEYZ1+1,NEXZ1+1
	END IF
        DO I=1,(2*NEXZ1+1)*(2*NEYZ1+2)
        IF(MDF(I) .EQ. 3)THEN
        PW=SOL(NOPP(I)+2)
        WRITE(444,1007)XPT(I),YPT(I),PW
        END IF
        END DO

        IF(TIME.eq.STEPT)THEN
	WRITE(555,1010)TIME,NEY+1,NEX+1
	ELSE
	WRITE(555,1012)TIME,NEY+1,NEX+1
	END IF
        DO I=1,(NEYZ1*2+2)*2*NEXZ1+1,NEYZ1*2+2
        IF(MDF(I) .EQ. 3)THEN
        PW=SOL(NOPP(I)+2)
        WRITE(555,1007)XPT(I),YPT(I),PW
        END IF
        END DO

        DO I=(2*NEYZ1+2)*2*NEXZ1+2,NHADD
        IF(MDF(I) .EQ. 3)THEN
        PW=SOL(NOPP(I)+2)
        WRITE(555,1007)XPT(I),YPT(I),PW
        END IF
        END DO
1001    FORMAT(//,10X,' Re =',E20.10,' We =',E20.10,' F =',E20.10,
     1                ' OME =',E20.10)
1002    FORMAT(//,10X,' MESH DATA',//)
1003    FORMAT(10X,' NEY =',I5,' NEX =',I5)
1004    FORMAT(10X,' Dimensionless length (x-direction) =',E20.10)
1005    FORMAT(//,10X,' SOLUCION FOR TIME =',E20.10,' ERROR =',
     1   E20.10,' ITE =',I5,//)
1006    FORMAT(10X,E20.10,E20.10,
     1   E20.10,E20.10)
1007    FORMAT(10X,E15.7,E15.7,E20.10)
1008    FORMAT(10X,' A0:',E20.10,' OMEGA:',E15.6)
1009    FORMAT('title = "velocity" ',/,'vari
     1ables = "x", "y","u","v"',/'zone T="',E20.10,'",
     2N=', I5, ',E=',I5,',DATAPACKING=POINT
     3,ZONETYPE=FEQUADRILATERAL')
c1009    FORMAT('title = "velocity" ',/,'vari
c     1ables = "x", "y","u","v"',/'zone T="',E20.10,'",
c     2i=', I5, ',j=',I5,',DATAPACKING=POINT')
1010    FORMAT('title = "pressure" ',/,'variables = "
     1x", "y","p"',/'zone T="',E20.10,'",i=', I5, ',
     2j=',I5,',DATAPACKING=POINT')
1011    FORMAT('zone T="',E20.10,'",
     1N=',I5, ',E=',I5,',DATAPACKING=POINT
     2,ZONETYPE=FEQUADRILATERAL')
c 1011    FORMAT('zone T="',E20.10,'",
c     1i=',I5, ',j=',I5,',DATAPACKING=POINT')
1012    FORMAT('zone T="',E20.10,'",
     1i=',I5, ',j=',I5,',DATAPACKING=POINT')
        RETURN
	END

C       *****************************************

C       *****************************************
	SUBROUTINE ABFIND

c	Last modification 03-Oct-2012-JK

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION W(3),GP(3)
	DIMENSION DTANX(3),DTANY(3),DPY(3),RHS(25) ! DTANX, DTANY were previous auxilary variable, not needed now
	DIMENSION V1(9),V2(9),V3(4)
        DIMENSION UPYH(3),VPYH(3),PIY(9)
        DIMENSION astemp1(25), astemp2(25) 
	COMMON/RES/R2(15999)
	COMMON/TFUN1/PHI(9),PSI(4),PHIC(9),PHIE(9),PHICC(9)
	COMMON/DER1/AXPT(9),AYPT(9),WYL(3),DIRL(3,2),WXL(3)
	COMMON/DER2/U,UC,UE,V,VC,VE,P,rs !Added on 02-Feb-2012 By: JK To evaulate the derivative of the primary variables for the Viscosity Calculation
	COMMON/DER3/XE,XC,YC,YE,DET,Y
        COMMON/MESPA/XW,YW
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/DERJA/DYEDH(3),DYCDH(3),DH(3),DHE(3),DHC(3) ! DYEDH DYCDH were previous auxilary variable, not needed now
	COMMON/PARAM/RE,F,OME,WE
	COMMON/FSLAB/IELFS(1700), IFZONE(1700)
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/MESH1/NNY,NNX,NNX1,NNX2
	COMMON/MESCO/YPT(7011),XPT(7011)
	COMMON/MESHCON/NEY,NEX,NEX1,NEX2
	COMMON/MESHWT/WY(41),WX(180)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
	COMMON/FRON1/NELL
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
	COMMON/FRON3/A(25,25)
	COMMON/TEMP1/UV,VV,UPV,VPV,DYT1,DYT2,DXT1,DXT2,DYH(3)
	COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME
	COMMON/NMESH/DIR(180,2)	!Added on Jun-2010 By:JK To store the sin and cos of each spine.
        COMMON/NONNEWPARA/alp,beta,pn !Added on Jun-2010 By:JK Parameters for the non-Newtonian fluid
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1 !JKMoving
        COMMON/ZMESH/wxz(15) !JKMoving
        COMMON/DIRTEMP/dirtemp(3,2),wyltemp(3)!JKMoving 
	INTEGER W1,W2,W3,W4,V1,V2,V3
	DATA W,GP/.27777777778,.44444444444,.27777777778,
     1         .1127016654,.5,.8872983346/
	DATA V1,V2,V3/1,4,6,9,11,13,15,18,20,  ! V1 defines the local equation number for U velocity
     1             2,5,7,10,12,14,16,19,21,    ! V2 defines the local equation number for V velocity
     1             3,8,17,22/                  ! V3 defines the local equation number for the Pressure. Last three equations are assigned to the H value on the interface.
C       *******************************
        if(nell .eq. 1)then
          vol=0.d0	! ??var
          ekin=0.d0	! ??var
          esup=0.d0	! ??var
        end if
c        np=nopp(nop((ney*nex),12))       ! 
C     Local coordinates
	DO I=1,9
	  J=NOP(NELL,I)
	  AXPT(I)=XPT(J)
	  AYPT(I)=YPT(J)
	END DO

C     Local Proportionality for coordinates within element
       IF(IFZONE(NELL).EQ.1)THEN
        Do i=1,3
          wyl(i)=0.
          wxl(i)=0.
        End Do
        icol=(nell-1)/neyz1
        irowp=nell-neyz1*icol

        Do i=1,3
          wxl(i)=1
          wyl(i)=wy(2*(irowp-1)+i)        ! Local Proportionality 
          dirl(i,1)=dir((icol*2)+i,1)     ! Cos for the local spine with respect to z coor. 
          dirl(i,2)=dir((icol*2)+i,2)     ! Sin for the local spine with respect to z coor.
        End Do

C     Local Residual Vector and Jacobian Matrix
	DO J=1,25
	  RHS(J)=0.!Nullify Local Residual Vector
	  DO I=1,25
	    A(I,J)=0.!Nullify local Jacobian Matrix
	  END DO
	END DO

C     Assigning Global Equation Number to Local Numbering Vector
	DO I=1,9
	 W1(I)=NOPP(NOP(NELL,I)) !W1(i)--U component 
	 W2(I)=W1(I)+1 !W2(i)--V component 
	END DO
	 W3(1)=W1(1)+2 !W3(i)--P component
	 W3(2)=W1(3)+2
	 W3(3)=W1(7)+2
	 W3(4)=W1(9)+2
	NELH=((NELL-1)/NEYZ1+1)*NEYZ1 ! NELH is the local var and take Integer
	DO I=1,3
	 W4(I)=NOPP(NOP(NELH,9+I)) !W4(i)--h component
	END DO
      END IF
       IF(IFZONE(NELL).EQ.2)THEN
        Do i=1,3
          wyl(i)=0.
          wxl(i)=0.
        End Do
        icol=(nell-1-NEZ1)/ney
        irowp=(nell-NEZ1)-ney*icol

        Do i=1,3
          wxl(i)=(wx((icol*2)+i)-wx(NNX1))/(wx(1)
     1     -wx(NNX1))
          wyl(i)=wxz(2*(irow-1)+i)     ! Local Proportionality 
c          dirl(i,1)=dir((icol*2)+i+NNXZ1-1,1)     ! Cos for the local spine with respect to z coor. 
c          dirl(i,2)=dir((icol*2)+i+NNXZ1-1,2)     ! Sin for the local spine with respect to z coor.
c          dirl(i,1)=-1.    ! Cos for the local spine with respect to z coor. 
c          dirl(i,2)=0.    ! Sin for the local spine with respect to z coor.
          dirl(i,1)=0.    ! Cos for the local spine with respect to z coor. 
          dirl(i,2)=1.    ! Sin for the local spine with respect to z coor.
        End Do

C     Local Residual Vector and Jacobian Matrix
	DO J=1,25
	  RHS(J)=0.!Nullify Local Residual Vector
	  DO I=1,25
	    A(I,J)=0.!Nullify local Jacobian Matrix
	  END DO
	END DO

C     Assigning Global Equation Number to Local Numbering Vector
	DO I=1,9
	 W1(I)=NOPP(NOP(NELL,I)) !W1(i)--U component 
	 W2(I)=W1(I)+1 !W2(i)--V component 
	END DO
	 W3(1)=W1(1)+2 !W3(i)--P component
	 W3(2)=W1(3)+2
	 W3(3)=W1(7)+2
	 W3(4)=W1(9)+2
	NELH=((NELL-NEZ1-1)/NEY+1)*NEY ! NELH is the local var and take Integer
	DO I=1,3
	 W4(I)=NOPP(nhadd) !W4(i)--h component
	END DO
      END IF
       IF(IFZONE(NELL).EQ.3)THEN
        Do i=1,3
          wyl(i)=0.
          wxl(i)=0.
        End Do
        icol=(nell-1-NEZ1)/ney
        irowp=(nell-NEZ1)-ney*icol

        Do i=1,3
          wxl(i)=(wx((icol*2)+i)-wx(NNX1))/(wx(1)
     1     -wx(NNX1))
          wyl(i)=wy(2*(irowp-1-NEXZ1)+i)        ! Local Proportionality 
          dirl(i,1)=dir((icol*2)+i+NNXZ1-1,1)     ! Cos for the local spine with respect to z coor. 
          dirl(i,2)=dir((icol*2)+i+NNXZ1-1,2)     ! Sin for the local spine with respect to z coor.
        End Do

C     Local Residual Vector and Jacobian Matrix
	DO J=1,25
	  RHS(J)=0.!Nullify Local Residual Vector
	  DO I=1,25
	    A(I,J)=0.!Nullify local Jacobian Matrix
	  END DO
	END DO

C     Assigning Global Equation Number to Local Numbering Vector
	DO I=1,9
	 W1(I)=NOPP(NOP(NELL,I)) !W1(i)--U component 
	 W2(I)=W1(I)+1 !W2(i)--V component 
	END DO
	 W3(1)=W1(1)+2 !W3(i)--P component
	 W3(2)=W1(3)+2
	 W3(3)=W1(7)+2
	 W3(4)=W1(9)+2
	NELH=((NELL-NEZ1-1)/NEY+1)*NEY ! NELH is the local var and take Integer
	DO I=1,3
	 W4(I)=NOPP(NOP(NELH+NEZ1,9+I)) !W4(i)--h component
	END DO
      END IF
C     Numerical integration (Gauss). Surface integrals.
	DO L=1,3
	DO M=1,3

          CALL TFUNCT(GP(L),GP(M)) ! Evaluate the shape function on Gaussian Points
          CALL DERIV		! Evaluate the primary variables and their derivatives and SOME auxiliary variables

C       auxiliary variables for time derivatives and the Determination of the Jacobian
    
          dyt=ct1*dyt1/stept+ct2*dyt2
	  dxt=ct1*dxt1/stept+ct2*dxt2 !Added Jun-2010 By: JK To account the movement for the spines in z direction
          upn=(u-uv)*ct1/stept+ct2*upv
          vpn=ct1*(v-vv)/stept+ct2*vpv
          tdt=Xc*Ye-Xe*Yc            !Jacobian for coordinate transformation
          Do i=1,9
	     piy(i)=-ye*phic(i)+yc*phie(i) ! Evaluate the phy(i) ?? May not be needed
          End Do
	  CALL DERJAC ! Evaluate SOME auxiliary variables
          
         
          c1=w(l)*w(m) !constant for 2-D guassian integration 

          strainrate=        (Ve**2*Y**2*(2*Xc**2 + Yc**2) 
     *      + Ue**2*Y**2*(Xc**2 + 2*Yc**2) + 
     *    Xe**2*((Uc**2 + 2*Vc**2)*Y**2 + 2*V**2*Yc**2) - 
     *    2*Xe*(Uc*Vc*Y**2 + 2*V**2*Xc*Yc)*Ye + 
     *    (2*V**2*Xc**2 + (2*Uc**2 + Vc**2)*Y**2)*Ye**2 - 
     *    2*Ue*Y**2*(Uc*Xc*Xe + Ve*Xc*Yc - Vc*Xc*Ye + 2*Uc*Yc*Ye) - 
     *    2*Ve*Y**2*(2*Vc*Xc*Xe - Uc*Xe*Yc + Vc*Yc*Ye))/
     *  (Y**2*(Xe*Yc - Xc*Ye)**2)

          eta= beta+(1-beta)*((1+alp*alp*strainrate)**((pn-1.d0)/2.d0))
  
          etap= ((pn-1.d0)/2.d0)*(1-beta)*(alp*alp)*((1
     *          +alp*alp*strainrate)**((pn-3.d0)/2))!eta prime with respcet to strainrate 
  
          ekin=ekin+c1*(u**2+v**2)*xc*ye      ! Kinetic Energy MAY not be needed??
c       Momentum balance (x and y components)
	DO I=1,9
	  K=V1(I)
	  N=V2(I)
c global auxilary variable

c        write(*,*)etap
c for x momentum          
          AA=RE*upn*phi(i)
          B1=RE*Uc*(dyt - V)*phi(i)
          B2=RE*Ue*(-dyt + V)*phi(i)
          B3=(RE*(-dxt + U)*Uc*We*phi(i) - P*phic(i))/We
          B4=RE*(dxt - U)*Ue*phi(i) + (P*phie(i))/We
          CC1=-(Ve*Xc*Yc) + Ue*(Xc**2 + 2*Yc**2) + Vc*Xc*Ye 
     *       - Uc*(Xc*Xe + 2*Yc*Ye)
          C2=Xe*(-(Ue*Xc) + Uc*Xe + Ve*Yc) 
     *       - (Vc*Xe + 2*Ue*Yc)*Ye + 2*Uc*Ye**2
          B=Xe*B1+Xc*B2+Ye*B3+Yc*B4
          C=phie(i)*CC1+phic(i)*C2
       
c for y momentum
          E=RE*vpn*phi(i)
          F1=RE*(dyt - V)*Vc*phi(i) + (P*phic(i))/We
          F2=(RE*(-dyt + V)*Ve*We*phi(i) - P*phie(i))/We
          F3=RE*(-dxt + U)*Vc*phi(i)
          F4=RE*(dxt - U)*Ve*phi(i)
          G1=-(Ue*Xc*Yc) + Uc*Xe*Yc + Ve*(2*Xc**2 + Yc**2) 
     *     - Vc*(2*Xc*Xe + Yc*Ye)
          G2=2*Vc*Xe**2 + Ue*Xc*Ye - Uc*Xe*Ye + Vc*Ye**2 
     *     - Ve*(2*Xc*Xe + Yc*Ye)
          PMR=-((P*phi(i))/We)
          TMR=(2*V*phi(i))/Y
          F=Xe*F1+Xc*F2+Ye*F3+Yc*F4
          G=phie(i)*G1+phic(i)*G2
	  
          DO J=1,9
	    KK=V1(J)
	    NN=V2(J)
c eta primes with respect to Uj, Vj
            etaup=etap*(0.0*phi(j)
     *            +((2*Xe*(-(Ue*Xc) + Uc*Xe + Ve*Yc) 
     *            - 2*(Vc*Xe + 2*Ue*Yc)*Ye + 4*Uc*Ye**2)/
     *            (Xe*Yc - Xc*Ye)**2)*phic(j)
     *            +((2*(-(Ve*Xc*Yc) + Ue*(Xc**2 + 2*Yc**2) 
     *            + Vc*Xc*Ye - Uc*(Xc*Xe + 2*Yc*Ye)))/
     *            (Xe*Yc - Xc*Ye)**2)*phie(j)) !d/dU*dU/dUj+d/dUc*dUc/dUj+dUe/dUj 

            etavp=etap*(((4*V)/Y**2)*phi(j)+((2*(-(Ue*Xc*Yc) 
     *           + Uc*Xe*Yc + Ve*(2*Xc**2 + Yc**2) - Vc*(2*Xc*Xe 
     *           + Yc*Ye)))/(Xe*Yc - Xc*Ye)**2)*phie(j)
     *           +((2*(2*Vc*Xe**2 + Ue*Xc*Ye - Uc*Xe*Ye 
     *           + Vc*Ye**2 - Ve*(2*Xc*Xe + Yc*Ye)))/(Xe*Yc 
     *           - Xc*Ye)**2)*phic(j)) !d/dV*dV/dVj+d/dVc*dVc/dVj+dVe/dVj 

c COLMz derivative with respect to Uj
            tdtup=0.0
            Aup=(ct1*RE*phi(i)*phi(j))/stept
            B1up=-(RE*(-dyt + V)*phi(i)*phic(j))
            B2up=RE*(-dyt + V)*phi(i)*phie(j)
            B3up=RE*Uc*phi(i)*phi(j) + RE*(-dxt 
     *      + U)*phi(i)*phic(j)
            B4up=-(RE*Ue*phi(i)*phi(j)) - RE*(-dxt 
     *      + U)*phi(i)*phie(j)
            Xeup=0.0
            Xcup=0.0
            Yeup=0.0
            Ycup=0.0
            Bup=Xe*B1up+Xc*B2up+Ye*B3up+Yc*B4up
            C1up=-(Xc*Xe*phic(j)) - 2*Yc*Ye*phic(j) 
     *      + (Xc**2 + 2*Yc**2)*phie(j)
            C2up=2*Ye**2*phic(j) - 2*Yc*Ye*phie(j) 
     *      + Xe*(Xe*phic(j) - Xc*phie(j))
            Cup=phie(i)*C1up+phic(i)*C2up

cCOLMr derivative with respect to Uj
             Eup=0.0
             F1up=0.0
             F2up=0.0
             F3up=RE*Vc*phi(i)*phi(j)
             F4up=-(RE*Ve*phi(i)*phi(j))
             G1up=Xe*Yc*phic(j) - Xc*Yc*phie(j)
             G2up=-(Xe*Ye*phic(j)) + Xc*Ye*phie(j)
             Fup=Xe*F1up+Xc*F2up+Ye*F3up+Yc*F4up
             Gup=phie(i)*G1up+phic(i)*G2up
             PMRup=0.0
             TMRup=0.0        

cCOLMz derivative with respect to Vj

            Avp=0.0
            B1vp=-(RE*Uc*phi(i)*phi(j))
            B2vp=RE*Ue*phi(i)*phi(j)
            B3vp=0.0
            B4vp=0.0
            Xevp=0.0
            Xcvp=0.0
            Yevp=0.0
            Ycvp=0.0
            Bvp=Xe*B1vp+Xc*B2vp+Ye*B3vp+Yc*B4vp
            C1vp=Xc*Ye*phic(j) - Xc*Yc*phie(j)
            C2vp=-(Xe*Ye*phic(j)) + Xe*Yc*phie(j)
            Cvp=phie(i)*C1vp+phic(i)*C2vp


cCOLMr derivative with respect to Vj
             Evp=(ct1*RE*phi(i)*phi(j))/stept
             F1vp=-(RE*Vc*phi(i)*phi(j)) - RE*(-dyt + V)*phi(i)*phic(j)
             F2vp=RE*Ve*phi(i)*phi(j) + RE*(-dyt + V)*phi(i)*phie(j)
             F3vp=RE*(-dxt + U)*phi(i)*phic(j)
             F4vp=-(RE*(-dxt + U)*phi(i)*phie(j))
             G1vp=-2*Xc*Xe*phic(j) - Yc*Ye*phic(j) 
     *        + (2*Xc**2 + Yc**2)*phie(j)
             G2vp=2*Xe**2*phic(j) + Ye**2*phic(j) 
     *        - (2*Xc*Xe + Yc*Ye)*phie(j)
             Fvp=Xe*F1vp+Xc*F2vp+Ye*F3vp+Yc*F4vp
             Gvp=phie(i)*G1vp+phic(i)*G2vp
             PMRvp=0.0
             TMRvp=(2*phi(i)*phi(j))/Y
c             write(*,*)etaup

c Jacobian for COLM _ z _ r

	    a(k,kk)=a(k,kk)-c1*((Aup*tdt+Bup+etaup*C/tdt+eta*Cup/tdt)*Y)
	    a(n,kk)=a(n,kk)-c1*((Eup*tdt+Fup+etaup*G/tdt+eta*Gup/tdt)*Y
     *              +tdt*(PMRup+etaup*TMR+eta*TMRup))
	    a(k,nn)=a(k,nn)-c1*((Avp*tdt+Bvp+etavp*C/tdt+eta*Cvp/tdt)*Y)
	    a(n,nn)=a(n,nn)-c1*((Evp*tdt+Fvp+etavp*G/tdt+eta*Gvp/tdt)*Y
     *              +tdt*(PMRvp+etavp*TMR+eta*TMRvp))


  	    	  
	  END DO
c COLMz&r derivative with respect to Pj
	  DO J=1,4
	    KK=V3(J)

            a(k,kk)=a(k,kk)-c1*((-(Ye*phic(i)) 
     *      + Yc*phie(i))/We)*psi(j)*Y
            a(n,kk)=a(n,kk)-c1*(((Xe*phic(i) - Xc*phie(i))/We)*Y
     *            +(-(phi(i)/We))*tdt)*psi(j)
	  END DO
c Residual for COLMz &r
	   rhs(k)=rhs(k)+c1*(AA*tdt+B+eta*C/tdt)*Y
	   rhs(n)=rhs(n)+c1*((E*tdt+F +eta*G/tdt)*Y+(pmr+eta*tmr)*tdt)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c COLMz&r derivative with respect to Hj
	  DO J=1,3
	    JK=22+J

              Xh=dh(j)*dirl(j,1)
              Yh=dh(j)*dirl(j,2)  !dY/dh Y=y(j)dh(j)dirl(j,2) j=1,3
              Xeh=dhe(j)*dirl(j,1)
              Xch=dhc(j)*dirl(j,1)
              Yeh=dhe(j)*dirl(j,2)
              Ych=dhc(j)*dirl(j,2)

              etady=(-4*V**2)/Y**3
              etadye= (2*(Ve**2*Xc*(2*Xc**2 + Yc**2) 
     -                + Ue**2*(Xc**3 + 2*Xc*Yc**2) + 
     -                Xe*(-(Uc*Vc*(Xe*Yc + Xc*Ye)) 
     -                + Vc**2*(2*Xc*Xe + Yc*Ye) + 
     -                Uc**2*(Xc*Xe + 2*Yc*Ye)) - 
     -                Ve*(-2*Uc*Xc*Xe*Yc + Vc*(4*Xc**2*Xe 
     -                + Xe*Yc**2 + Xc*Yc*Ye)) + 
     -                Ue*(Xc*(-2*Ve*Xc*Yc + Vc*Xe*Yc 
     -                + Vc*Xc*Ye) -2*Uc*(Xe*(Xc**2 + Yc**2) 
     -                + Xc*Yc*Ye))))/(Xe*Yc - Xc*Ye)**3
              etadyc=(-2*(Ve**2*Xc*(2*Xc*Xe + Yc*Ye) 
     -              + Ue**2*Xc*(Xc*Xe + 2*Yc*Ye) + 
     -              Xe*(-2*Uc*Vc*Xe*Ye + Vc**2*(2*Xe**2 + Ye**2) + 
     -              Uc**2*(Xe**2 + 2*Ye**2)) - 
     -              Ue*(Xc*(Ve*Xe*Yc + Ve*Xc*Ye - 2*Vc*Xe*Ye) + 
     -              2*Uc*(Xe*Yc*Ye + Xc*(Xe**2 + Ye**2))) + 
     -              Ve*(Uc*Xe*(Xe*Yc + Xc*Ye) - Vc*(Xe*Yc*Ye 
     -              + Xc*(4*Xe**2 + Ye**2)))))/(Xe*Yc - Xc*Ye)**3
              etadx=0.0
              etadxe= (-2*(Ue**2*Yc*(Xc**2 + 2*Yc**2) 
     -              + Ve**2*(2*Xc**2*Yc + Yc**3) + 
     -              Ye*(-(Uc*Vc*(Xe*Yc + Xc*Ye)) 
     -              + Vc**2*(2*Xc*Xe + Yc*Ye) + 
     -              Uc**2*(Xc*Xe + 2*Yc*Ye)) - 
     -              Ue*(2*Xc*Yc*(Ve*Yc - Vc*Ye) + 
     -              Uc*(Xc*Xe*Yc + Xc**2*Ye + 4*Yc**2*Ye)) + 
     -              Ve*(Uc*Yc*(Xe*Yc + Xc*Ye) - 2*Vc*(Xc*Xe*Yc 
     -              + (Xc**2 + Yc**2)*Ye))))/(Xe*Yc - Xc*Ye)**3
              etadxc=(2*(2*Ve*Xe*(Ve*Xc - Vc*Xe)*Yc + 
     -              (-2*Vc*Ve*Xc*Xe + 2*Vc**2*Xe**2 
     -              + (Uc*Xe + Ve*Yc)**2)*Ye - 
     -              2*Vc*(Uc*Xe + Ve*Yc)*Ye**2 
     -              + (2*Uc**2 + Vc**2)*Ye**3 + 
     -              Ue**2*Yc*(Xc*Xe + 2*Yc*Ye) - 
     -              Ue*((Ve*Yc - Vc*Ye)*(Xe*Yc + Xc*Ye) + 
     -              Uc*(Xe**2*Yc + Xc*Xe*Ye 
     -              + 4*Yc*Ye**2))))/(Xe*Yc - Xc*Ye)**3
              etah=etap*(etady*Yh+etadye*Yeh+etadyc*Ych
     *                   etadx*Xh+etadxe*Xeh+etadxc*Ych)
 
              Ah=0.0

              B1h=((ct1*RE*Uc*phi(i))/stept)*Yh
              B2h=-((ct1*RE*Ue*phi(i))/stept)*Yh
              B3h=-((ct1*RE*Uc*phi(i))/stept)*Xh
              B4h=(ct1*RE*Ue*phi(i))/stept*Xh
              Bh=Xeh*B1+B1h*Xe+Xch*B2+B2h*Xc
     *        +Yeh*B3+Ye*B3h+Ych*B4+Yc*B4h


              C1h=(Vc*Xc - 2*Uc*Yc)*Yeh 
     *           + (-(Ve*Xc) + 4*Ue*Yc - 2*Uc*Ye)*Ych 
     *           +(-(Uc*Xc))*Xeh
     *           +(2*Ue*Xc - Uc*Xe - Ve*Yc + Vc*Ye)*Xch
              C2h=(-(Vc*Xe) - 2*Ue*Yc + 4*Uc*Ye)*Yeh 
     *           + (Ve*Xe - 2*Ue*Ye)*Ych 
     *           +(-(Ue*Xc) + 2*Uc*Xe + Ve*Yc - Vc*Ye)*Xeh
     *           +(-(Ue*Xe))*Xch
              Ch=phie(i)*C1h+phic(i)*C2h

             Eh=0.0

             F1h=((ct1*RE*Vc*phi(i))/stept)*Yh
             F2h=(-((ct1*RE*Ve*phi(i))/stept))*Yh
             F3h=(-((ct1*RE*Vc*phi(i))/stept))*Xh
             F4h=((ct1*RE*Ve*phi(i))/stept)*Xh
             Fh=Xeh*F1+F1h*Xe+Xch*F2+F2h*Xc+Yeh*F3+Ye*F3h+Ych*F4+Yc*F4h

             G1h=(-(Vc*Yc))*Yeh 
     *           + (-(Ue*Xc) + Uc*Xe + 2*Ve*Yc - Vc*Ye)*Ych 
     *           +(-2*Vc*Xc + Uc*Yc)*Xeh
     *           +(4*Ve*Xc - 2*Vc*Xe - Ue*Yc)*Xch
             G2h=(Ue*Xc - Uc*Xe - Ve*Yc + 2*Vc*Ye)*Yeh
     *           +(-(Ve*Ye))*Ych
     *           +(-2*Ve*Xc + 4*Vc*Xe - Uc*Ye)*Xeh
     *           +(-2*Ve*Xe + Ue*Ye)*Xch
             Gh=phie(i)*G1h+phic(i)*G2h


             PMRh=0.0
             TMRh=((-2*V*phi(i))/Y**2)*Yh
              
             tdth=Xc*Yeh-Xe*Ych-Yc*Xeh+Ye*Xch

             a(k,jk)= a(k,jk)-c1*((Ah*tdt+AA*tdth+Bh+etah*C/tdt
     *               +eta*Ch/tdt-eta*C*tdth/(tdt*tdt))*Y
     *               +Yh*(AA*tdt+B+eta*C/tdt))

c              a(n,jk)=a(n,jk)-c1*((E*tdt+F +eta*G/tdt)*Yh
c     *               +Y*(Eh*tdt+E*tdth+Fh+etah*G/tdt+eta*Gh/tdt
c     *               -eta*G*tdth/(tdt*tdt))+(pmr+eta*tmr)*tdth*(pmrh
c     *               +eta*tmrh*etah*tmr)*tdt)
              a(n,jk)=a(n,jk)-c1*((E*tdt+F +eta*G/tdt)*Yh
     *               +Y*(Eh*tdt+E*tdth+Fh+etah*G/tdt+eta*Gh/tdt
     *               -eta*G*tdth/(tdt*tdt))+(PMR+eta*TMR)*tdth
     *               +tdt*(PMRh+etah*TMR+eta*TMRh))
	  END DO

	END DO
C       Mass balance (Continuity equation)
	DO I=1,4
	  K=V3(I)
	  DO J=1,9
	    KK=V1(J)
	    NN=V2(J)

            a(k,kk)=a(k,kk)-(c1*Y*(Ye*phic(j) - Yc*phie(j))*psi(i))     !JK_Full_Expression

			a(k,nn)=a(k,nn)-(-(c1*det*Y*(phi(j)/Y 
     1             - (-(Xe*phic(j)) + Xc*phie(j))/det)*psi(i)))         !JK_Full_Expression 
	  END DO
	  DO J=1,3
	    JK=J+22

           a(k,jk)=a(k,jk)-(c1*psi(i)*((Ve*dh(j) 
     *             + V*dhe(j))*dirl(j,2)*Xc - (Vc*dh(j) 
     *             + V*dhc(j))*dirl(j,2)*Xe 
     *             - Ue*dirl(j,2)*Yc*phi(j)*wyl(1) 
     *             + Uc*dirl(j,2)*Ye*phi(j)*wyl(1) 
     *             + Ve*dirl(j,1)*Y*phic(j)*wyl(1) 
     *             - Ue*dirl(j,2)*Y*phic(j)*wyl(1) 
     *             + V*dirl(j,1)*Ye*phic(j)*wyl(1) 
     *             - Vc*dirl(j,1)*Y*phie(j)*wyl(1) 
     *             + Uc*dirl(j,2)*Y*phie(j)*wyl(1) 
     *             - V*dirl(j,1)*Yc*phie(j)*wyl(1) 
     *             - Ue*dirl(j,2)*Yc*phi(3 + j)*wyl(2) 
     *             + Uc*dirl(j,2)*Ye*phi(3 + j)*wyl(2) 
     *             + Ve*dirl(j,1)*Y*phic(3 + j)*wyl(2) 
     *             - Ue*dirl(j,2)*Y*phic(3 + j)*wyl(2) 
     *             + V*dirl(j,1)*Ye*phic(3 + j)*wyl(2) 
     *             - Vc*dirl(j,1)*Y*phie(3 + j)*wyl(2) 
     *             + Uc*dirl(j,2)*Y*phie(3 + j)*wyl(2) 
     *             - V*dirl(j,1)*Yc*phie(3 + j)*wyl(2) 
     *             + (Ye*(Uc*dirl(j,2)*phi(6 + j) 
     *             + V*dirl(j,1)*phic(6 + j)) 
     *             - Yc*(Ue*dirl(j,2)*phi(6 + j) 
     *             + V*dirl(j,1)*phie(6 + j)) 
     *             + Y*(Ve*dirl(j,1)*phic(6 + j) 
     *             - Vc*dirl(j,1)*phie(6 + j) 
     *             + dirl(j,2)*(-(Ue*phic(6 + j)) 
     *             + Uc*phie(6 + j))))*wyl(3)))                  !JK_Full_Expression
	  END DO

	  rhs(k)=rhs(k)+(-(c1*(det*V + Y*(-(Ve*Xc) + Vc*Xe 
     1          + Ue*Yc - Uc*Ye))*psi(i)))  ! JK_det_test
	END DO
	
	END DO ! fin <M>
	END DO ! fin <L>
c	IF(Nell.eq.ney)then
c	rstemp=rhs(V2(7))
c	astemp1=a(V2(7),23)
c	astemp2=a(V2(7),24)
c	astemp3=a(V2(7),25)
c        rstemp2=rhs(V1(7))
c        astemp4=a(V1(7),23)
c        astemp5=a(V1(7),24)
c        astemp6=a(V1(7),25)
c	End IF
c         IF(Nell.eq.ney)then
c           rstemp=rhs(V2(7))
c           Do i=1,25
c             astemp1(i)=a(V2(7),i)
c             astemp2(i)=a(V1(7),i)
c           End Do
c           rstemp2=rhs(V1(7))
c	 End IF
C     Line integrals
	IF(IELFS(NELL) .EQ. 1)THEN 

	DO K=1,3
	  CALL TFUNCT(GP(K),1.d0)
	  CALL DERIV
	  CALL DERJAC
          com=1./dsqrt(xc**2+yc**2)
          com3=com**3
          tanx=xc*com
          tany=yc*com
          dyt=ct1*dyt1/stept+ct2*dyt2
	  dxt=ct1*dxt1/stept+ct2*dxt2
          esup=esup+w(k)*dsqrt(xc**2+yc**2)
	  DO J=1,3
            dtanx(j)=-xc*yc*dycdh(j)*com3
            dtany(j)=xc**2*dycdh(j)*com3
	  END DO
	  DO I=7,9
            m=v1(i)
            n=v2(i)
            rhs(m)=rhs(m)+((Xc*Y*phic(i)*w(k))/(We*dsqrt(Xc**2 
     1            + Yc**2)))                      !JK_Full_Expression
	     rhs(n)=rhs(n)+((((Xc**2 + Yc**2)*phi(i) 
     1             + Y*Yc*phic(i))*w(k))/(We*dsqrt(Xc**2 
     2             + Yc**2)))              !JK_Full_Expression
	    DO J=1,3
	      JK=J+22
              a(m,jk)=a(m,jk)-(((dh(j)*dirl(j,2)*Xc**3 
     1               + dhc(j)*dirl(j,1)*Y*Yc**2 
     2               + dirl(j,2)*Xc*Yc*(-(dhc(j)*Y) 
     3               + dh(j)*Yc))*phic(i)*w(k))/(We*(Xc**2 
     4               + Yc**2)**1.5))                  !JK_Full_Expression

              a(n,jk)=a(n,jk)-(((dhc(j)*dirl(j,1)*Xc**3*phi(i) 
     1                + dirl(j,2)*Yc**3*(dhc(j)*phi(i) 
     2                + dh(j)*phic(i)) 
     3                + dhc(j)*dirl(j,1)*Xc*Yc*(Yc*phi(i) 
     4                - Y*phic(i)) 
     5                + dirl(j,2)*Xc**2*(dhc(j)*Y*phic(i) 
     6                + Yc*(dhc(j)*phi(i) 
     7                + dh(j)*phic(i))))*w(k))/(We*(Xc**2 
     8                + Yc**2)**1.5))
	    END DO
	  END DO
        
c         IF(Nell.eq.ney)then
c           rhs(V2(7))=rstemp
c           Do i=1,25
c            a(V2(7),i)= astemp1(i)
c            a(V1(7),i)= astemp2(i)
c           End Do
c         rhs(V1(7))= rstemp2
c	 End IF
c        IF(Nell.eq.ney)then
c       rhs(V2(7))=rstemp
c       a(V2(7),23)=astemp1
c       a(V2(7),24)=astemp2
c       a(V2(7),25)=astemp3
c        rhs(V1(7))=rstemp2
c        a(V1(7),23)=astemp4
c        a(V1(7),24)=astemp5
c        a(V1(7),25)=astemp6
c        End IF
C         Kinematic equation
	  DO I=1,3
	    IK=22+I
	    DO J=1,3
	      JK=J+22
         ypr=0.0d0
         ytpr=0.0d0
         xtpr=0.0d0
         xcpr=0.0d0
         ycpr=0.0d0
         xyxy=0.0d0
          ypr=dirl(j,2)*(phi(j)*wyl(1)
     1    + phi(3 + j)*wyl(2) + phi(6 + j)*wyl(3))
          ytpr=(ct1*dirl(j,2)*(phi(j)*wyl(1)
     1    + phi(3 + j)*wyl(2) + phi(6 + j)*wyl(3)))/stept
          xtpr=(ct1*dirl(j,1)*(phi(j)*wyl(1)
     1    + phi(3 + j)*wyl(2) + phi(6 + j)*wyl(3)))/stept
          ycpr=dirl(j,2)*(phic(j)*wyl(1)
     1    + phic(3 + j)*wyl(2) + phic(6 + j)*wyl(3))
          xcpr=dirl(j,1)*(phic(j)*wyl(1)
     1    + phic(3 + j)*wyl(2) + phic(6 + j)*wyl(3))
          xyxy=ypr*((-(dyt*Xc) + V*Xc + (dxt- U)*Yc))
     1    +Y*(ycpr*dxt+xtpr*Yc-xcpr*dyt-ytpr*Xc+Xcpr*V-Ycpr*U)
          a(ik,jk)=a(ik,jk)-(w(k)*phi(6+i)*xyxy)


	    END DO
	    DO J=7,9
	      KK=V1(J)
	      NN=V2(J)
              a(ik,kk)=a(ik,kk)+w(k)*phi(i+6)*phi(j)*yc*y
              a(ik,nn)=a(ik,nn)-w(k)*phi(i+6)*phi(j)*xc*y      ! JK_Full_Expression(same as before)
	    END DO                                                   
             rhs(ik)=rhs(ik)+(Y*(-(dyt*Xc) + V*Xc + (dxt 
     1              - U)*Yc)*phi(i+6)*w(k))         !JK_Full_Expression
	  END DO
          yfs=phi(7)*aypt(7)+phi(8)*aypt(8)+phi(9)*aypt(9)
          vol=vol+w(k)*xc*yfs
	END DO ! fin <K>
      	IF(NELL .EQ. (NEX*NEY))then
          ekin=ekin/2.
c          write(36,*)vol,ekin,esup
        end if
        
	END IF ! fin <line integrals>
	
        IF(Nell.eq.ney)then
	  CALL TFUNCT(0.0,1.0)
	  CALL DERIV
	  CALL DERJAC
c            rhs(V2(7))=rhs(V2(7))+Y*Yc/(dsqrt(Xc**2+Yc**2))
c            rhs(V1(7))=rhs(V1(7))+Y*Xc/(dsqrt(Xc**2+Yc**2))
	  DO J=1,3
	    JK=22+J
              Xh=dh(j)*dirl(j,1)
              Yh=dh(j)*dirl(j,2)  !dY/dh Y=y(j)dh(j)dirl(j,2) j=1,3
              Xeh=dhe(j)*dirl(j,1)
              Xch=dhc(j)*dirl(j,1)
              Yeh=dhe(j)*dirl(j,2)
              Ych=dhc(j)*dirl(j,2)
c              a(V2(7),JK)=a(V2(7),JK)-((Yc/(Xc**2 + Yc**2)**0.5)*Yh
c     *                    +(Y*((-1.*Yc**2)/(Xc**2 + Yc**2)**1.5 
c     *                    + (Xc**2 + Yc**2)**(-0.5)))*Ych
c     *                    +((-1.*Xc*Y*Yc)/(Xc**2 + Yc**2)**1.5)*Xch)
c              a(V1(7),JK)=a(V1(7),JK)-((Xc/(Xc**2 + Yc**2)**0.5)*Yh
c     *                    + ((-1.*Xc*Y*Yc)/(Xc**2 + Yc**2)**1.5)*Ych
c     *                    + (Y*((-1.*Xc**2)/(Xc**2 + Yc**2)**1.5 
c     *                    + (Xc**2 + Yc**2)**(-0.5)))*Xch)
          END DO    

c            rhs(V2(7))=rhs(V2(7))-Y
c            rhs(V1(7))=rhs(V1(7))+Y*Xc/(dsqrt(Xc**2+Yc**2))
	  DO J=1,3
	    JK=22+J
              Xh=dh(j)*dirl(j,1)
              Yh=dh(j)*dirl(j,2)  !dY/dh Y=y(j)dh(j)dirl(j,2) j=1,3
              Xeh=dhe(j)*dirl(j,1)
              Xch=dhc(j)*dirl(j,1)
              Yeh=dhe(j)*dirl(j,2)
              Ych=dhc(j)*dirl(j,2)
c              a(V2(7),JK)=a(V2(7),JK)-(-Yh)
c     *                    +(Y*((-1.*Yc**2)/(Xc**2 + Yc**2)**1.5 
c     *                    + (Xc**2 + Yc**2)**(-0.5)))*Ych
c     *                    +((-1.*Xc*Y*Yc)/(Xc**2 + Yc**2)**1.5)*Xch)
c              a(V1(7),JK)=a(V1(7),JK)-((Xc/(Xc**2 + Yc**2)**0.5)*Yh
c     *                    + ((-1.*Xc*Y*Yc)/(Xc**2 + Yc**2)**1.5)*Ych
c     *                    + (Y*((-1.*Xc**2)/(Xc**2 + Yc**2)**1.5 
c     *                    + (Xc**2 + Yc**2)**(-0.5)))*Xch)
          END DO    
        End IF
CcccccccccccccccFREE BOUNDARY CONDITION 040913JK-Add free boundary condition
	IF(NELL .le. 0)THEN 

         IF(Nell.eq.ney)then
           rstemp=rhs(V2(7))
           Do i=1,25
             astemp1(i)=a(V2(7),i)
             astemp2(i)=a(V1(7),i)
           End Do
           rstemp2=rhs(V1(7))
	 End IF
	DO K=1,3
	  CALL TFUNCT(0.0,GP(K))
	  CALL DERIV
	  CALL DERJAC
	  DO I=1,7,3
            m=v1(i)
            n=v2(i)

          strainrate=        (Ve**2*Y**2*(2*Xc**2 + Yc**2) 
     *      + Ue**2*Y**2*(Xc**2 + 2*Yc**2) + 
     *    Xe**2*((Uc**2 + 2*Vc**2)*Y**2 + 2*V**2*Yc**2) - 
     *    2*Xe*(Uc*Vc*Y**2 + 2*V**2*Xc*Yc)*Ye + 
     *    (2*V**2*Xc**2 + (2*Uc**2 + Vc**2)*Y**2)*Ye**2 - 
     *    2*Ue*Y**2*(Uc*Xc*Xe + Ve*Xc*Yc - Vc*Xc*Ye + 2*Uc*Yc*Ye) - 
     *    2*Ve*Y**2*(2*Vc*Xc*Xe - Uc*Xe*Yc + Vc*Yc*Ye))/
     *  (Y**2*(Xe*Yc - Xc*Ye)**2)

          eta= beta+(1-beta)*((1+alp*alp*strainrate)**((pn-1.d0)/2.d0))
  
          etap= ((pn-1.d0)/2.d0)*(1-beta)*(alp*alp)*((1
     *          +alp*alp*strainrate)**((pn-3.d0)/2))!eta prime with respcet to strainrate
 
            rhs(m)=rhs(m)+w(k)*(eta*phi(i)*(Xc*Ue-Xe*Uc-Yc*Ve)*Y/Yc)
	    rhs(n)=rhs(n)+w(k)*(2*eta*phi(i)*(Xc*Ve-Xe*Vc)*Y/Yc
     *          +phi(i)*P*Y*Xe/We)

          DO J=1,9
	    KK=V1(J)
	    NN=V2(J)
c eta primes with respect to Uj, Vj
            etaup=etap*(0.0*phi(j)
     *            +((2*Xe*(-(Ue*Xc) + Uc*Xe + Ve*Yc) 
     *            - 2*(Vc*Xe + 2*Ue*Yc)*Ye + 4*Uc*Ye**2)/
     *            (Xe*Yc - Xc*Ye)**2)*phic(j)
     *            +((2*(-(Ve*Xc*Yc) + Ue*(Xc**2 + 2*Yc**2) 
     *            + Vc*Xc*Ye - Uc*(Xc*Xe + 2*Yc*Ye)))/
     *            (Xe*Yc - Xc*Ye)**2)*phie(j)) !d/dU*dU/dUj+d/dUc*dUc/dUj+dUe/dUj 

            etavp=etap*(((4*V)/Y**2)*phi(j)+((2*(-(Ue*Xc*Yc) 
     *           + Uc*Xe*Yc + Ve*(2*Xc**2 + Yc**2) - Vc*(2*Xc*Xe 
     *           + Yc*Ye)))/(Xe*Yc - Xc*Ye)**2)*phie(j)
     *           +((2*(2*Vc*Xe**2 + Ue*Xc*Ye - Uc*Xe*Ye 
     *           + Vc*Ye**2 - Ve*(2*Xc*Xe + Yc*Ye)))/(Xe*Yc 
     *           - Xc*Ye)**2)*phic(j)) !d/dV*dV/dVj+d/dVc*dVc/dVj+dVe/dVj

             a(m,kk)=a(m,kk)-w(k)*phi(i)*(etaup*(phi(i)*(Xc*Ue
     *               -Xe*Uc-Yc*Ve)*Y/Yc)+eta*(((Xc*Y)/Yc)*phie(j)
     *               +(-((Xe*Y)/Yc))*phic(j)))
             a(m,nn)=a(m,nn)-w(k)*phi(i)*(etavp*(phi(i)*(Xc*Ue
     *               -Xe*Uc-Yc*Ve)*Y/Yc)+eta*((-Y)*phie(j)
     *               +(0.0*phic(j))))
             a(n,kk)=a(n,kk)-w(k)*phi(i)*2*etaup*((Xc*Ve-Xe*Vc)*Y/Yc)
             a(n,nn)=a(n,nn)-w(k)*phi(i)*2*(etavp*((Xc*Ve-Xe*Vc)*Y/Yc)
     *               +eta*(((Xc*Y)/Yc)*phie(j)+(-((Xe*Y)/Yc))*phic(j)))

          END DO 
	  DO J=1,4
	    KK=V3(J)
            a(n,kk)=a(n,kk)-w(k)*psi(j)*phi(i)*Y*Xe/We
          END DO
	    DO J=1,3
	      JK=J+22
              Xh=dh(j)*dirl(j,1)
              Yh=dh(j)*dirl(j,2)  !dY/dh Y=y(j)dh(j)dirl(j,2) j=1,3
              Xeh=dhe(j)*dirl(j,1)
              Xch=dhc(j)*dirl(j,1)
              Yeh=dhe(j)*dirl(j,2)
              Ych=dhc(j)*dirl(j,2)

              etady=(-4*V**2)/Y**3
              etadye= (2*(Ve**2*Xc*(2*Xc**2 + Yc**2) 
     -                + Ue**2*(Xc**3 + 2*Xc*Yc**2) + 
     -                Xe*(-(Uc*Vc*(Xe*Yc + Xc*Ye)) 
     -                + Vc**2*(2*Xc*Xe + Yc*Ye) + 
     -                Uc**2*(Xc*Xe + 2*Yc*Ye)) - 
     -                Ve*(-2*Uc*Xc*Xe*Yc + Vc*(4*Xc**2*Xe 
     -                + Xe*Yc**2 + Xc*Yc*Ye)) + 
     -                Ue*(Xc*(-2*Ve*Xc*Yc + Vc*Xe*Yc 
     -                + Vc*Xc*Ye) -2*Uc*(Xe*(Xc**2 + Yc**2) 
     -                + Xc*Yc*Ye))))/(Xe*Yc - Xc*Ye)**3
              etadyc=(-2*(Ve**2*Xc*(2*Xc*Xe + Yc*Ye) 
     -              + Ue**2*Xc*(Xc*Xe + 2*Yc*Ye) + 
     -              Xe*(-2*Uc*Vc*Xe*Ye + Vc**2*(2*Xe**2 + Ye**2) + 
     -              Uc**2*(Xe**2 + 2*Ye**2)) - 
     -              Ue*(Xc*(Ve*Xe*Yc + Ve*Xc*Ye - 2*Vc*Xe*Ye) + 
     -              2*Uc*(Xe*Yc*Ye + Xc*(Xe**2 + Ye**2))) + 
     -              Ve*(Uc*Xe*(Xe*Yc + Xc*Ye) - Vc*(Xe*Yc*Ye 
     -              + Xc*(4*Xe**2 + Ye**2)))))/(Xe*Yc - Xc*Ye)**3
              etadx=0.0
              etadxe= (-2*(Ue**2*Yc*(Xc**2 + 2*Yc**2) 
     -              + Ve**2*(2*Xc**2*Yc + Yc**3) + 
     -              Ye*(-(Uc*Vc*(Xe*Yc + Xc*Ye)) 
     -              + Vc**2*(2*Xc*Xe + Yc*Ye) + 
     -              Uc**2*(Xc*Xe + 2*Yc*Ye)) - 
     -              Ue*(2*Xc*Yc*(Ve*Yc - Vc*Ye) + 
     -              Uc*(Xc*Xe*Yc + Xc**2*Ye + 4*Yc**2*Ye)) + 
     -              Ve*(Uc*Yc*(Xe*Yc + Xc*Ye) - 2*Vc*(Xc*Xe*Yc 
     -              + (Xc**2 + Yc**2)*Ye))))/(Xe*Yc - Xc*Ye)**3
              etadxc=(2*(2*Ve*Xe*(Ve*Xc - Vc*Xe)*Yc + 
     -              (-2*Vc*Ve*Xc*Xe + 2*Vc**2*Xe**2 
     -              + (Uc*Xe + Ve*Yc)**2)*Ye - 
     -              2*Vc*(Uc*Xe + Ve*Yc)*Ye**2 
     -              + (2*Uc**2 + Vc**2)*Ye**3 + 
     -              Ue**2*Yc*(Xc*Xe + 2*Yc*Ye) - 
     -              Ue*((Ve*Yc - Vc*Ye)*(Xe*Yc + Xc*Ye) + 
     -              Uc*(Xe**2*Yc + Xc*Xe*Ye 
     -              + 4*Yc*Ye**2))))/(Xe*Yc - Xc*Ye)**3
              etah=etap*(etady*Yh+etadye*Yeh+etadyc*Ych
     *                   etadx*Xh+etadxe*Xeh+etadxc*Ych)
              a(m,jk)=a(m,jk)-w(k)*phi(i)*((etah*(Xc*Ue
     *               -Xe*Uc-Yc*Ve)*Y/Yc)
     *               +eta*((-((Uc*Y)/Yc))*Xeh
     *               +((Ue*Y)/Yc)*Xch
     *               +((Ue*Xc - Uc*Xe - Ve*Yc)/Yc)*Yh
     *               +((-(Ue*Xc*Y) + Uc*Xe*Y)/Yc**2)*Ych))

              a(n,jk)=a(n,jk)-w(k)*phi(i)*(2*etah*((Xc*Ve
     *               -Xe*Vc)*Y/Yc)+2*eta*((-((Vc*Y)/Yc))*Xeh
     *               +((Ve*Y)/Yc)*Xch
     *               +((Ve*Xc - Vc*Xe)/Yc)*Yh
     *               +((-(Ve*Xc*Y) + Vc*Xe*Y)/Yc**2)*Ych)
     *               +((P*Xe)/We)*Yh
     *               +((P*Y)/We)*Xeh)

	    END DO
	  END DO

	END DO ! fin <K>
        
         IF(Nell.eq.ney)then
           rhs(V2(7))=rstemp
           Do i=1,25
            a(V2(7),i)= astemp1(i)
            a(V1(7),i)= astemp2(i)
           End Do
          rhs(V1(7))= rstemp2
	 End IF
	END IF ! fin FBC
C     ensamble residuos globales	
	DO I=1,9
	  R2(W1(I))=R2(W1(I))+RHS(V1(I))
	  R2(W2(I))=R2(W2(I))+RHS(V2(I))
	END DO
	DO I=1,4
	  R2(W3(I))=R2(W3(I))+RHS(V3(I))
	END DO
	IF(IELFS(NELL) .EQ. 1)THEN
	  DO I=1,3
	    R2(W4(I))=R2(W4(I))+RHS(22+I)
	  END DO
	END IF
        RETURN
        END

       SUBROUTINE CALVELOT
       !This subroutine calculates the interfacial velocities, i.e. tangential and normal velocity to the surface. And the curvatures
       !Author: JK
       !Last modified: 01-17-2013
        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        INTEGER W1,W2,W3,W4
	DIMENSION JMOD(500),DSOLV(15999),NCOD(15999)
	COMMON/FRON/NP,NH,NHADD,NE,NBN(1700),NCN(1700)
        DIMENSION W(3),GP(3)
	COMMON/TFUN1/PHI(9),PSI(4),PHIC(9),PHIE(9),PHICC(9)
        COMMON/DER1/AXPT(9),AYPT(9),WYL(3),DIRL(3,2),WXL(3)
        COMMON/DER2/U,UC,UE,V,VC,VE,P
        COMMON/DER3/XE,XC,YC,YE,DET,Y
	COMMON/FRON2/NOP(1700,12),NOPP(7182),MDF(7182)
        COMMON/GEN1/SOL(15999),W1(9),W2(9),W5(9),W3(4),W4(3)
	COMMON/MESCO/YPT(7011),XPT(7011)
        COMMON/MESH1/NNY,NNX,NNX1,NNX2                  !JK_new_mesh
        COMMON/MESHCON/NEY,NEX,NEX1,NEX2                !JK_new_mesh
	COMMON/FSARRAY/H(180),BCP(180,2)
	COMMON/TEMP2/SOLV(15999),DSOL(15999)
	COMMON/TEMP3/CT1,CT2,STEPT,TIME,TIMEF,STEPTV,NTIME  !moidified for continuation run: JK
        COMMON/ZONES/NEXZ1,NEYZ1,NEZ1,NNXZ1,NNYZ1
       DO NEELEM=NEYZ1,NEXZ1*NEYZ1,NEYZ1
          DO I=1,9
             J=NOP(NEELEM,I)
             AXPT(I)=XPT(J)
             AYPT(I)=YPT(J)
          END DO
          IF (NEELEM.EQ.NEYZ1) THEN
             write(39,371) time
             cero=0.0d0
             CALL TFUNCT(cero,1.0d0)
             ZCC=0.
             RCC=0.
             ZC=0.
             RC=0.
             Rsup=0.! modificado
             VeloU=0.
             VeloV=0.
             VeloUC=0.
             VeloVC=0.
             DO J=7,9
                ZCC=ZCC+AXPT(J)*PHICC(J)
                RCC=RCC+AYPT(J)*PHICC(J)
                ZC=ZC+AXPT(J)*PHIC(J)
                RC=RC+AYPT(J)*PHIC(J)
                Rs =Rs +AYPT(J)*PHI(J) ! modificado
                VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
                VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
                VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
                VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
              END DO
              RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
              RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
              velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
              velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5
              write(39,310) AXPT(7),AYPT(7),RKappa1,RKappa2,RKappa1
     1 +RKappa2,velonorm,velotan

         END IF

         CALL TFUNCT(0.5d0,1.d0)
         ZCC=0.
         RCC=0.
         ZC=0.
         RC=0.
         Rs=0.! modificado
         VeloUC=0.
         VeloVC=0.
         VeloU=0.
         VeloV=0.
!         ROSC=0.
!         ROS=0.
         DO J=7,9
            ZCC=ZCC+AXPT(J)*PHICC(J)
            RCC=RCC+AYPT(J)*PHICC(J)
            ZC=ZC+AXPT(J)*PHIC(J)
            RC=RC+AYPT(J)*PHIC(J)
            Rs =Rs +AYPT(J)*PHI(J) ! modificado
            VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
            VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
            VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
            VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
         END DO
         RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
         RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
         velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
         velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5

         write(39,310) AXPT(8),AYPT(8),RKappa1,RKappa2,RKappa1
     1 +RKappa2,velonorm,velotan

         CALL TFUNCT(1.d0,1.d0)
         ZCC=0.
         RCC=0.
         ZC=0.
         RC=0.
         Rs=0.! modificado
         VeloUC=0.
         VeloVC=0.
         VeloU=0.
         VeloV=0.
         ROSC=0.
         ROS=0.
         DO J=7,9
            ZCC=ZCC+AXPT(J)*PHICC(J)
            RCC=RCC+AYPT(J)*PHICC(J)
            ZC=ZC+AXPT(J)*PHIC(J)
            RC=RC+AYPT(J)*PHIC(J)
            Rs =Rs +AYPT(J)*PHI(J) ! modificado
            VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
            VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
            VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
            VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
         END DO
         RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
         RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
         velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
         velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5

         write(39,310) AXPT(9),AYPT(9),RKappa1,RKappa2,RKappa1
     1 +RKappa2 ,velonorm,velotan
        END Do

       DO NEELEM=NEY+NEXZ1*NEYZ1,NEX*NEY+NEXZ1*NEYZ1,NEY
          DO I=1,9
             J=NOP(NEELEM,I)
             AXPT(I)=XPT(J)
             AYPT(I)=YPT(J)
          END DO
c          IF (NEELEM.EQ.NEY) THEN
c             write(39,371) time
c             cero=0.0d0
c             CALL TFUNCT(cero,1.0d0)
c             ZCC=0.
c             RCC=0.
c             ZC=0.
c             RC=0.
c             Rsup=0.! modificado
c             VeloU=0.
c             VeloV=0.
c             VeloUC=0.
c             VeloVC=0.
c             DO J=7,9
c                ZCC=ZCC+AXPT(J)*PHICC(J)
c                RCC=RCC+AYPT(J)*PHICC(J)
c                ZC=ZC+AXPT(J)*PHIC(J)
c                RC=RC+AYPT(J)*PHIC(J)
c                Rs =Rs +AYPT(J)*PHI(J) ! modificado
c                VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
c                VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
c                VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
c                VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
c              END DO
c              RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
c              RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
c              velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
c              velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5
c              write(39,310) AXPT(7),AYPT(7),RKappa1,RKappa2,RKappa1
c     1 +RKappa2,velonorm,velotan
c
c         END IF

         CALL TFUNCT(0.5d0,1.d0)
         ZCC=0.
         RCC=0.
         ZC=0.
         RC=0.
         Rs=0.! modificado
         VeloUC=0.
         VeloVC=0.
         VeloU=0.
         VeloV=0.
!         ROSC=0.
!         ROS=0.
         DO J=7,9
            ZCC=ZCC+AXPT(J)*PHICC(J)
            RCC=RCC+AYPT(J)*PHICC(J)
            ZC=ZC+AXPT(J)*PHIC(J)
            RC=RC+AYPT(J)*PHIC(J)
            Rs =Rs +AYPT(J)*PHI(J) ! modificado
            VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
            VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
            VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
            VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
         END DO
         RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
         RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
         velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
         velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5

         write(39,310) AXPT(8),AYPT(8),RKappa1,RKappa2,RKappa1
     1 +RKappa2,velonorm,velotan

         CALL TFUNCT(1.d0,1.d0)
         ZCC=0.
         RCC=0.
         ZC=0.
         RC=0.
         Rs=0.! modificado
         VeloUC=0.
         VeloVC=0.
         VeloU=0.
         VeloV=0.
         ROSC=0.
         ROS=0.
         DO J=7,9
            ZCC=ZCC+AXPT(J)*PHICC(J)
            RCC=RCC+AYPT(J)*PHICC(J)
            ZC=ZC+AXPT(J)*PHIC(J)
            RC=RC+AYPT(J)*PHIC(J)
            Rs =Rs +AYPT(J)*PHI(J) ! modificado
            VeloU=VeloU + SOL(NOPP(NOP(NEELEM,J)))*PHI(J)
            VeloV=VeloV + SOL(NOPP(NOP(NEELEM,J))+1)*PHI(J)
            VeloUC=VeloUC + SOL(NOPP(NOP(NEELEM,J)))*PHIC(J)
            VeloVC=VeloVC + SOL(NOPP(NOP(NEELEM,J))+1)*PHIC(J)
         END DO
         RKappa1=(-zcc*rc+zc*rcc)/(RC**2+ZC**2)**1.5
         RKappa2=-zc/(Rs*(RC**2+ZC**2)**0.5)
         velonorm=(ZC*VeloV-RC*VeloU)/(RC**2+ZC**2)**0.5
         velotan=(RC*VeloV+ZC*VeloU)/(RC**2+ZC**2)**0.5

         write(39,310) AXPT(9),AYPT(9),RKappa1,RKappa2,RKappa1
     1 +RKappa2 ,velonorm,velotan
        END Do




309    format('#',9X,'X:',15X,'Y:',15X,'Rkappa1',16X,'RKappa2',15X
     1,'RKappa',15X,'Vn',15X,'Vt')
310    format(e20.10,e20.10,e20.10,e20.10,e20.10,e20.10,e20.10)
311    format(e15.5,e20.10,e20.10,e20.10,e20.10,e20.10)
312    format('#',9X,'X:',15X,'Y:',15X,'time',16X,'Rkappa2',15X
     1,'RKappa1',15X,'RKappa')
371    FORMAT('ZONE T="t=',E15.5,'", F=POINT')	


	RETURN
	END

