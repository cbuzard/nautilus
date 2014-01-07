! -----------------------------------------------------------------------
!   #     #     #     #     #  #######  ###  #        #     #   #####   
!   ##    #    # #    #     #     #      #   #        #     #  #     #  
!   # #   #   #   #   #     #     #      #   #        #     #  #        
!   #  #  #  #     #  #     #     #      #   #        #     #   #####   
!   #   # #  #######  #     #     #      #   #        #     #        #  
!   #    ##  #     #  #     #     #      #   #        #     #  #     #  
!   #     #  #     #   #####      #     ###  #######   #####    #####   
! -----------------------------------------------------------------------
!                            ?
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~|^"~~~~~~~~~~~~~~~~~~~~~~~~~o~~~~~~~~~~~
!        o                   |                  o      __o
!         o                  |                 o     |X__>
!       ___o                 |                __o
!     (X___>--             __|__            |X__>     o
!                         |     \                   __o
!                         |      \                |X__>
!  _______________________|_______\________________
! <                                                \____________   _
!  \                                                            \ (_)
!   \    O       O       O                                       >=)
!    \__________________________________________________________/ (_)
!
!                            ___
!                           / o \
!                      __   \   /   _
!                        \__/ | \__/ \
!                       \___//|\\___/\
!                        ___/ | \___
!                             |     \
!                            /
! -----------------------------------------------------------------------
! 
! A fast 1D gas-grain chemical model by FH (2008)
! Based upon the OSU gas-grain chemical model
! Uses the OSU chemical network
! Updated from gg_osu_2006v1d by RTG/VW
! Rate equations from Hasegawa & Herbst (1992)
! Modified rates following Caselli et al. (1998)
! Stiff solver for sparse Jacobians: LSODES/ODEPACK (Hindmarsh 1983)
! Turbulent mixing implemented through first order operator splitting
! 
! Files:
! nautilus.f90          The main program
! nls_control.d         Control parameters
! nls_header_mod.f90    Declarations and constants
! nls_io.f90            Input/output
! nls_odes.f90          Rate equations and jacobian
! nls_phys_1D.f90       1D structure (overwrites some parameters)
! nls_mesh.f90          1D mesh
! nls_diffusion.f90     Turbulent diffusion routines
!
! opkd*.f               Odepack files
! 
! nls_surf_input.d      Surface parameters (desorption energies, etc..)
! gg_react_osu_03_2008  Chemical network
! gg_photodiss*         H2 and CO photodissociation
! nls_init.d            Initial abundances (must be present even if IREAD=0)
!
! nlso_srates.d         Output of some surface quantities (incomplete)
! nlso_tail.d           Chemical compasition at the last timestep (if used in 0D)
! output_1D.******      1D outputs (binaries) 
! rates1D.******        Outputs of the rate coefficients for 1 specific mesh point
! nlso_mod.d            Output of the modified rate coefficients (not checked)
! nlso_spec.d           List of species with their names and numbers
! 
! April 2011 VW  An appoximative calculation of the X-ray ionization rate 
! have been added according to Glasgold at al. (1999)
! 
! -----------------------------------------------------------------------

      PROGRAM Gasgrain

      use header
      use unitvar
      IMPLICIT none

      integer :: lrw, liw
 
      real(kind=8), dimension(NSMAX) :: Y

      integer :: itol, itask, istate, iopt, mf
      real(kind=8) :: atol
      real(kind=8) :: T, TOUT, TIN
      real(kind=8) :: tdeb,tfin

      data itol, itask, istate, iopt, mf, atol/2,1,1,1,021,1.d-99/

      call CPU_TIME(tdeb) ! To get computation time if needed
      CALL FILESET
      CALL READINPUT

! Dimension of the work arrays for the solver 
! The number of non zero values is checked with the testjac flag
! NJAC should be around the largest printed value

      if (testjac.eq.1) then
      print *,'--------------------------------------------'
      print *,'Dummy run to check'
      print *,'the number of non zero elements per column !'
      print *,'--------------------------------------------'
      lrw = 20 + 9*NSMAX*NPTMAX + 10000*NSMAX*NPTMAX
      liw = 31 + NSMAX*NPTMAX + 10000*NSMAX*NPTMAX
      else
      lrw = 20 + 3 * NJAC*NSMAX + 21 * NSMAX
      liw = 31 + 3 * NJAC*NSMAX + 21 * NSMAX
      endif
      
! Build spatial mesh 
      call mesh

! Initialization of elemental/chemical quantities

      CALL INITIAL

      IT=0
      CALL START(TOUT)

! 1D physical structure (nls_phys_1D)
      call phys_1D

! Write species name/index correspondance
      CALL WRITESPEC

! Global variable for the current spatial point
      iptstore=1

! Initialize indices of reactants and products 

      call chemsetup

! Initializing ZXN
      do ipts=1,nptmax
      ZXN(:,ipts) = XN(:)
      enddo

! Initializing T
! T = local, TIME = global
      T = 0.d0
      TIME = 0.d0

! Allocate the JA array
      allocate(JA(1:liw))

! The real time loop
      do while (t.lt.0.9*tfinal)

      call ztimestep(zdt) ! Determination of the diffusive timestep 
      TIME=T+ZDT ! Final time for chemistry T -> T + ZDT

      it=it+1

! Store the current time in TIN (for 1D calculation)
      TIN = T

         WRITE(NERR,5) 'IT=',IT,', TIME=',TIME/TYEAR,' yrs'
    5    FORMAT (A3,I5,A7,1PD10.3,A4)

      do ipts=1,nptmax ! Start of the spatial loop for chemistry

      iptstore = ipts

! T being changed in dlsode, needs to be defined again

      T=TIN

! Feed 1D physical structure
      TAU=TAU1D(ipts)
      ZETAX=ZETAX1D(ipts)
! XNT is the total density of H atoms (both H and H2)
! But DENS1D is the real number density
! Problem with the sound speed if NH > NH2 (cf phys_1D)
      XNT=2.*DENS1D(ipts)
      TEMP=TEMP1D(ipts)
      DTEMP=DTEMP1D(ipts)

! Chemical evolution for each spatial point

      Y(:nsmax) = ZXN(:,ipts)
      XN(:nsmax) = ZXN(:,ipts)

      call evolve (T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

! Output of the rates once every 10 chemical outputs
!      if ((mod(it,wstepr).eq.0).and.(ipts.eq.irateout)) then
      call rates1D
!      endif

      if (istate.eq.-3) stop

      ZXN(:,ipts) = XN(:) ! putting back

      enddo ! end of the spatial loop for chemistry 

! Generic call to zdiffusion, a subroutine calling the chosen numerical scheme 
      call zdiffusion ! diffusion 

      if (mod(it,wstep).eq.0) then
      call write1D
      endif

      ENDDO
 
      if (nptmax.eq.1) call writetail

      call CPU_TIME(tfin)
      print*,'CPU time=',(tfin-tdeb)/60.0 ,'min'
      END

! ======================================================================
! ======================================================================
      SUBROUTINE FILESET
      USE unitvar
      use header
      implicit none


      OPEN (UNIT=NCON,FILE='nls_control.d',STATUS='OLD')
      OPEN (UNIT=NMOD,FILE='nlso_mod.d',STATUS='UNKNOWN')
      OPEN (UNIT=NSP, FILE='nlso_spec.d',STATUS='UNKNOWN')
      OPEN (UNIT=NGR, FILE='nls_surf_fev2012.dat',STATUS='OLD')
      OPEN (UNIT=NJR, FILE='nls_gas_fev2012.dat',STATUS='OLD')
!      OPEN (UNIT=NJR, FILE='nls_gas_update.dat',STATUS='OLD')
      OPEN (UNIT=NJR2, FILE='nls_grain_fev2012.dat',STATUS='OLD')
      OPEN (UNIT=NTAI,FILE='nlso_tail.d',STATUS='UNKNOWN')
      OPEN (UNIT=NINI,FILE='nls_init.d',STATUS='OLD') 
      OPEN (UNIT=CODIS,FILE='gg_CO_Photodiss.d',STATUS='OLD')
      OPEN (UNIT=H2DIS,FILE='gg_H2_Photodiss.d',STATUS='OLD')

      RETURN 
      END

! ======================================================================
! ======================================================================
      SUBROUTINE INITIAL
      use header
      use unitvar
      IMPLICIT none

      REAL(kind=8), dimension(NEMAX) :: MASS
      real(kind=8) :: MSUM
      integer :: ILAB, KSUM, j, k, i, isptemp

! Set elements' characteristics=========================================
! --- Find the atomic species associated with a given element

      ILAB=1
      DO J=1,NSMAX
         KSUM=0
! ------ Calculate species' mass
         DO K=1,NEMAX
            KSUM=KSUM+IELM(K,J)
         ENDDO
! ------ Check for atomic species
         IF ((KSUM.EQ.1).AND.(ICG(J).EQ.0).AND.&
         (SPEC(J)(:1).NE.'J          ').AND.(SPEC(J)(:1).NE.'X          ')) THEN
            IF (ILAB.GT.NEMAX) then
                  STOP '***More elements than NEMAX***'
           endif       
! --------- Save species number
            ISPELM(ILAB)=J
            ILAB=ILAB+1
         ENDIF
       
! ------ Check for electron species number
         IF (SPEC(J).EQ.'e-         ') ISPE=J
      ENDDO

! --- Re-arrange order of elements to match IELM columns (reactions file)
      DO J=1,NEMAX-1
         IF (IELM(J,ISPELM(J)).NE.1) THEN
            DO K=J+1,NEMAX
               IF (IELM(J,ISPELM(K)).EQ.1) THEN
                  ISPTEMP=ISPELM(K)
                  ISPELM(K)=ISPELM(J)
                  ISPELM(J)=ISPTEMP
               ENDIF
            ENDDO
         ENDIF
      ENDDO

! --- Set elements' masses
      DO I=1,NEMAX
         IF (SPEC(ISPELM(I)).EQ.'H          ') MASS(I)=1
         IF (SPEC(ISPELM(I)).EQ.'D          ') MASS(I)=2
         IF (SPEC(ISPELM(I)).EQ.'He         ') MASS(I)=4
         IF (SPEC(ISPELM(I)).EQ.'C          ') MASS(I)=12
         IF (SPEC(ISPELM(I)).EQ.'N          ') MASS(I)=14
         IF (SPEC(ISPELM(I)).EQ.'O          ') MASS(I)=16
         IF (SPEC(ISPELM(I)).EQ.'Na         ') MASS(I)=23
         IF (SPEC(ISPELM(I)).EQ.'Mg         ') MASS(I)=24
         IF (SPEC(ISPELM(I)).EQ.'Si         ') MASS(I)=28
         IF (SPEC(ISPELM(I)).EQ.'P          ') MASS(I)=31
         IF (SPEC(ISPELM(I)).EQ.'S          ') MASS(I)=32
         IF (SPEC(ISPELM(I)).EQ.'Cl         ') MASS(I)=35
         IF (SPEC(ISPELM(I)).EQ.'Fe         ') MASS(I)=56
         IF (SPEC(ISPELM(I)).EQ.'F          ') MASS(I)=19
      ENDDO

! Set species' characteristics==========================================
! --- Set special species labels
      YH     = 'H          '
      YJH    = 'JH         '
      YH2    = 'H2         '
      YJH2   = 'JH2        '
      YHE    = 'He         '
      YHEP   = 'He+        '
      YE     = 'e-         '
      YGRAIN = 'GRAIN0     '
      YCO    = 'CO         '

! --- Set reference species
      DO I=1,NSMAX 
! ------ Calculate masses
         MSUM=0
         DO K=1,NEMAX 
            MSUM=MSUM+MASS(K)*IELM(K,I) 
         ENDDO 
         AWT(I)=MSUM
         IF (SPEC(I).EQ.YE) AWT(I)=1.D+0/1836.D+0 
         IF (SPEC(I).EQ.YGRAIN .OR. SPEC(I).EQ.'GRAIN-      ')&
                    AWT(I)=4.0*PI*RD*RD*RD*RHOD/3.0/AMU
      ENDDO

! Initialize the Av/NH ratio
! Can be scaled for different dust/gas ratios
! Here DTOGM is the original dust/gas mass ratio (from nls_control.d)
! DTOGM is changed later into the dust/hydrogen mass ratio

      RAVNH=5.34d-22*(DTOGM/1.d-2)

! Find ITYPE first and last reactions===================================
      DO I=0,NITYPE-1
         IRXSTA(I)=0
         IRXFIN(I)=0
         DO J=1,NKMAX
            IF ((ITYPE(J).EQ.I).AND.(IRXSTA(I).EQ.0)) IRXSTA(I)=J
            IF (ITYPE(J).EQ.I) IRXFIN(I)=J
         ENDDO
      ENDDO

! Find the index of CO and H2
      do i=1,nsmax
      if (SPEC(i).eq.YH2) INDH2=i
      if (SPEC(i).eq.YCO) INDCO=i
      if (SPEC(i).eq.YHE) INDHE=i
      enddo

! Compute TNS = number of sites per grain
      TNS = SNS*4.d0*PI*RD**2

! Initialise reaction rates=============================================
      CALL GRAINRATE

      RETURN 
      END

! ======================================================================
! ======================================================================
      SUBROUTINE START(TOUT)

      use header
      use unitvar
      IMPLICIT none

      CHARACTER(len=11), dimension(NSMAX) :: SREAD
      REAL(kind=8), dimension(NSMAX) :: Y

      integer :: i, j, k
      real(kind=8) :: TOUT

! Initialise times======================================================
      TOUT=0.0D+0
      TIME=0.0D+0

! Set initial abundances================================================
      DO I=1,NSMAX
         XN(I)=XNMIN
         DO K=1,NS0
            IF (SPEC(I).EQ.XS0(K)) THEN
               XN(I)=XN0(K)
            ENDIF
         ENDDO
      ENDDO



! Compute elemental abundances

      DO J=1,NEMAX
         ELEMS(J)=0.0D+0
      ENDDO
      DO I=1,NSMAX 
         DO J=1,NEMAX
            ELEMS(J)=ELEMS(J)+IELM(J,I)*XN(I)
         ENDDO
      ENDDO

! Recompute DTOGM to remove He
! In the following, DTOGM is used as a H/dust mass ratio
      do i=1,nemax
      IF (SPEC(ISPELM(I)).EQ.YHE) then
      DTOGM = DTOGM*(1.d0+4*ELEMS(I))
! Mean molecular weight (cgs) 
! Approximated here (the exact calculus would require a sume over AWT
! Used for the diffusion in disks, not for chemistry
      MEANW = 2.d0 + 4.d0*ELEMS(I)/(1.d0+ELEMS(I))
      endif
      enddo

! Compute the grain abundance

      GTODN=(4.D+0*PI*RHOD*RD*RD*RD)/(3.D+0*DTOGM*AMU)

      where(SPEC.EQ.YGRAIN) XN=1.0/GTODN


! Read initial abundances if IREAD is switched on=======================
      IF (IREAD.NE.0) THEN

! Skip the first line on nls_init.d
         READ (NINI,*) 

! ------ Read species and abundances
         READ (NINI,110) (SREAD(I),XN(I),I=1,NSMAX)
  110 FORMAT (5(A11,2X,1PE12.5,2X))
         READ (NINI,*)
! ------ Check if species in nls_init.d correspond to the reaction file
         DO I=1,NSMAX
            IF (SREAD(I).NE.SPEC(I)) THEN
               WRITE (NERR,*) 'Input species in init file ',&
                              'do not match those in reaction file'
               STOP
            ENDIF
         ENDDO
      ENDIF


! Set the electron abundance via conservation===========
! And check at the same time that nls_init has the same elemental abundance
! as nls_control

      Y(:)=XN(:)
      CALL CONSERVE(Y)
      XN(:)=Y(:) 

      RETURN 
      END 

! ======================================================================
! ======================================================================

      SUBROUTINE EVOLVE (T,Y,TOUT,itol,atol,itask,istate,iopt,mf,liw,lrw)

      use header
      use unitvar
      IMPLICIT none

      integer :: liw,lrw

      INTEGER, dimension(liw) :: IWORK
      REAL(kind=8), dimension(NSMAX) :: Y
      REAL(kind=8), dimension(lrw) :: RWORK
      REAL(KIND=8), dimension(NSMAX) :: DUMMYPDJ, DUMMYY
      integer IDUMMY
      integer :: i
      integer :: itol, itask, istate, iopt, mf
      real(kind=8) :: atol
      real(kind=8), dimension(nsmax) :: satol

      real(kind=8) :: T, TOUT, TIN

      integer :: NNZ

      EXTERNAL FCHEM,JAC


! Initialize work arrays

      iwork(:) = 0
      rwork(:) = 0.
      IWORK(5) = 5
      RWORK(6) = 3.154D14
      IWORK(6) = 10000
      IWORK(7) = 2

      if (it.eq.1) IWORK(6)=2000

! Changing the time to avoid T + DT = T 

      TOUT=TIME

      TIN = T
      TOUT = TOUT - TIN
      T = 0.d0      

      Y(:) = XN(:)

! If testjac is 1, print non zero elements per column of the Jacobian
! Done in odes/JACVW

      if (TESTJAC.eq.1) then
      DUMMYY=1.d-5
      call ratcon(DUMMYY)
      do IDUMMY=1,NSMAX
      call JACVW(DUMMYY,IDUMMY,DUMMYPDJ)
      enddo
      STOP
      endif

      call shieldingsetup

! Don't stop running, Forrest

      do while (t.lt.tout)

      istate = 1

! Adaptive absolute tolerance to avoid too high precision on very abundant species,
! H2 for instance. Helps running a bit faster

      do i=1,nsmax
      satol(i)=max(atol,1.d-16*y(i))
      enddo

! ratcond is already called in compute IAJA
!      call ratcon(Y)

! Feed IWORK with IA and JA

      call computeIAJA(Y)
      NNZ=IA(NSMAX+1)-1
      iwork(30+1:30+NSMAX+1)=IA(1:NSMAX+1)
      iwork(31+NSMAX+1:31+NSMAX+NNZ)=JA(1:NNZ)

      call dlsodes(fchem,nsmax,y,t,tout,itol,rtol,satol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)       

! Whenever the solver fails converging, print the reason.
! cf odpkdmain.f for translation
      if (istate.ne.2) print *, 'IPTS = ', ipts, 'ISTATE = ', ISTATE

      CALL CONSERVE(Y)

! Stop, Forrest
      enddo

      T = TOUT + TIN 

      XN=Y

      RETURN 
      END 
 
! ======================================================================
! ======================================================================
      SUBROUTINE FCHEM (N,T,Y,YP)
      use header
      use unitvar
      IMPLICIT none

      REAL(kind=8), dimension(NSMAX) :: Y,YP
      real(kind=8) :: T
      INTEGER :: N

      CALL RATCON2(Y)
      CALL FCHEMVW(N,Y,YP) 

      RETURN
      END


! ======================================================================
! ======================================================================
      subroutine JAC(N, T, Y, J, IAN, JAN, PDJ)
      use header
      implicit none
      integer N,J
      real(kind=8) :: T
      REAL(KIND=8), dimension(N) :: IAN, JAN
      REAL(kind=8), dimension(NSMAX) :: Y,PDJ

!      CALL RATCON2(Y)

      call JACVW(Y,J,PDJ) 

      return
      end

! ======================================================================
! ======================================================================
      SUBROUTINE CONSERVE(Y)
      use header
      use unitvar
      IMPLICIT none

      REAL(kind=8), dimension(NSMAX) :: Y
      real(kind=8), dimension(NEMAX) :: ELMSUM
      real(kind=8) :: CHASUM

      integer :: i, k

! Prevent too low abundances
      do i=1,NSMAX
      if (Y(i).le.1.d-99) Y(i)=1.d-99
      enddo

! --- Conserve electrons
      CHASUM=0.0D+0
      DO I=1,NSMAX
         IF (I.NE.ISPE) CHASUM=CHASUM+ICG(I)*Y(I)
      ENDDO
      IF (CHASUM.LE.0.0D+0) CHASUM=XNMIN
      Y(ISPE)=CHASUM

! --- Conserve other elements if selected
      IF (ICONS.GT.0) THEN
         DO K=1,ICONS
            ELMSUM(K)=0.0D+0
         ENDDO
         DO I=1,NSMAX
            DO K=1,ICONS
               IF (I.NE.ISPELM(K)) ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
            ENDDO
         ENDDO
         DO K=1,ICONS
            Y(ISPELM(K))=ELEMS(K)-ELMSUM(K)
            IF (Y(ISPELM(K)).LE.0.0D+0) Y(ISPELM(K))=XNMIN
         ENDDO
      ENDIF

! Check for conservation
       ELMSUM(:)=0.0D+0
       DO I=1,NSMAX
          DO K=1,NEMAX
             ELMSUM(K)=ELMSUM(K)+IELM(K,I)*Y(I)
          ENDDO
       ENDDO

       do k=1,nemax
       if (abs(ELEMS(K)-ELMSUM(K))/ELEMS(K).ge.0.01) then 
       print *, 'CAUTION : Element ',SPEC(ISPELM(K)), 'is not conserved'
       print *, 'Relative difference: ', abs(ELEMS(K)-ELMSUM(K))/ELEMS(K)
       endif
       if (SPEC(ISPELM(K)).eq.YH) then
                     if (abs(ELEMS(K)-Y(INDH2)*2.D0)/ELEMS(K).ge.0.01) print *,'H is too depleted on the grains !!!!'
       endif
       if (SPEC(ISPELM(K)).eq.YHE) then
                     if (abs(ELEMS(K)-Y(INDHE))/ELEMS(K).ge.0.01) print *,'He is too depleted on the grains !!!!'
       endif       
       enddo

! VW fev 2012 add a test for the helium and H2 abundance in the gas phase
! prevent excessive depletion
       

      RETURN
      END

! ======================================================================
! Dummy jacobians when the code is run with mf=222
! Not to use unless, the solver has big problems converging
! ======================================================================
!      SUBROUTINE DUMMY 
!      IMPLICIT none
!      integer :: N,J
!      real(kind=8) :: T,Y,IAN, JAN, PDJ 

!      entry jac (N,T,Y,ML,MU,PD,NROWPD)
!      entry jac (N, T, Y, J, IAN, JAN, PDJ)
!      return

!      END
 
! ======================================================================
! ======================================================================
      SUBROUTINE MODIF(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
      use header
      use unitvar
      IMPLICIT none

      integer :: J,IMOD1,IMOD2
      real(kind=8) :: BARR,YMOD1,YMOD2, PICK, TESTREF1, TESTREF2, TESTNUM

      EX1(J)=0.0D+0
      EX2(J)=0.0D+0

! --- Check value of x = t_acc/t_evap
! TINEVA = 1/t_evap
! TINACC = 1/t_acc
      IF (TINACC(JSP1(J)).GT.0.0D+0) THEN
         EX1(J)=TINEVA(JSP1(J))/TINACC(JSP1(J))
      ENDIF
      IF (TINACC(JSP2(J)).GT.0.0D+0) THEN
         EX2(J)=TINEVA(JSP2(J))/TINACC(JSP2(J))
      ENDIF
! Hence x = 0 if t_evap or t_acc = 0

! --- Assign max rates

      IF (BARR.EQ.1.0D+0) THEN
       IF (IMOD1.NE.0) THEN
          IF (EX1(J).LT.1.0D+0) THEN
             IF (RDIF1(J).GT.TINACC(JSP1(J))) RDIF1(J)=TINACC(JSP1(J))
          ENDIF
          IF (EX1(J).GE.1.0D+0) THEN
             IF (RDIF1(J).GT.TINEVA(JSP1(J))) RDIF1(J)=TINEVA(JSP1(J))
          ENDIF
       ENDIF

       IF (IMOD2.NE.0) THEN
          IF (EX2(J).LT.1.0D+0) THEN
             IF (RDIF2(J).GT.TINACC(JSP2(J))) RDIF2(J)=TINACC(JSP2(J))
          ENDIF
          IF (EX2(J).GE.1.0D+0) THEN
             IF (RDIF2(J).GT.TINEVA(JSP2(J))) RDIF2(J)=TINEVA(JSP2(J))
          ENDIF
       ENDIF
      ENDIF

! --- Species rate to compare chosen by fastest diffusion rate
      IF (BARR.NE.1.0D+0) THEN
      PICK=0

       TESTREF1=TINACC(JSP1(J))
       IF (EX1(J).GE.1.0D+0) TESTREF1=TINEVA(JSP1(J))
       TESTREF2=TINACC(JSP2(J))
       IF (EX2(J).GE.1.0D+0) TESTREF2=TINEVA(JSP2(J))

       IF (RDIF1(J).GE.RDIF2(J)) THEN
         TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD2*GTODN
         IF (YMOD2*GTODN.LT.1.0D+0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
         IF (TESTNUM.GT.TESTREF1) PICK=1
       ENDIF
       IF (RDIF2(J).GT.RDIF1(J)) THEN
         TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD1*GTODN
         IF (YMOD1*GTODN.LT.1.0D+0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
         IF (TESTNUM.GT.TESTREF2) PICK=2
       ENDIF

       IF (PICK.EQ.1) THEN
                RDIF1(J)=TESTREF1/BARR/YMOD2/GTODN
                IF (YMOD2*GTODN.LT.1.0D+0) RDIF1(J)=TESTREF1/BARR
                RDIF2(J)=0.0D+0
       ENDIF

       IF (PICK.EQ.2) THEN
               RDIF2(J)=TESTREF2/BARR/YMOD1/GTODN
               IF (YMOD1*GTODN.LT.1.0D+0) RDIF2(J)=TESTREF2/BARR
               RDIF1(J)=0.0D+0
       ENDIF

      ENDIF

      RETURN
      END

! ======================================================================
! ======================================================================
      SUBROUTINE RATCON(Y)

! Reactions coefficient formally dependent on the abundances Y are 
! computed in a companion subroutine: RATCON2
! Grain surface reactions, self shielding, etc...
! VW modification of everything for the new gas-phase network
! Fev 2012

      use header 
      use unitvar
      IMPLICIT none

      real(kind=8) :: T300, TI, TSQ
      REAL(kind=8), dimension(NSMAX) :: Y
      integer :: nsta, nfin
      integer :: k, j, w, m, n
      integer, dimension(10) :: indice
      real(kind=8), dimension(10) :: distmin, distmax

      T300=TEMP/300.D+0
      TI=1.0D+00/TEMP
      TSQ=SQRT(TEMP)
       
! ====== Rxn ITYPE 0
! ITYPE 0: Gas phase reactions with GRAINS =) 
         DO J=IRXSTA(0),IRXFIN(0)
            XK(J)=A(J)*(T300**B(J))
         ENDDO
        
! In the case idust eq 0, we still need the formation of H on the grains
! this is done with the XH species, reaction types 10 and 11

! ====== Rxn ITYPE 10 and 11
! ITYPE 10 and 11: H2 formation on the grains when idust eq 0
       if (idust.eq.0) then
               DO J=IRXSTA(10),IRXFIN(10)
                       XK(J)=A(J)*1.186D7*exp(225.D0/TEMP)**(-1)*GTODN/XNT
!                     print *,'GTODN', GTODN
              enddo
              DO J=IRXSTA(11),IRXFIN(11)       
                              XK(J)=A(J)*(T300**B(J))*XNT/GTODN
              enddo       
       endif
       

! ====== Rxn ITYPE 1
! ITYPE 1: Photodissoc/ionisation with cosmic rays
! Add X-rays in case ZETAX is not 0
        DO J=IRXSTA(1),IRXFIN(1)
            XK(J)=A(J)*(ZETA0+ZETAX)
       ENDDO
        DO J=IRXSTA(2),IRXFIN(2)
            XK(J)=A(J)*(ZETA0+ZETAX)
       ENDDO

! ====== Rxns ITYPE 4 - 8
! Bimolecular gas phase reactions - several possible formula 
         W=1
         NSTA=0
         NFIN=0
         distmin(:)=9999.
         distmax(:)=9999.
         DO J=4,8
                IF ((IRXSTA(J).NE.0).AND.(NSTA.EQ.0)) NSTA=J
                IF (IRXFIN(J).NE.0) NFIN=J
         ENDDO
         DO J=IRXSTA(NSTA),IRXFIN(NFIN)
        
!---------------  KOOIJ FORMULA
            IF (FORMULA(J).eq.3) THEN
              
               XK(J)=A(J)*(T300**B(J))*EXP(-C(J)/TEMP)

!              Check for temperature bounderies
               if (TEMP.LT.Tmin(J)) XK(J)=A(J)*((Tmin(J)/300.D0)**B(J))*EXP(-C(J)/Tmin(J))
               if (TEMP.GT.Tmax(J)) XK(J)=A(J)*((Tmax(J)/300.D0)**B(J))*EXP(-C(J)/Tmax(J))
                 
!              Check for the presence of several rate coefficients present in the network for the
!              the same reaction
               IF (NUM(J+1).EQ.NUM(J)) THEN
                  INDICE(W)=J
                  distmin(w)=tmin(j)-temp
                  distmax(w)=temp-tmax(j)
                  W = W + 1
               ENDIF

               IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

                  INDICE(W)=J
                  distmin(w)=tmin(j)-temp
                  distmax(w)=temp-tmax(j)

                  DO M=1,W
                     N=INDICE(M)
                     !IF(IT==1) PRINT*,N,M, SYMBOL(:,N), tmin(N), tmax(N),distmin(M),distmax(M)
                     IF (TEMP.LT.Tmin(N)) XK(N)=0.d0
                     IF (TEMP.GT.Tmax(N)) XK(N)=0.d0
                  ENDDO

                  IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN

                     IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
                        N=indice(minloc(abs(distmin),dim=1))
                        XK(N)=A(N)*((Tmin(N)/300.D0)**B(N))*EXP(-C(N)/Tmin(N))
                     ELSE
                        N=indice(minloc(abs(distmax),dim=1))
                        XK(N)=A(N)*((Tmax(N)/300.D0)**B(N))*EXP(-C(N)/Tmax(N))
                     ENDIF
                  ENDIF

                  W=1
                  INDICE(:)=0
                  distmin(:)=9999.
                  distmax(:)=9999.
               ENDIF
            ENDIF
                     
!---------------  IONPOL1 FORMULA
            IF (FORMULA(J).EQ.4) THEN
               XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TEMP)**0.5))

!              Check for temperature bounderies
               IF (TEMP.LT.Tmin(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/Tmin(J))**0.5))
               IF (TEMP.GT.Tmax(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TMAX(J))**0.5))
                 
!              Check for the presence of several rate coefficients present in the network for the
!              the same reaction
               IF (NUM(J+1).EQ.NUM(J)) THEN
                  INDICE(W)=J
                  distmin(w)=tmin(j)-temp
                  distmax(w)=temp-tmax(j)
                  W = W + 1
               ENDIF

               IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

                  INDICE(W)=J
                  distmin(w)=tmin(j)-temp
                  distmax(w)=temp-tmax(j)

                  DO M=1,W
                     N=INDICE(M)
                     IF (TEMP.LT.Tmin(N)) XK(N)= 0.d0
                     IF (TEMP.GT.Tmax(N)) XK(N)= 0.d0
                  ENDDO
              
                  IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN
                     IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
                        N=indice(minloc(abs(distmin),dim=1))
                        XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/Tmin(N))**0.5))
                     ELSE
                        N=indice(minloc(abs(distmax),dim=1))
                        XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/TMAX(N))**0.5))
                     ENDIF
                  ENDIF

                  W=1
                  INDICE(:)=0
                  distmin(:)=9999.
                  distmax(:)=9999.
                ENDIF
             ENDIF

!---------------  IONPOL2 FORMULA
             IF (FORMULA(J).EQ.5) then
                XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TEMP)**0.5)+(C(J)**2*300.d0/(10.526*TEMP)))

!               Check for temperature bounderies
                IF (TEMP.LT.Tmin(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMIN(J))**0.5)+(C(J)**2*300.d0/(10.526*TMIN(J))))
                IF (TEMP.GT.Tmax(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMAX(J))**0.5)+(C(J)**2*300.d0/(10.526*TMAX(J))))
                 
!               Check for the presence of several rate coefficients present in the network for the
!               the same reaction
                IF (NUM(J+1).EQ.NUM(J)) THEN
                   INDICE(W)=J
                   distmin(w)=tmin(j)-temp
                   distmax(w)=temp-tmax(j)
                   W = W + 1
                ENDIF

                IF ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) THEN

                   INDICE(W)=J
                   distmin(w)=tmin(j)-temp
                   distmax(w)=temp-tmax(j)

                   DO M=1,W
                      N=INDICE(M)
                      IF (TEMP.LT.Tmin(N)) XK(N)=       0.d0
                      IF (TEMP.GT.Tmax(N)) XK(N)=       0.d0
                   ENDDO

                   IF (maxval(XK(indice(1:w))).lt.1.d-99) THEN
                      IF (minval(abs(distmin)).lt.minval(abs(distmax))) THEN
                         N=indice(minloc(abs(distmin),dim=1))
                         XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMIN(N))**0.5)+(C(N)**2*300.d0/(10.526*TMIN(N))))
                      ELSE
                         N=indice(minloc(abs(distmax),dim=1))
                         XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMAX(N))**0.5)+(C(N)**2*300.d0/(10.526*TMAX(N))))
                      ENDIF
                   ENDIF

                   W=1
                   INDICE(:)=0
                   distmin(:)=9999.
                   distmax(:)=9999.
                ENDIF
             ENDIF
         ENDDO

! === Grain surfaces
         IF (IDUST.NE.0) THEN

! ========= Set diffusion and evaporation rates (s-1)
            DO K=1,NSMAX
               TINDIF(K)=CHF(K)*EXP(-EB(K)/DTEMP)/TNS
               TINEVA(K)=CHF(K)*EXP(-ED(K)/DTEMP)
            ENDDO

! ========= Rxn ITYPE 15 - thermal evaporation
! ITYPE 15: Thermal evaporation
            DO J=IRXSTA(15),IRXFIN(15)
               XK(J)=A(J)*XJ(J)*TINEVA(JSP1(J))
            ENDDO

! ========= Rxn ITYPE 16
! ITYPE 16: Cosmic-ray evaporation
            DO J=IRXSTA(16),IRXFIN(16)
               XK(J)=A(J)*XJ(J)*((ZETA0+ZETAX)/1.3D-17)&
                        *CHF(JSP1(J))*CRFE*CRT*EXP(-ED(JSP1(J))/TSMAX)
            ENDDO
        

! Photodesorption, when used appears through ITYPES 66 and 67
      if (irxsta(66).ne.0) then
! ========= Rxn ITYPE 66
! ITYPE 66: CO photodesorption by external UV
! 1.d8 is I_ISRF-FUV from Oberg et al. 2007, ApJ, 662, 23
      do j=irxsta(66),irxfin(66)
      XK(J)=A(J)/SNS*UVGAS*1.d8*exp(-2.*TAU) 
      enddo

! ========= Rxn ITYPE 67
! ITYPE 67: CO photodesorption by CR generated UV
      do j=irxsta(67),irxfin(67)
      XK(J)=A(J)/SNS*1.d4
      enddo
     
      endif

      if (irxsta(98).ne.0) then
! ========= Rxn ITYPE 98 test the storage of H2S under a refractory form
! ITYPE 98: storage of H2S under a refractory form
      do j=irxsta(98),irxfin(98)
      XK(J)=A(J)*(T300**B(J))*EXP(-C(J)*TI)
      enddo
      endif
      
     
! ====== Rxn ITYPE 17
! ITYPE 17: Photodissociations by Cosmic rays on grain surfaces
! Add X-rays in case ZETAX is not 0
         DO J=IRXSTA(17),IRXFIN(17)
            XK(J)=A(J)*(ZETA0 + ZETAX)
!            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
        ENDDO

! ====== Rxn ITYPE 18
! ITYPE 18: Photodissociations by Cosmic rays on grain surfaces
         DO J=IRXSTA(18),IRXFIN(18)
            XK(J)=A(J)*(ZETA0+ZETAX)
!            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
        ENDDO
       
        ENDIF

! When dust is turned off, zero all dust rates==========================
      IF ((IDUST.EQ.0).AND.(IT.EQ.1)) THEN
         DO J=IRXSTA(14),IRXFIN(99)
            XK(J)=0.0D+0
            XJ(J)=0.0D+0
         ENDDO
      ENDIF


      RETURN
      END

! ======================================================================
! Computes the H2 and CO column density for their self-shielding
! for each ipts, the column density is incremented recursively
! from the results for NH2 and NCO computed by RATCON2 the spatial step before
! NB: ZN is shifted with respect to N
! ======================================================================
      SUBROUTINE SHIELDINGSETUP
 
      use header
      implicit none
      real(kind=8) :: XNH2,XNCO
      
      if (iptstore.eq.1) then
      ZNH2(iptstore)=0.
      ZNCO(iptstore)=0.
      else

      XNH2=ZXN(indH2,iptstore-1)
      XNCO=ZXN(indCO,iptstore-1)

      ZNH2(iptstore)=ZNH2(iptstore-1)+XNT*zstepsize*XNH2
      ZNCO(iptstore)=ZNCO(iptstore-1)+XNT*zstepsize*XNCO
      endif

      return
      end

! ======================================================================
! ======================================================================
      SUBROUTINE RATCON2(Y)

      use header 
      use unitvar
      IMPLICIT none

      REAL(kind=8) :: ACTIV,BARR,MONLAY,DIFF
      REAL(kind=8) :: XNH2,XNCO
      REAL(kind=8) :: TETABIS,TETABIS1,TETABIS2,TETABIS3
      real(kind=8) :: T300, TI, TSQ
      REAL(kind=8), dimension(NSMAX) :: Y
      real(kind=8) :: YMOD1, YMOD2
      INTEGER IMOD1,IMOD2
      integer :: j, l

      T300=TEMP/300.D+0
      TI=1.0D+00/TEMP
      TSQ=SQRT(TEMP)
      MONLAY=LAYERS*TNS/GTODN

      XNH2=Y(indH2)
      XNCO=Y(indCO)

! Density/Av-dependent==================================================
!
! VW 02/07 the treatment for the H2 and CO photodissociation was corrected
! for larger and smaller values of H2, CO column densities and Av not in the 
! data tables. I also added the choice of using this approximation in gg_control.d
!
! ====== Rxn ITYPE 13
! ITYPE 2: Gas phase photodissociations/ionisations by UV
         DO J=IRXSTA(3),IRXFIN(3)
            XK(J)=A(J)*EXP(-C(J)*TAU)*UVGAS

! MODIFY THE H2 AND CO PHOTODISSOCIATION IF ISABS EQ 1
         IF (ISABS.EQ.1) THEN 
    
! ====== Compute the H2 self-shielding
         IF (SYMBOL(1,J).EQ.YH2) THEN
                TETABIS=1.D0
         if (iptstore.eq.1) then
                ZNH2(iptstore) = TAU/RAVNH * XNH2
         endif
         NH2=ZNH2(iptstore)

! ======= Linear extrapolation of the shielding factors
            DO L=1,NL1-1
               IF ((N1H2(L).LE.NH2).AND.(N1H2(L+1).GE.NH2)) THEN
               TETABIS=T1H2(L)+(NH2-N1H2(L))*(T1H2(L+1)-T1H2(L))/(N1H2(L+1)-N1H2(L))
               ENDIF
            ENDDO
               IF (NH2.GT.N1H2(NL1)) TETABIS = T1H2(NL1)

          XK(J)=2.54D-11*TETABIS

          XK(J)=XK(J)*UVGAS
          ENDIF

! ====== Compute the CO self-shielding
          IF (SYMBOL(1,J).EQ.YCO) THEN

          TETABIS1=1.D0
          TETABIS2=1.D0
          TETABIS3=1.D0

          if (iptstore.eq.1) then
             ZNH2(iptstore) = TAU/RAVNH * XNH2
             ZNCO(iptstore) = TAU/RAVNH * XNCO
          endif
          NH2=ZNH2(iptstore)
          NCO=ZNCO(iptstore)

! ======= Linear extrapolation of the three shileding factors
           DO L=1,NL2-1
               IF ((N2CO(L).LE.NCO).AND.(N2CO(L+1).GE.NCO))  &
               TETABIS2=T2CO(L)+(NCO-N2CO(L))*(T2CO(L+1)-T2CO(L))&
               /(N2CO(L+1)-N2CO(L))
           ENDDO

           DO L=1,NL3-1
               IF ((N2H2(L).LE.NH2).AND.(N2H2(L+1).GE.NH2)) &
               TETABIS1=T2H2(L)+(NH2-N2H2(L))*(T2H2(L+1)-T2H2(L))&
               /(N2H2(L+1)-N2H2(L))
               IF ((AV2(L).LE.TAU).AND.(AV2(L+1).GE.TAU)) &
               TETABIS3=T2AV(L)+(TAU-AV2(L))*(T2AV(L+1)-T2AV(L))&
               /(AV2(L+1)-AV2(L))
           ENDDO

! Saturate the rate coefficient if necessary (when density or Av are out of 
! the shielding array, from the photodiss files)

           IF (NCO.GT.N2CO(NL2)) TETABIS2 = T2CO(NL2)
           IF (NH2.GT.N2H2(NL3)) TETABIS1 = T2H2(NL3)
           IF (TAU.GT.AV2(NL3))  TETABIS3 = T2AV(NL3)

           XK(J)=1.03D-10*TETABIS1*TETABIS2*TETABIS3
           XK(J)=XK(J)*UVGAS

            ENDIF
       
         ENDIF

         ENDDO

! Continually time-dependent grain rates================================
      IF (IDUST.NE.0) THEN

! ====== Rxn ITYPE 99
! ITYPE 99: Adsorption on grains
         DO J=IRXSTA(99),IRXFIN(99)
! ========= Set accretion rates
            TINACC(JSP1(J))=CONDSP(JSP1(J))*TSQ*Y(JSP1(J))*XNT
            TINACC(JSP2(J))=TINACC(JSP1(J))
            XK(J)=A(J)*XJ(J)*TINACC(JSP1(J))/Y(JSP1(J))/GTODN
        ENDDO

! ====== Rxn ITYPE 14
! ITYPE 14: Grain surface reactions
         DO J=IRXSTA(14),IRXFIN(14)
            IMOD1=0
            IMOD2=0
            BARR=1.0D+0
! --------- Calculate activation energy barrier multiplier
            IF (EA(J).GE.1.0D-40) THEN
               ACTIV=EA(J)/DTEMP
! ------------ Choose fastest of classical or tunnelling
               IF (ACTIV.GT.ACT1(J)) ACTIV=ACT1(J)
               BARR=EXP(-ACTIV)
            ENDIF

! --------- Thermal hopping diffusion method
            RDIF1(J)=TINDIF(JSP1(J))
            RDIF2(J)=TINDIF(JSP2(J))

! --------- Check for JH,JH2
            IF (SYMBOL(1,J).EQ.YJH)  IMOD1=1
            IF (SYMBOL(1,J).EQ.YJH2) IMOD1=2
            IF (SYMBOL(2,J).EQ.YJH)  IMOD2=1
            IF (SYMBOL(2,J).EQ.YJH2) IMOD2=2

! --------- QM for JH,JH2 only - others are too heavy
            IF (IMOD1+IMOD2.NE.0) THEN
! ------------ QM1 - Tunnelling (if it's faster than thermal)
               IF (IGRQM.EQ.1) THEN
                  IF ((IMOD1.NE.0).AND.&
                    (RQ1(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ1(JSP1(J))
                  IF ((IMOD2.NE.0).AND.&
                    (RQ1(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ1(JSP2(J))
               ENDIF
! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
               IF (IGRQM.EQ.2) THEN
                  IF ((IMOD1.NE.0).AND.&
                      (RQ2(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ2(JSP1(J))
                  IF ((IMOD2.NE.0).AND.&
                      (RQ2(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ2(JSP2(J))
               ENDIF
! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
               IF (IGRQM.EQ.3) THEN
                  IF (IMOD1.NE.0) THEN
                     IF (RQ1(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ1(JSP1(J))
                     IF (RQ2(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ2(JSP1(J))
                  ENDIF
                  IF (IMOD2.NE.0) THEN
                     IF (RQ1(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ1(JSP2(J))
                     IF (RQ2(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ2(JSP2(J))
                  ENDIF
               ENDIF
            ENDIF

! --------- Modify according to IMODH switch:
            IF (IMODH.NE.0) THEN
! ------------ If H+H->H2 is only modified rxn:
               IF ((IMODH.EQ.-1).AND.(IMOD1.NE.1.OR.IMOD2.NE.1)) THEN
                  IMOD1=0
                  IMOD2=0
               ENDIF

! ------------ If only H is modified:
               IF ((IMODH.EQ.1).AND.(IMOD1.NE.1)) IMOD1=0
               IF ((IMODH.EQ.1).AND.(IMOD2.NE.1)) IMOD2=0

! ------------ Set to modify all rates, if selected (just atoms)
               IF (IMODH.EQ.3) THEN
                       IF ((SYMBOL(1,J).EQ.YJH).OR.&
                              (SYMBOL(1,J).EQ.'JHe        ').OR.&
                              (SYMBOL(1,J).EQ.'JC         ').OR.&
                              (SYMBOL(1,J).EQ.'JN         ').OR.&
                              (SYMBOL(1,J).EQ.'JO         ').OR.&
                              (SYMBOL(1,J).EQ.'JS         ').OR.&
                              (SYMBOL(1,J).EQ.'JSi        ').OR.&
                              (SYMBOL(1,J).EQ.'JFe        ').OR.&
                              (SYMBOL(1,J).EQ.'JNa        ').OR.&
                              (SYMBOL(1,J).EQ.'JMg        ').OR.&
                              (SYMBOL(1,J).EQ.'JP         ').OR.&
                              (SYMBOL(1,J).EQ.'JCl        ')) IMOD1=3
                       IF ((SYMBOL(2,J).EQ.YJH).OR.&
                              (SYMBOL(2,J).EQ.'JHe        ').OR.&
                              (SYMBOL(2,J).EQ.'JC         ').OR.&
                              (SYMBOL(2,J).EQ.'JN         ').OR.&
                              (SYMBOL(2,J).EQ.'JO         ').OR.&
                              (SYMBOL(2,J).EQ.'JS         ').OR.&
                              (SYMBOL(2,J).EQ.'JSi        ').OR.&
                              (SYMBOL(2,J).EQ.'JFe        ').OR.&
                              (SYMBOL(2,J).EQ.'JNa        ').OR.&
                              (SYMBOL(2,J).EQ.'JMg        ').OR.&
                              (SYMBOL(2,J).EQ.'JP         ').OR.&
                              (SYMBOL(2,J).EQ.'JCl        ')) IMOD2=3
               ENDIF

! ------------ Modify rates (RDIF1 & RDIF2) according to their own evap/acc rates
                   YMOD1=Y(JSP1(J))
                   YMOD2=Y(JSP2(J))
               CALL MODIF(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
            ENDIF

            DIFF=RDIF1(J)+RDIF2(J)

            XK(J)=A(J)*XJ(J)*BARR*DIFF*GTODN/XNT
!              XK(J)=0.D0
! --------- Allow only 1 monolayer of each to react
! Not used for the time being
!            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*TNS/GTODN/Y(JSP1(J))
!            IF (Y(JSP2(J)).GT.MONLAY) XK(J)=XK(J)*TNS/GTODN/Y(JSP2(J))

         ENDDO

! ====== Rxn ITYPE 19 - 20
! ITYPE 19: Photodissociations by UV photons on grain surfaces
! ITYPE 20: Photodissociations by UV photons on grain surfaces
         DO J=IRXSTA(19),IRXFIN(20)
            XK(J)=A(J)*EXP(-C(J)*TAU)*UVGAS
!            IF (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
         ENDDO

! Useful for testing
! To disable some reaction types
!do j=irxsta(14),irxfin(15)
!xk(j)=0.
!enddo

!do j=irxsta(17),irxfin(20)
!xk(j)=0.
!enddo
      ENDIF

! Continually time-dependent gas phase rates============================
! H2 formation
! XJ(1) and XJ(2) are zero if IDUST=1
! cf GRAINRATE
! VW Fev 2012 - this process has been removed
!      do j=irxsta(0),irxfin(0)
!      if ((SYMBOL(1,J).eq.YH).and.(SYMBOL(2,j).eq.YH)) then
!      XK(j)=XJ(j)*A(j)*(T300**B(j))*GTODN/XNT/Y(JSP1(j))
!      endif
!      enddo

! If rate acoefficients are too small, put them to 0
      where (XK.lt.RXNMIN) XK=0.d0

      RETURN
      END

! ======================================================================
! ======================================================================
      SUBROUTINE GRAINRATE
      use header
      use unitvar
      IMPLICIT none

      REAL(kind=8), dimension(nsmax) :: REA1,REA2,REA3,REA4
      REAL(kind=8), dimension(nkmax) :: REA5
      REAL(kind=8), dimension(nsmax) :: SMASS
      real(kind=8) :: SMA,REDMAS,STICK,EVFRAC,DHFSUM,SUM1,SUM2
      INTEGER, dimension(NKMAX) :: INT1
      integer :: NGS,NEA,NPATH,NEVAP,BADFLAG,ATOMS
      CHARACTER(len=11), dimension(5,NKMAX) :: GSREAD
      character(len=11), dimension(nsmax) :: GSPEC

      real(kind=8) :: cond
      integer :: i, j,k,l,n4, n5, n6

! Set accretion rate====================================================

! COND is used to calculate R_acc = (sigma_d) * <v_i> * n_i * n_d
! Uses 'Mean' Maxwellian speed, rather than RMS (is this justifiable?)

      COND=PI*RD*RD*SQRT(8.0D+0*BOLTZ/PI/AMU)

! --- Evaluate sticking coeff and accretion rate factor for each species
      STICK=0.0D+0
      DO I=1,NSMAX
         IF (ICG(I).EQ.0) THEN
            STICK=STICK0
         ENDIF 
          IF (ICG(I).GT.0) THEN 
            STICK=STICKP 
         ENDIF
          IF (ICG(I).LT.0) THEN 
            STICK=STICKN 
         ENDIF
!         IF (SPEC(I).EQ.YH2)      STICK=0.D+0 
!         IF (SPEC(I).EQ.YHE)      STICK=0.D+0 
!         IF (SPEC(I).EQ.YH)       STICK=0.D+0 
         IF (SPEC(I).EQ.YHEP)     STICK=0.D+0 
         IF (SPEC(I).EQ.'e-         ')      STICK=0.D+0
         IF (SPEC(I).EQ.'H+         ')     STICK=0.D+0
         IF (SPEC(I).EQ.YGRAIN)   STICK=0.D+0
         IF (SPEC(I).EQ.'GRAIN-     ') STICK=0.D+0
!         IF (SPEC(I).EQ.'H-')     STICK=0.D+0
!         IF (SPEC(I).EQ.'H2+')    STICK=0.D+0

         IF (I.GT.NSGAS) STICK=0.D+0
         CONDSP(I)=COND*STICK/SQRT(AWT(I))
      ENDDO

! Read in molecular information for surface rates=======================
      READ(NGR,*)

! --- Read info into dummy arrays
      NGS=0
      DO I=1,NSMAX
  700    CONTINUE
         READ(NGR,705) GSPEC(I),INT1(I),REA1(I),REA2(I),REA3(I),REA4(I)
  705    FORMAT(A11,I4,F7.0,F6.0,D8.1,27X,F8.2)
         IF (GSPEC(I).EQ.'X          ') GO TO 700
         IF (GSPEC(I).EQ.'           ') EXIT
         NGS=NGS+1
      ENDDO
       

! --- Read activation energies into dummy arrays
      READ(NGR,720) NEA
  720 FORMAT(I4)
      DO J=1,NEA
        READ(NGR,730) (GSREAD(L,J),L=1,5),REA5(J)
  730   FORMAT(5A11,D9.2)
      ENDDO

! --- Transfer from dummies to arrays with correct species numbers
      DO I=1,NSMAX
         SMASS(I)=0
         ED(I)=0.0D+0
         EB(I)=0.0D+0
         DEB(I)=0.0D+0
         DHF(I)=0.0D+0
         DO J=1,NGS
            IF (SPEC(I).EQ.GSPEC(J)) THEN
               SMASS(I)=INT1(J)
               ED(I)=REA1(J)
               EB(I)=REA2(J)
               DEB(I)=REA3(J)
               DHF(I)=REA4(J)
               IF ((SPEC(I).NE.YJH).AND.(SPEC(I).NE.YJH2).AND.&
                  (EBFAC.GE.0.0D+0)) EB(I)=EBFAC*ED(I)
            ENDIF
         ENDDO
         !IF(SPEC(I) == 'JN2O2      ') PRINT*,ED(I)
      ENDDO

      DO I=1,NKMAX
         EA(I)=0.0D+0
         DO J=1,NEA
            IF (SYMBOL(4,I)(:1).EQ.'J          ') THEN
               IF ((SYMBOL(1,I).EQ.GSREAD(1,J)).AND.&
                  (SYMBOL(2,I).EQ.GSREAD(2,J)).AND.&
                  (SYMBOL(4,I).EQ.GSREAD(3,J)).AND.&
                  (SYMBOL(5,I).EQ.GSREAD(4,J)).AND.&
                  (SYMBOL(6,I).EQ.GSREAD(5,J))) EA(I)=REA5(J)
            ENDIF
            IF (SYMBOL(4,I)(:1).NE.'J          ') THEN
               IF ((SYMBOL(1,I).EQ.GSREAD(1,J)).AND.&
                  (SYMBOL(2,I).EQ.GSREAD(2,J)).AND.&
                  (SYMBOL(4,I).EQ.GSREAD(3,J)(2:)).AND.&
                  (SYMBOL(5,I).EQ.GSREAD(4,J)(2:)).AND.&
                  (SYMBOL(6,I).EQ.GSREAD(5,J)(2:))) EA(I)=REA5(J)
            ENDIF
         ENDDO
	 !IF(symbol(4,i) == 'JO2H       ') PRINT*, symbol(:,i), Ea(i)
      ENDDO

! Set up constants, quantum rate info===================================
      DO I=1,NSMAX
         CHF(I)=0.0D+0
         RQ1(I)=0.0D+0
         RQ2(I)=0.0D+0
! ------ For species which have been assigned surface info, SMASS=/=0
         IF (SMASS(I).NE.0) THEN
            SMA=REAL(SMASS(I))
! --------- Set characteristic frequency
            CHF(I)=SQRT(2.0D+0*BOLTZ/PI/PI/AMU * SNS*ED(I)/SMA)
! --------- Set quantum rates
            IF (DEB(I).GE.1.0D-38) THEN
               RQ1(I)=DEB(I)*BOLTZ/4.0D+0/HBAR/TNS
            ELSE
               RQ1(I)=0.0D+0
            ENDIF
            RQ2(I)=CHF(I)/TNS*&
                 EXP(-2.0D+0*ACM/HBAR*SQRT(2.0D+0*AMU*SMA*BOLTZ*EB(I)))
         ENDIF
      ENDDO

! === Cycle all reactions
      DO J=1,NKMAX

! ------ Initialise all XJ rate factors, and get species 1 & 2
         XJ(J)=1.0D+0
         JSP1(J)=0
         JSP2(J)=0
         DO I=1,NSMAX
            IF (SYMBOL(1,J).EQ.SPEC(I)) JSP1(J)=I
            IF (SYMBOL(2,J).EQ.SPEC(I)) JSP2(J)=I
         ENDDO

! === ITYPE 14 - SURFACE REACTIONS
      IF (ITYPE(J).EQ.14) THEN
         NPATH=0

! ------ Check for branching
         DO K=1,NKMAX
            IF (((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
                (SYMBOL(2,J).EQ.SYMBOL(2,K))).OR.&
               ((SYMBOL(2,J).EQ.SYMBOL(1,K)).AND.&
                (SYMBOL(1,J).EQ.SYMBOL(2,K)))) THEN
               IF (SYMBOL(4,K)(:1).EQ.'J          ') NPATH=NPATH+1
            ENDIF
         ENDDO

! ------ Branching ratio
         IF (NPATH.EQ.0) THEN
            XJ(J)=0.0D+0
         ELSE
            XJ(J)=XJ(J)/REAL(NPATH)
         ENDIF

! ------ Factor of 2 for same species reactions
         IF (JSP1(J).EQ.JSP2(J)) XJ(J)=XJ(J)/2.0D+0
!        print *,SYMBOL(1,J)
!       print *,SYMBOL(2,J)
!       print *,XJ(J)
       
!       stop
       
! ------ Calculate evaporation fraction
         NEVAP=0
         DO K=1,NKMAX
            IF ((SYMBOL(4,J)(:1).EQ.'J          ').AND.(A(K).NE.0.d0)) THEN
               IF ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
                  (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
                  (SYMBOL(4,J)(2:).EQ.SYMBOL(4,K)).AND.&
                  (SYMBOL(5,J)(2:).EQ.SYMBOL(5,K)).AND.&
                  (SYMBOL(6,J)(2:).EQ.SYMBOL(6,K)).AND.&
                  (SYMBOL(4,K)(:1).NE.'J          ')) NEVAP=NEVAP+1
            ENDIF
            IF ((SYMBOL(4,J)(:1).NE.'J          ').AND.(A(J).NE.0.d0)) THEN
               IF ((SYMBOL(1,J).EQ.SYMBOL(1,K)).AND.&
                  (SYMBOL(2,J).EQ.SYMBOL(2,K)).AND.&
                  (SYMBOL(4,J).EQ.SYMBOL(4,K)(2:)).AND.&
                  (SYMBOL(5,J).EQ.SYMBOL(5,K)(2:)).AND.&
                  (SYMBOL(6,J).EQ.SYMBOL(6,K)(2:)).AND.&
                  (SYMBOL(4,K)(:1).EQ.'J          ')) NEVAP=NEVAP+1
            ENDIF
         ENDDO

         N4=0
         N5=0
         N6=0
         DO I=NSGAS+1,NSMAX
            IF (SYMBOL(4,J)(:1).EQ.'J          ') THEN
               IF (SYMBOL(4,J).EQ.SPEC(I)) N4=I
               IF (SYMBOL(5,J).EQ.SPEC(I)) N5=I
               IF (SYMBOL(6,J).EQ.SPEC(I)) N6=I
            ENDIF
            IF ((SYMBOL(4,J)(:1).NE.'J          ').AND.&
           (SYMBOL(4,J)(:1).NE.'X          ')) THEN
               IF (SYMBOL(4,J).EQ.SPEC(I)(2:)) N4=I
               IF (SYMBOL(5,J).EQ.SPEC(I)(2:)) N5=I
               IF (SYMBOL(6,J).EQ.SPEC(I)(2:)) N6=I
            ENDIF
         ENDDO

         DHFSUM=DHF(JSP1(J))+DHF(JSP2(J))-DHF(N4)
         IF (N5.NE.0) DHFSUM=DHFSUM-DHF(N5)
         IF (N6.NE.0) DHFSUM=DHFSUM-DHF(N6)
! ------ Convert from kcal to J, from J to K
         DHFSUM=DHFSUM*4.184D+03/1.38054D-23
! ------ Convert from #moles-1 to #reactions-1
         DHFSUM=DHFSUM/AVOGADRO

         DHFSUM=DHFSUM+EA(J)

         SUM1=ED(N4)
         IF (N5.NE.0) SUM1=MAX(ED(N4),ED(N5))
         IF (N6.NE.0) SUM1=MAX(ED(N4),ED(N5),ED(N6))

         ATOMS=0
         DO K=1,NEMAX
            ATOMS=ATOMS+IELM(K,N4)
         ENDDO

         SUM2=1.0D+0-(SUM1/DHFSUM)
         IF (ATOMS.EQ.2) SUM2=SUM2**(3*ATOMS-5)
         IF (ATOMS.GT.2) SUM2=SUM2**(3*ATOMS-6)
!         SUM2=SUM2**(3*ATOMS-6)
         SUM2=ARRK*SUM2
         EVFRAC=SUM2/(1+SUM2)
         
!        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
!         IF (SYMBOL(4,J).EQ.'H2O     ') THEN
!                 EVFRAC=0.009
!                EVFRAC_H2O=0.009
!         ENDIF
         
         BADFLAG=0
         IF (DHF(JSP1(J)).LE.-999.0) THEN
            EVFRAC=0.0D+0
            BADFLAG=BADFLAG+1
         ENDIF
         IF (DHF(JSP2(J)).LE.-999.0) THEN
            EVFRAC=0.0D+0
            BADFLAG=BADFLAG+1
         ENDIF
         IF (DHF(N4).LE.-999.0) THEN
            EVFRAC=0.0D+0
            BADFLAG=BADFLAG+1
         ENDIF
         IF (N5.NE.0) THEN
            EVFRAC=0.0D+0
            BADFLAG=BADFLAG+1
         ENDIF
         IF (N6.NE.0) THEN
            EVFRAC=0.0D+0
            BADFLAG=BADFLAG+1
         ENDIF

         IF (EVFRAC.GE.1.0D+0) EVFRAC=1.0D+0
         IF (EVFRAC.LE.0.0D+0) EVFRAC=0.0D+0
         IF (NEVAP.EQ.0) EVFRAC=0.0D+0
         IF (DHFSUM.LE.0.0D+0) EVFRAC=0.0D+0

!         EVFRAC=0.0D+0

         IF (SYMBOL(4,J)(:1).EQ.'J          ') THEN
            EVFRAC=1.0D+0-EVFRAC
         ENDIF
         
!        V.W. Jul 2006 f=evfrac=0.009 for H2O (Kroes & Anderson 2006) 
!        IF (SYMBOL(4,J).EQ.'JH2O    ') THEN
!            EVFRAC=1.0D+0-EVFRAC_H2O
!         ENDIF

         XJ(J)=XJ(J)*EVFRAC
!        print *,SYMBOL(1,J)
!       print *,SYMBOL(2,J)
!       print *,XJ(J)
       
!       stop

! ------ Calculate quantum activation energy
         REDMAS=REAL(SMASS(JSP1(J))*SMASS(JSP2(J)))/&
                REAL(SMASS(JSP1(J))+SMASS(JSP2(J)))
         ACT1(J)=2.0D+0 * ACT/HBAR * SQRT(2.0D+0*AMU*REDMAS*BOLTZ*EA(J))
      ENDIF

! === ITYPE 16 - C.R. DESORPTION
      IF (ITYPE(J).EQ.16) THEN
         IF (SMASS(JSP1(J)).EQ.0) XJ(J)=0.0D+0
      ENDIF

! === ITYPE 99 - ACCRETION
      IF (ITYPE(J).EQ.99) THEN
! ------ Save tag of resultant grain surface species
         DO I=1,NSMAX
            IF (SYMBOL(4,J).EQ.SPEC(I)) JSP2(J)=I
         ENDDO
      ENDIF

      ENDDO

! === Zero dummy H2 formation rxns, if necc.
!      IF (IDUST.NE.0) THEN
!         XJ(1)=0.0D+0
!         XJ(2)=0.0D+0
!      ENDIF

      RETURN
      END

! ======================================================================
! ======================================================================
