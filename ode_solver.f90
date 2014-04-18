module ode_solver

use shielding

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Initialize reactants and products of all reactions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_chemical_reactants()
use global_variables

implicit none

! Locals
integer :: I,J,L

integer :: no_species

no_species = nb_species + 1

! By default, non existing reactants (dummy species) will be assigned (nb_species+1)
reaction_substances(1:7, 1:nb_reactions) = no_species

do I=1,nb_reactions
  do J=1,nb_species

    do L=1,7
      if (SYMBOL(L,I).EQ.species_name(J)) then
        reaction_substances(L,I) = J
      endif
    enddo

  enddo
enddo   

return
end subroutine set_chemical_reactants

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Count the number of non-zeros elements in each line of the jacobian
!! to dimension the arrays used in ODEPACK.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine count_nonzeros()
use global_variables
implicit none

! Locals
integer :: i

! Dummy parameters for restricted call of get_jacobian
real(double_precision), dimension(nb_species) :: DUMMYPDJ, DUMMYY
integer IDUMMY
integer, parameter :: dummy_n = 3
real(double_precision), parameter :: dummy_t = 0.d0
real(double_precision), dimension(dummy_n) :: dummy_ian, dummy_jan

integer :: max_nonzeros, NUMBERJAC

! TODO comment not needed anymore maybe include it in the routine doxygen documentation
! Dimension of the work arrays for the solver 
! The number of non zero values is checked with the testjac flag
! nb_nonzeros_values should be around the largest printed value

! Forced initialisation of global variables that will be needed, especially for the 'set_constant_rates' part. We donc care about specific values, 
!! all that counts is that we can retrieve the number of non-zeros elements.
iptstore = 1

max_nonzeros = 0

dummyy(1:nb_species) = 1.d-5

call set_constant_rates()
call set_dependant_rates(dummyy)

do IDUMMY=1,nb_species
  call get_jacobian(n=dummy_n, t=dummy_t, y=dummyy,j=idummy,ian=dummy_ian, jan=dummy_jan, pdj=dummypdj)
    
  NUMBERJAC=0
  do i=1,nb_species
    if (dummypdj(i).ne.0.d0) then
      NUMBERJAC = NUMBERJAC + 1
    endif
  enddo

  if (NUMBERJAC.gt.max_nonzeros) then
    max_nonzeros = NUMBERJAC
  endif
  
enddo

nb_nonzeros_values = max_nonzeros

return
end subroutine count_nonzeros

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou & Franck Hersant
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Computes a column (for the J-th species) of the chemical jacobian. 
!
!> @warning Even if N, T, IAN, and JAN are not actually used, they are needed
!! because ODEPACK need a routine with a specific format, and specific inputs
!! and outputs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_jacobian(N, T, Y, J, IAN, JAN, PDJ)
use global_variables
implicit none

! Inputs
integer, intent(in) :: N
integer, intent(in) :: J
real(double_precision), intent(in) :: T
real(double_precision), intent(in), dimension(N) :: IAN
real(double_precision), intent(in), dimension(N) :: JAN
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances

! Outputs
real(double_precision), intent(out), dimension(nb_species) :: PDJ

! Locals
integer :: no_species

real(double_precision), dimension(nb_species+1) :: PDJ2
integer :: i
integer :: reactant1_idx, reactant2_idx, reactant3_idx, product1_idx, product2_idx, product3_idx, product4_idx
integer :: reaction_idx ! The index of a given reaction

! Temp values to increase speed
real(double_precision) :: H_number_density_squared ! H_number_density*H_number_density, to gain speed
real(double_precision) :: tmp_value ! To optimize speed, temporary variable is created to avoid multiple calculation of the same thing

no_species=nb_species+1 ! Index corresponding to no species (meaning that there is no 3rd reactant for instance

H_number_density_squared = H_number_density * H_number_density

PDJ2(1:nb_species+1) = 0.d0

do i=1,nb_reactions_using_species(j)
  reaction_idx = relevant_reactions(i, j) ! j being the species index, given as a parameter

  reactant1_idx = reaction_substances(1, reaction_idx)
  reactant2_idx = reaction_substances(2, reaction_idx)
  reactant3_idx = reaction_substances(3, reaction_idx)
  product1_idx = reaction_substances(4, reaction_idx)
  product2_idx = reaction_substances(5, reaction_idx)
  product3_idx = reaction_substances(6, reaction_idx)
  product4_idx = reaction_substances(7, reaction_idx)
  
  ! if statements are written in a specific order to increase speed. The goal is to test first the most probable event, and 
  !! then always go to 'else' statement, not to test if we have already found our case. One, then two bodies reactions are the most 
  !! abundants reactions. 

  ! One reactant only
  if (reactant2_idx.eq.no_species) then
    if (reactant1_idx.eq.J) then 
      tmp_value = XK(reaction_idx)
      PDJ2(product1_idx)   = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx)   = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx)   = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx)   = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
    endif
  
  ! Two bodies reaction
  else if (reactant3_idx.eq.no_species) then
    if (reactant1_idx.eq.J) then 
      tmp_value = XK(reaction_idx) * Y(reactant2_idx) * H_number_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
    endif

    if (reactant2_idx.eq.J) then 
      tmp_value = XK(reaction_idx) * Y(reactant1_idx) * H_number_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
    endif
  
  ! Three bodies reaction
  else
    if (reactant1_idx.eq.J) then 
      tmp_value = XK(reaction_idx) * Y(reactant2_idx) * Y(reactant3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

    if (reactant2_idx.eq.J) then 
      tmp_value = XK(reaction_idx) * Y(reactant1_idx) * Y(reactant3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

    if (reactant3_idx.eq.J) then 
      tmp_value = XK(reaction_idx) * Y(reactant1_idx) * Y(reactant2_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reactant1_idx) = PDJ2(reactant1_idx) - tmp_value
      PDJ2(reactant2_idx) = PDJ2(reactant2_idx) - tmp_value
      PDJ2(reactant3_idx) = PDJ2(reactant3_idx) - tmp_value
    endif

  endif

enddo

PDJ(1:nb_species)=PDJ2(1:nb_species)

return
end subroutine get_jacobian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read the full list of reactions and retrieve, for each species
!! the list of reactions involving it. This will be used in get_jacobian 
!! to increase speed. 
!!\n\n
!! max_reactions_same_species is set here. nb_reactions_using_species and relevant_reactions array are set here.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_relevant_reactions()
use global_variables
implicit none

! Locals
integer, dimension(nb_reactions, nb_species+1) :: is_species_used ! For each species, tell which reactions use it or not 
!! (0 if not used, 1 if used at least once)

integer :: reaction, species, idx
integer :: reactant1_idx, reactant2_idx, reactant3_idx

is_species_used(1:nb_reactions, 1:nb_species+1) = 0

do reaction=1,nb_reactions

  reactant1_idx = reaction_substances(1, reaction)
  reactant2_idx = reaction_substances(2, reaction)
  reactant3_idx = reaction_substances(3, reaction)
  
  is_species_used(reaction, reactant1_idx) = 1
  is_species_used(reaction, reactant2_idx) = 1
  is_species_used(reaction, reactant3_idx) = 1

enddo

! We get the total number of reactions in which each species can be involved
! We skip the 'nb_species+1' species that is only a fake species for "no reactant"
nb_reactions_using_species(1:nb_species) = sum(is_species_used(1:nb_reactions, 1:nb_species), 1)

! What is the maximum number of reactions involving one particular species?
max_reactions_same_species = maxval(nb_reactions_using_species(1:nb_species))

allocate(relevant_reactions(max_reactions_same_species, nb_species))

relevant_reactions(1:max_reactions_same_species, 1:nb_species) = 0 ! For the extra elements (because not all species 
!! will have 'max_reactions' reactions involving it).

! For each species, we get the references of reactions that have it as a reactant. The number of reactions is different for each species 
!! Thus, at least one species will have a full line of meaningfull indexes. The other will have the rest of their line completed by zeros.
do species=1,nb_species
  idx = 1
  do reaction=1, nb_reactions
    if (is_species_used(reaction, species).eq.1) then
      relevant_reactions(idx, species) = reaction
      idx = idx + 1
    endif
  enddo
enddo


return
end subroutine init_relevant_reactions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Initialize working arrays used for ODEPACK. We must know 
!! in advance the maximum number of non-zeros elements. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialize_work_arrays()
use global_variables
implicit none

lrw = 20 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species
liw = 31 + 3 * nb_nonzeros_values*nb_species + 21 * nb_species

allocate(iwork(liw))
allocate(rwork(lrw))


return
end subroutine initialize_work_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2008
!
! DESCRIPTION: 
!> @brief Calculates the position of non-zero values in the jacobian. 
!! Set the global variables iwork and rwork
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_work_arrays(Y)
use global_variables
implicit none

! Inputs
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances

! Locals
integer :: i,j,k
real(double_precision), dimension(nb_species) :: PDJ
integer :: NNZ

! Dummy parameters for restricted call of get_jacobian
integer, parameter :: dummy_n = 3
real(double_precision), parameter :: dummy_t = 0.d0
real(double_precision), dimension(dummy_n) :: dummy_ian, dummy_jan


! For IA and JA
integer, dimension(nb_species+1) :: IA !< For each species, the number (minus 1) of non-zeros values for the previous species
integer, dimension(liw) :: JA !< List the non-zeros values. For each non-zeros values, tell to what species it correspond 

call set_constant_rates()
call set_dependant_rates(Y)

! Initialize work arrays

iwork(1:liw) = 0
rwork(1:lrw) = 0.d0
IWORK(5) = 5
RWORK(6) = 3.154D14
IWORK(6) = 10000
IWORK(7) = 2

if (timestep.eq.1) then
  IWORK(6)=2000
endif

k=1

do j=1,nb_species
  call get_jacobian(n=dummy_n, t=dummy_t, y=Y,j=J,ian=dummy_ian, jan=dummy_jan, pdj=PDJ)

  IA(j)=k

  do i=1,nb_species
    if (abs(PDJ(i)).gt.1.d-99) then
      JA(k)=i
      k=k+1
    endif
  enddo

enddo

IA(nb_species+1)=k

NNZ=IA(nb_species+1)-1
iwork(30+1:30+nb_species+1)=IA(1:nb_species+1)
iwork(31+nb_species+1:31+nb_species+NNZ)=JA(1:NNZ)

return
end subroutine set_work_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2008
!
! DESCRIPTION: 
!> @brief Computes the chemical evolution. In particular, calculate all the
!! derivatives for abundances.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_temporal_derivatives(N,T,Y,YDOT)
use global_variables

implicit none

! Inputs
integer, intent(in) :: N
real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances
real(double_precision), intent(in) :: T ! Not used by the code, but needed for the ODEPACK call format expected for FCHEM in dlsodes

! Outputs
real(double_precision), intent(out), dimension(nb_species) :: YDOT

! Locals
integer :: no_species
real(double_precision), dimension(nb_species+1) :: YD2
!REAL(KIND=16), dimension(nb_species+1) :: YD2
integer :: i
integer :: reactant1_idx, reactant2_idx, reactant3_idx, product1_idx, product2_idx, product3_idx, product4_idx
real(double_precision) :: rate

call set_dependant_rates(y)


no_species = nb_species + 1

ydot(1:nb_species) = 0.d0
yd2(1:nb_species) = 0.d0

! The differential equations are calculated in a loop here
do I=1,nb_reactions

  reactant1_idx = reaction_substances(1, i)
  reactant2_idx = reaction_substances(2, i)
  reactant3_idx = reaction_substances(3, i)

  product1_idx = reaction_substances(4, i)
  product2_idx = reaction_substances(5, i)
  product3_idx = reaction_substances(6, i)
  product4_idx = reaction_substances(7, i)

  ! One reactant only
  if (reactant2_idx.eq.no_species) then
    RATE = XK(I) * Y(reactant1_idx)  
  else
    if (reactant3_idx.eq.no_species) then
      ! Two bodies reactions
      RATE = XK(I) * Y(reactant1_idx) * Y(reactant2_idx) * H_number_density
    else 
      ! Three bodies reactions
      RATE = XK(I)*Y(reactant1_idx) * Y(reactant2_idx) * Y(reactant3_idx) * H_number_density * H_number_density
    endif
  endif

  YD2(product1_idx) = YD2(product1_idx) + RATE
  YD2(product2_idx) = YD2(product2_idx) + RATE
  YD2(product3_idx) = YD2(product3_idx) + RATE
  YD2(product4_idx) = YD2(product4_idx) + RATE
  YD2(reactant1_idx) = YD2(reactant1_idx) - RATE
  YD2(reactant2_idx) = YD2(reactant2_idx) - RATE
  YD2(reactant3_idx) = YD2(reactant3_idx) - RATE
enddo   

YDOT(1:nb_species) = YD2(1:nb_species)

return
end subroutine get_temporal_derivatives

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Valentine Wakelam
!
!> @date 2012
!
! DESCRIPTION: 
!> @brief Reactions coefficient formally dependent on the abundances Y are 
!! computed in a companion subroutine: set_dependant_rates
!! Grain surface reactions, self shielding, etc...
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine set_constant_rates()

  use global_variables 
  
  implicit none

  ! Locals
  real(double_precision) :: T300, TI, TSQ
  integer :: nsta, nfin
  integer :: k, j, w, m, n
  integer, dimension(10) :: indice
  real(double_precision), dimension(10) :: distmin, distmax

  T300=gas_temperature/300.D+0
  TI=1.0d00/gas_temperature
  TSQ=SQRT(gas_temperature)

  ! ====== Rxn ITYPE 0
  ! ITYPE 0: Gas phase reactions with GRAINS =) 
  do J=IRXSTA(0),IRXFIN(0)
    XK(J)=A(J)*(T300**B(J))
  enddo

  ! In the case IS_GRAIN_REACTIONS eq 0, we still need the formation of H on the grains
  ! this is done with the XH species, reaction types 10 and 11

  ! ====== Rxn ITYPE 10 and 11
  ! ITYPE 10 and 11: H2 formation on the grains when IS_GRAIN_REACTIONS eq 0
  if (IS_GRAIN_REACTIONS.eq.0) then
    do J=IRXSTA(10),IRXFIN(10)
      XK(J)=A(J)*1.186D7*exp(225.D0/gas_temperature)**(-1)*GTODN/H_number_density
    enddo
    do J=IRXSTA(11),IRXFIN(11)       
      XK(J)=A(J)*(T300**B(J))*H_number_density/GTODN
    enddo       
  endif


  ! ====== Rxn ITYPE 1
  ! ITYPE 1: Photodissoc/ionisation with cosmic rays
  ! Add X-rays in case X_IONISATION_RATE is not 0
  do J=IRXSTA(1),IRXFIN(1)
    XK(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
  enddo
  do J=IRXSTA(2),IRXFIN(2)
    XK(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
  enddo

  ! ====== Rxns ITYPE 4 - 8
  ! Bimolecular gas phase reactions - several possible formula 
  W=1
  NSTA=0
  NFIN=0
  distmin(:)=9999.d0
  distmax(:)=9999.d0
  do J=4,8
    if ((IRXSTA(J).NE.0).AND.(NSTA.EQ.0)) NSTA=J
    if (IRXFIN(J).NE.0) NFIN=J
  enddo
  do J=IRXSTA(NSTA),IRXFIN(NFIN)

    !---------------  KOOIJ FORMULA
    if (FORMULA(J).eq.3) then

      XK(J)=A(J)*(T300**B(J))*EXP(-C(J)/gas_temperature)

      !              Check for temperature bounderies
      if (gas_temperature.LT.Tmin(J)) XK(J)=A(J)*((Tmin(J)/300.D0)**B(J))*EXP(-C(J)/Tmin(J))
      if (gas_temperature.GT.Tmax(J)) XK(J)=A(J)*((Tmax(J)/300.D0)**B(J))*EXP(-C(J)/Tmax(J))

      !              Check for the presence of several rate coefficients present in the network for the
      !              the same reaction
      if (NUM(J+1).EQ.NUM(J)) then
        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)
        W = W + 1
      endif

      if ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)

        do M=1,W
          N=INDICE(M)
          !IF(IT==1) write(*,*) N,M, SYMBOL(:,N), tmin(N), tmax(N),distmin(M),distmax(M)
          if (gas_temperature.LT.Tmin(N)) XK(N)=0.d0
          if (gas_temperature.GT.Tmax(N)) XK(N)=0.d0
        enddo

        if (maxval(XK(indice(1:w))).lt.1.d-99) then

          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*((Tmin(N)/300.D0)**B(N))*EXP(-C(N)/Tmin(N))
          else
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*((Tmax(N)/300.D0)**B(N))*EXP(-C(N)/Tmax(N))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  IONPOL1 FORMULA
    if (FORMULA(J).EQ.4) then
      XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/gas_temperature)**0.5))

      !              Check for temperature bounderies
      if (gas_temperature.LT.Tmin(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/Tmin(J))**0.5))
      if (gas_temperature.GT.Tmax(J)) XK(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TMAX(J))**0.5))

      !              Check for the presence of several rate coefficients present in the network for the
      !              the same reaction
      if (NUM(J+1).EQ.NUM(J)) then
        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)
        W = W + 1
      endif

      if ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)

        do M=1,W
          N=INDICE(M)
          if (gas_temperature.LT.Tmin(N)) XK(N)= 0.d0
          if (gas_temperature.GT.Tmax(N)) XK(N)= 0.d0
        enddo

        if (maxval(XK(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/Tmin(N))**0.5))
          else
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/TMAX(N))**0.5))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  IONPOL2 FORMULA
    if (FORMULA(J).EQ.5) then
      XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/gas_temperature)**0.5)+(C(J)**2*300.d0/(10.526*gas_temperature)))

      !               Check for temperature bounderies
      if (gas_temperature.LT.Tmin(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMIN(J))**0.5)+(C(J)**2*300.d0/(10.526*TMIN(J))))
      if (gas_temperature.GT.Tmax(J)) XK(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMAX(J))**0.5)+(C(J)**2*300.d0/(10.526*TMAX(J))))

      !               Check for the presence of several rate coefficients present in the network for the
      !               the same reaction
      if (NUM(J+1).EQ.NUM(J)) then
        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)
        W = W + 1
      endif

      if ((NUM(J+1).NE.NUM(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=tmin(j)-gas_temperature
        distmax(w)=gas_temperature-tmax(j)

        do M=1,W
          N=INDICE(M)
          if (gas_temperature.LT.Tmin(N)) XK(N)=       0.d0
          if (gas_temperature.GT.Tmax(N)) XK(N)=       0.d0
        enddo

        if (maxval(XK(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMIN(N))**0.5)+(C(N)**2*300.d0/(10.526*TMIN(N))))
          else
            N=indice(minloc(abs(distmax),dim=1))
            XK(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMAX(N))**0.5)+(C(N)**2*300.d0/(10.526*TMAX(N))))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif
  enddo

  ! === Grain surfaces
  if (IS_GRAIN_REACTIONS.NE.0) then

    ! ========= Set diffusion and evaporation rates (s-1)
    do K=1,nb_species
      TINDIF(K)=CHF(K)*EXP(-EB(K)/dust_temperature)/nb_sites_per_grain
      TINEVA(K)=CHF(K)*EXP(-ED(K)/dust_temperature)
    enddo

    ! ========= Rxn ITYPE 15 - thermal evaporation
    ! ITYPE 15: Thermal evaporation
    do J=IRXSTA(15),IRXFIN(15)
      XK(J)=A(J)*XJ(J)*TINEVA(JSP1(J))
    enddo

    ! ========= Rxn ITYPE 16
    ! ITYPE 16: Cosmic-ray evaporation
    do J=IRXSTA(16),IRXFIN(16)
      XK(J)=A(J)*XJ(J)*((CR_IONISATION_RATE+X_IONISATION_RATE)/1.3D-17)&
      *CHF(JSP1(J))*CRFE*CRT*EXP(-ED(JSP1(J))/TSMAX)
    enddo


    ! Photodesorption, when used appears through ITYPES 66 and 67
    if (irxsta(66).ne.0) then
      ! ========= Rxn ITYPE 66
      ! ITYPE 66: CO photodesorption by external UV
      ! 1.d8 is I_ISRF-FUV from Oberg et al. 2007, ApJ, 662, 23
      do j=irxsta(66),irxfin(66)
        XK(J)=A(J)/SITE_DENSITY*UV_FLUX*1.d8*exp(-2.*visual_extinction) 
      enddo

      ! ========= Rxn ITYPE 67
      ! ITYPE 67: CO photodesorption by CR generated UV
      do j=irxsta(67),irxfin(67)
        XK(J)=A(J)/SITE_DENSITY*1.d4
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
    ! Add X-rays in case X_IONISATION_RATE is not 0
    do J=IRXSTA(17),IRXFIN(17)
      XK(J)=A(J)*(CR_IONISATION_RATE + X_IONISATION_RATE)
      !            if (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    enddo

    ! ====== Rxn ITYPE 18
    ! ITYPE 18: Photodissociations by Cosmic rays on grain surfaces
    do J=IRXSTA(18),IRXFIN(18)
      XK(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
      !            if (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    enddo

  endif

  ! When dust is turned off, zero all dust rates==========================
  if ((IS_GRAIN_REACTIONS.EQ.0).AND.(timestep.EQ.1)) then
    do J=IRXSTA(14),IRXFIN(99)
      XK(J)=0.d0
      XJ(J)=0.d0
    enddo
  endif


  return
  end subroutine set_constant_rates
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Valentine Wakelam
!
!> @date 2012
!
! DESCRIPTION: 
!> @brief Set Reactions coefficient formally dependent on the abundances Y. 
!! Grain surface reactions, self shielding, etc...
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine set_dependant_rates(Y)

  use global_variables 
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in), dimension(nb_species) :: Y !< [in] abundances

  ! Locals
  real(double_precision) :: ACTIV,BARR,MONLAY,DIFF
  real(double_precision) :: XNH2,XNCO
  real(double_precision) :: TETABIS,TETABIS1,TETABIS2,TETABIS3
  real(double_precision) :: T300, TI, TSQ
  real(double_precision) :: YMOD1, YMOD2
  integer IMOD1,IMOD2
  integer :: j, l

  T300=gas_temperature/300.D+0
  TI=1.0d00/gas_temperature
  TSQ=SQRT(gas_temperature)
  MONLAY=LAYERS*nb_sites_per_grain/GTODN

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
  do J=IRXSTA(3),IRXFIN(3)
    XK(J)=A(J)*EXP(-C(J)*visual_extinction)*UV_FLUX

    ! MODIFY THE H2 AND CO PHOTODISSOCIATION if IS_ABSORPTION EQ 1
    if (IS_ABSORPTION.EQ.1) then 

      ! ====== Compute the H2 self-shielding
      if (SYMBOL(1,J).EQ.YH2) then
        TETABIS=1.D0

        NH2 = visual_extinction/AV_NH_ratio * XNH2

        ! ======= Linear extrapolation of the shielding factors
        do L=1,NL1-1
          if ((N1H2(L).LE.NH2).AND.(N1H2(L+1).GE.NH2)) then
            TETABIS=T1H2(L)+(NH2-N1H2(L))*(T1H2(L+1)-T1H2(L))/(N1H2(L+1)-N1H2(L))
          endif
        enddo
        if (NH2.GT.N1H2(NL1)) TETABIS = T1H2(NL1)

        XK(J)=2.54D-11*TETABIS

        XK(J)=XK(J)*UV_FLUX
      endif

      ! ====== Compute the CO self-shielding
      if (SYMBOL(1,J).EQ.YCO) then

        TETABIS1=1.D0
        TETABIS2=1.D0
        TETABIS3=1.D0

        NH2 = visual_extinction/AV_NH_ratio * XNH2
        NCO = visual_extinction/AV_NH_ratio * XNCO

        ! ======= Linear extrapolation of the three shileding factors
        do L=1,NL2-1
          if ((N2CO(L).LE.NCO).AND.(N2CO(L+1).GE.NCO))  &
          TETABIS2=T2CO(L)+(NCO-N2CO(L))*(T2CO(L+1)-T2CO(L))&
          /(N2CO(L+1)-N2CO(L))
        enddo

        do L=1,NL3-1
          if ((N2H2(L).LE.NH2).AND.(N2H2(L+1).GE.NH2)) &
          TETABIS1=T2H2(L)+(NH2-N2H2(L))*(T2H2(L+1)-T2H2(L))&
          /(N2H2(L+1)-N2H2(L))
          if ((AV2(L).LE.visual_extinction).AND.(AV2(L+1).GE.visual_extinction)) &
          TETABIS3=T2AV(L)+(visual_extinction-AV2(L))*(T2AV(L+1)-T2AV(L))&
          /(AV2(L+1)-AV2(L))
        enddo

        ! Saturate the rate coefficient if necessary (when density or Av are out of 
        ! the shielding array, from the photodiss files)

        if (NCO.GT.N2CO(NL2)) TETABIS2 = T2CO(NL2)
        if (NH2.GT.N2H2(NL3)) TETABIS1 = T2H2(NL3)
        if (visual_extinction.GT.AV2(NL3))  TETABIS3 = T2AV(NL3)

        XK(J)=1.03D-10*TETABIS1*TETABIS2*TETABIS3
        XK(J)=XK(J)*UV_FLUX

      endif

    endif

  enddo

  ! Continually time-dependent grain rates================================
  if (IS_GRAIN_REACTIONS.NE.0) then

    ! ====== Rxn ITYPE 99
    ! ITYPE 99: Adsorption on grains
    do J=IRXSTA(99),IRXFIN(99)
      ! ========= Set accretion rates
      TINACC(JSP1(J))=CONDSP(JSP1(J))*TSQ*Y(JSP1(J))*H_number_density
      TINACC(JSP2(J))=TINACC(JSP1(J))
      XK(J)=A(J)*XJ(J)*TINACC(JSP1(J))/Y(JSP1(J))/GTODN
    enddo

    ! ====== Rxn ITYPE 14
    ! ITYPE 14: Grain surface reactions
    do J=IRXSTA(14),IRXFIN(14)
      IMOD1=0
      IMOD2=0
      BARR=1.0d0
      ! --------- Calculate activation energy barrier multiplier
      if (EA(J).GE.1.0D-40) then
        ACTIV=EA(J)/dust_temperature
        ! ------------ Choose fastest of classical or tunnelling
        if (ACTIV.GT.ACT1(J)) ACTIV=ACT1(J)
        BARR=EXP(-ACTIV)
      endif

      ! --------- Thermal hopping diffusion method
      RDIF1(J)=TINDIF(JSP1(J))
      RDIF2(J)=TINDIF(JSP2(J))

      ! --------- Check for JH,JH2
      if (SYMBOL(1,J).EQ.YJH)  IMOD1=1
      if (SYMBOL(1,J).EQ.YJH2) IMOD1=2
      if (SYMBOL(2,J).EQ.YJH)  IMOD2=1
      if (SYMBOL(2,J).EQ.YJH2) IMOD2=2

      ! --------- QM for JH,JH2 only - others are too heavy
      if (IMOD1+IMOD2.NE.0) then
        ! ------------ QM1 - Tunnelling (if it's faster than thermal)
        if (GRAIN_TUNNELING_DIFFUSION.EQ.1) then
          if ((IMOD1.NE.0).AND.&
          (RQ1(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ1(JSP1(J))
          if ((IMOD2.NE.0).AND.&
          (RQ1(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ1(JSP2(J))
        endif
        ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
        if (GRAIN_TUNNELING_DIFFUSION.EQ.2) then
          if ((IMOD1.NE.0).AND.&
          (RQ2(JSP1(J)).GT.RDIF1(J))) RDIF1(J)=RQ2(JSP1(J))
          if ((IMOD2.NE.0).AND.&
          (RQ2(JSP2(J)).GT.RDIF2(J))) RDIF2(J)=RQ2(JSP2(J))
        endif
        ! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
        if (GRAIN_TUNNELING_DIFFUSION.EQ.3) then
          if (IMOD1.NE.0) then
            if (RQ1(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ1(JSP1(J))
            if (RQ2(JSP1(J)).GT.RDIF1(J)) RDIF1(J)=RQ2(JSP1(J))
          endif
          if (IMOD2.NE.0) then
            if (RQ1(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ1(JSP2(J))
            if (RQ2(JSP2(J)).GT.RDIF2(J)) RDIF2(J)=RQ2(JSP2(J))
          endif
        endif
      endif

      ! --------- Modify according to IMODH switch:
      if (IMODH.NE.0) then
        ! ------------ if H+H->H2 is only modified rxn:
        if ((IMODH.EQ.-1).AND.(IMOD1.NE.1.OR.IMOD2.NE.1)) then
          IMOD1=0
          IMOD2=0
        endif

        ! ------------ if only H is modified:
        if ((IMODH.EQ.1).AND.(IMOD1.NE.1)) IMOD1=0
        if ((IMODH.EQ.1).AND.(IMOD2.NE.1)) IMOD2=0

        ! ------------ Set to modify all rates, if selected (just atoms)
        if (IMODH.EQ.3) then
          if ((SYMBOL(1,J).EQ.YJH).OR.&
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
          if ((SYMBOL(2,J).EQ.YJH).OR.&
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
        endif

        ! ------------ Modify rates (RDIF1 & RDIF2) according to their own evap/acc rates
        YMOD1=Y(JSP1(J))
        YMOD2=Y(JSP2(J))
        call modify_specific_rates(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
      endif

      DIFF=RDIF1(J)+RDIF2(J)

      XK(J)=A(J)*XJ(J)*BARR*DIFF*GTODN/H_number_density
      !              XK(J)=0.D0
      ! --------- Allow only 1 monolayer of each to react
      ! Not used for the time being
      !            if (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*nb_sites_per_grain/GTODN/Y(JSP1(J))
      !            if (Y(JSP2(J)).GT.MONLAY) XK(J)=XK(J)*nb_sites_per_grain/GTODN/Y(JSP2(J))

    enddo

    ! ====== Rxn ITYPE 19 - 20
    ! ITYPE 19: Photodissociations by UV photons on grain surfaces
    ! ITYPE 20: Photodissociations by UV photons on grain surfaces
    do J=IRXSTA(19),IRXFIN(20)
      XK(J)=A(J)*EXP(-C(J)*visual_extinction)*UV_FLUX
      !            if (Y(JSP1(J)).GT.MONLAY) XK(J)=XK(J)*MONLAY/Y(JSP1(J))
    enddo

    ! Useful for testing
    ! To disable some reaction types
    !do j=irxsta(14),irxfin(15)
    !xk(j)=0.
    !enddo

    !do j=irxsta(17),irxfin(20)
    !xk(j)=0.
    !enddo
  endif

  ! Continually time-dependent gas phase rates============================
  ! H2 formation
  ! XJ(1) and XJ(2) are zero if IS_GRAIN_REACTIONS=1
  ! cf GRAINRATE
  ! VW Fev 2012 - this process has been removed
  !      do j=irxsta(0),irxfin(0)
  !      if ((SYMBOL(1,J).eq.YH).and.(SYMBOL(2,j).eq.YH)) then
  !      XK(j)=XJ(j)*A(j)*(T300**B(j))*GTODN/H_number_density/Y(JSP1(j))
  !      endif
  !      enddo

  ! if rate acoefficients are too small, put them to 0
  where (XK.lt.RXNMIN) XK=0.d0

    return
    end subroutine set_dependant_rates
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Modify some reaction rates on the grain surface (itype=14) in
!! various conditions. Test to estimates the fastest process and replace them.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine modify_specific_rates(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
  use global_variables
  
  implicit none

  ! Inputs
  integer, intent(in) :: J
  integer, intent(in) :: IMOD1
  integer, intent(in) :: IMOD2
  real(double_precision), intent(in) :: BARR
  real(double_precision), intent(in) :: YMOD1
  real(double_precision), intent(in) :: YMOD2
  
  ! Locals
  real(double_precision) :: PICK, TESTREF1, TESTREF2, TESTNUM

  EX1(J)=0.d0
  EX2(J)=0.d0

  ! --- Check value of x = t_acc/t_evap
  ! TINEVA = 1/t_evap
  ! TINACC = 1/t_acc
  if (TINACC(JSP1(J)).GT.0.d0) then
    EX1(J)=TINEVA(JSP1(J))/TINACC(JSP1(J))
  endif
  if (TINACC(JSP2(J)).GT.0.d0) then
    EX2(J)=TINEVA(JSP2(J))/TINACC(JSP2(J))
  endif
  ! Hence x = 0 if t_evap or t_acc = 0

  ! --- Assign max rates

  if (BARR.EQ.1.0d0) then
    if (IMOD1.NE.0) then
      if (EX1(J).LT.1.0d0) then
        if (RDIF1(J).GT.TINACC(JSP1(J))) RDIF1(J)=TINACC(JSP1(J))
      endif
      if (EX1(J).GE.1.0d0) then
        if (RDIF1(J).GT.TINEVA(JSP1(J))) RDIF1(J)=TINEVA(JSP1(J))
      endif
    endif

    if (IMOD2.NE.0) then
      if (EX2(J).LT.1.0d0) then
        if (RDIF2(J).GT.TINACC(JSP2(J))) RDIF2(J)=TINACC(JSP2(J))
      endif
      if (EX2(J).GE.1.0d0) then
        if (RDIF2(J).GT.TINEVA(JSP2(J))) RDIF2(J)=TINEVA(JSP2(J))
      endif
    endif
  endif

  ! --- Species rate to compare chosen by fastest diffusion rate
  if (BARR.NE.1.0d0) then
    PICK=0.d0

    TESTREF1=TINACC(JSP1(J))
    if (EX1(J).GE.1.0d0) TESTREF1=TINEVA(JSP1(J))
    TESTREF2=TINACC(JSP2(J))
    if (EX2(J).GE.1.0d0) TESTREF2=TINEVA(JSP2(J))

    if (RDIF1(J).GE.RDIF2(J)) then
      TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD2*GTODN
      if (YMOD2*GTODN.LT.1.0d0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
      if (TESTNUM.GT.TESTREF1) PICK=1.d0
    endif
    if (RDIF2(J).GT.RDIF1(J)) then
      TESTNUM=(RDIF1(J)+RDIF2(J))*BARR*YMOD1*GTODN
      if (YMOD1*GTODN.LT.1.0d0) TESTNUM=(RDIF1(J)+RDIF2(J))*BARR
      if (TESTNUM.GT.TESTREF2) PICK=2.d0
    endif

    if (PICK.EQ.1) then
      RDIF1(J)=TESTREF1/BARR/YMOD2/GTODN
      if (YMOD2*GTODN.LT.1.0d0) RDIF1(J)=TESTREF1/BARR
      RDIF2(J)=0.d0
    endif

    if (PICK.EQ.2) then
      RDIF2(J)=TESTREF2/BARR/YMOD1/GTODN
      if (YMOD1*GTODN.LT.1.0d0) RDIF2(J)=TESTREF2/BARR
      RDIF1(J)=0.d0
    endif

  endif

  return
  end subroutine modify_specific_rates

end module ode_solver