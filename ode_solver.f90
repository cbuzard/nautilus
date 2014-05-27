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
!> @brief Initialize reagents and products of all reactions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine set_chemical_reagents()
use global_variables

implicit none

! Locals
integer :: I,J,L

integer :: no_species

no_species = nb_species + 1

! By default, non existing reagents (dummy species) will be assigned (nb_species+1)
REACTION_SUBSTANCES_ID(1:7, 1:nb_reactions) = no_species

do I=1,nb_reactions
  do J=1,nb_species

    do L=1,7
      if (REACTION_SUBSTANCES_NAMES(L,I).EQ.species_name(J)) then
        REACTION_SUBSTANCES_ID(L,I) = J
      endif
    enddo

  enddo
enddo   

return
end subroutine set_chemical_reagents

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
integer :: reagent1_idx, reagent2_idx, reagent3_idx, product1_idx, product2_idx, product3_idx, product4_idx
integer :: reaction_idx ! The index of a given reaction

! Temp values to increase speed
real(double_precision) :: H_number_density_squared ! H_number_density*H_number_density, to gain speed
real(double_precision) :: tmp_value ! To optimize speed, temporary variable is created to avoid multiple calculation of the same thing

no_species=nb_species+1 ! Index corresponding to no species (meaning that there is no 3rd reagent for instance

H_number_density_squared = H_number_density * H_number_density

PDJ2(1:nb_species+1) = 0.d0

do i=1,nb_reactions_using_species(j)
  reaction_idx = relevant_reactions(i, j) ! j being the species index, given as a parameter

  reagent1_idx = REACTION_SUBSTANCES_ID(1, reaction_idx)
  reagent2_idx = REACTION_SUBSTANCES_ID(2, reaction_idx)
  reagent3_idx = REACTION_SUBSTANCES_ID(3, reaction_idx)
  product1_idx = REACTION_SUBSTANCES_ID(4, reaction_idx)
  product2_idx = REACTION_SUBSTANCES_ID(5, reaction_idx)
  product3_idx = REACTION_SUBSTANCES_ID(6, reaction_idx)
  product4_idx = REACTION_SUBSTANCES_ID(7, reaction_idx)
  
  ! if statements are written in a specific order to increase speed. The goal is to test first the most probable event, and 
  !! then always go to 'else' statement, not to test if we have already found our case. One, then two bodies reactions are the most 
  !! abundants reactions. 

  ! One reagent only
  if (reagent2_idx.eq.no_species) then
    if (reagent1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx)
      PDJ2(product1_idx)   = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx)   = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx)   = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx)   = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
    endif
  
  ! Two bodies reaction
  else if (reagent3_idx.eq.no_species) then
    if (reagent1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reagent2_idx) * H_number_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
      PDJ2(reagent2_idx) = PDJ2(reagent2_idx) - tmp_value
    endif

    if (reagent2_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reagent1_idx) * H_number_density
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
      PDJ2(reagent2_idx) = PDJ2(reagent2_idx) - tmp_value
    endif
  
  ! Three bodies reaction
  else
    if (reagent1_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reagent2_idx) * Y(reagent3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
      PDJ2(reagent2_idx) = PDJ2(reagent2_idx) - tmp_value
      PDJ2(reagent3_idx) = PDJ2(reagent3_idx) - tmp_value
    endif

    if (reagent2_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reagent1_idx) * Y(reagent3_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
      PDJ2(reagent2_idx) = PDJ2(reagent2_idx) - tmp_value
      PDJ2(reagent3_idx) = PDJ2(reagent3_idx) - tmp_value
    endif

    if (reagent3_idx.eq.J) then 
      tmp_value = reaction_rates(reaction_idx) * Y(reagent1_idx) * Y(reagent2_idx) * H_number_density_squared
      PDJ2(product1_idx) = PDJ2(product1_idx) + tmp_value
      PDJ2(product2_idx) = PDJ2(product2_idx) + tmp_value
      PDJ2(product3_idx) = PDJ2(product3_idx) + tmp_value
      PDJ2(product4_idx) = PDJ2(product4_idx) + tmp_value
      PDJ2(reagent1_idx) = PDJ2(reagent1_idx) - tmp_value
      PDJ2(reagent2_idx) = PDJ2(reagent2_idx) - tmp_value
      PDJ2(reagent3_idx) = PDJ2(reagent3_idx) - tmp_value
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
integer :: reagent1_idx, reagent2_idx, reagent3_idx

is_species_used(1:nb_reactions, 1:nb_species+1) = 0

do reaction=1,nb_reactions

  reagent1_idx = REACTION_SUBSTANCES_ID(1, reaction)
  reagent2_idx = REACTION_SUBSTANCES_ID(2, reaction)
  reagent3_idx = REACTION_SUBSTANCES_ID(3, reaction)
  
  is_species_used(reaction, reagent1_idx) = 1
  is_species_used(reaction, reagent2_idx) = 1
  is_species_used(reaction, reagent3_idx) = 1

enddo

! We get the total number of reactions in which each species can be involved
! We skip the 'nb_species+1' species that is only a fake species for "no reagent"
nb_reactions_using_species(1:nb_species) = sum(is_species_used(1:nb_reactions, 1:nb_species), 1)

! What is the maximum number of reactions involving one particular species?
max_reactions_same_species = maxval(nb_reactions_using_species(1:nb_species))

allocate(relevant_reactions(max_reactions_same_species, nb_species))

relevant_reactions(1:max_reactions_same_species, 1:nb_species) = 0 ! For the extra elements (because not all species 
!! will have 'max_reactions' reactions involving it).

! For each species, we get the references of reactions that have it as a reagent. The number of reactions is different for each species 
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

iwork(1:liw) = 0
rwork(1:lrw) = 0.d0

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
integer :: reagent1_idx, reagent2_idx, reagent3_idx, product1_idx, product2_idx, product3_idx, product4_idx
real(double_precision) :: rate

call set_dependant_rates(y)


no_species = nb_species + 1

ydot(1:nb_species) = 0.d0
yd2(1:nb_species) = 0.d0

! The differential equations are calculated in a loop here
do I=1,nb_reactions

  reagent1_idx = REACTION_SUBSTANCES_ID(1, i)
  reagent2_idx = REACTION_SUBSTANCES_ID(2, i)
  reagent3_idx = REACTION_SUBSTANCES_ID(3, i)

  product1_idx = REACTION_SUBSTANCES_ID(4, i)
  product2_idx = REACTION_SUBSTANCES_ID(5, i)
  product3_idx = REACTION_SUBSTANCES_ID(6, i)
  product4_idx = REACTION_SUBSTANCES_ID(7, i)

  ! One reagent only
  if (reagent2_idx.eq.no_species) then
    RATE = reaction_rates(I) * Y(reagent1_idx)  
  else
    if (reagent3_idx.eq.no_species) then
      ! Two bodies reactions
      RATE = reaction_rates(I) * Y(reagent1_idx) * Y(reagent2_idx) * H_number_density
    else 
      ! Three bodies reactions
      RATE = reaction_rates(I)*Y(reagent1_idx) * Y(reagent2_idx) * Y(reagent3_idx) * H_number_density * H_number_density
    endif
  endif

  YD2(product1_idx) = YD2(product1_idx) + RATE
  YD2(product2_idx) = YD2(product2_idx) + RATE
  YD2(product3_idx) = YD2(product3_idx) + RATE
  YD2(product4_idx) = YD2(product4_idx) + RATE
  YD2(reagent1_idx) = YD2(reagent1_idx) - RATE
  YD2(reagent2_idx) = YD2(reagent2_idx) - RATE
  YD2(reagent3_idx) = YD2(reagent3_idx) - RATE
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
!> @brief Constant reactions rates with respect to abundances. 
!! Reactions coefficient formally dependent on the abundances Y are 
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
  do J=type_id_start(0),type_id_stop(0)
    reaction_rates(J)=A(J)*(T300**B(J))
  enddo

  ! In the case IS_GRAIN_REACTIONS eq 0, we still need the formation of H on the grains
  ! this is done with the XH species, reaction types 10 and 11

  ! ====== Rxn ITYPE 10 and 11
  ! ITYPE 10 and 11: H2 formation on the grains when IS_GRAIN_REACTIONS eq 0
  if (IS_GRAIN_REACTIONS.eq.0) then
    do J=type_id_start(10),type_id_stop(10)
      reaction_rates(J)=A(J)*1.186D7*exp(225.D0/gas_temperature)**(-1)*GTODN/H_number_density
    enddo
    do J=type_id_start(11),type_id_stop(11)       
      reaction_rates(J)=A(J)*(T300**B(J))*H_number_density/GTODN
    enddo       
  endif


  ! ====== Rxn ITYPE 1
  ! ITYPE 1: Photodissoc/ionisation with cosmic rays
  ! Add X-rays in case X_IONISATION_RATE is not 0
  do J=type_id_start(1),type_id_stop(1)
    reaction_rates(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
  enddo
  do J=type_id_start(2),type_id_stop(2)
    reaction_rates(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
  enddo

  ! ====== Rxns ITYPE 4 - 8
  ! Bimolecular gas phase reactions - several possible formula 
  W=1
  NSTA=0
  NFIN=0
  distmin(:)=9999.d0
  distmax(:)=9999.d0
  do J=4,8
    if ((type_id_start(J).NE.0).AND.(NSTA.EQ.0)) NSTA=J
    if (type_id_stop(J).NE.0) NFIN=J
  enddo
  do J=type_id_start(NSTA),type_id_stop(NFIN)

    !---------------  KOOIJ FORMULA
    if (RATE_FORMULA(J).eq.3) then

      reaction_rates(J)=A(J)*(T300**B(J))*EXP(-C(J)/gas_temperature)

      ! Check for temperature bounderies
      if (gas_temperature.LT.REACTION_TMIN(J)) reaction_rates(J)=A(J)*((REACTION_TMIN(J)/300.D0)**B(J))*EXP(-C(J)/REACTION_TMIN(J))
      if (gas_temperature.GT.REACTION_TMAX(J)) reaction_rates(J)=A(J)*((REACTION_TMAX(J)/300.D0)**B(J))*EXP(-C(J)/REACTION_TMAX(J))

      ! Check for the presence of several rate coefficients present in the network for the
      ! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)

        do M=1,W
          N=INDICE(M)
          !IF(IT==1) write(*,*) N,M, REACTION_SUBSTANCES_NAMES(:,N), REACTION_TMIN(N), REACTION_TMAX(N),distmin(M),distmax(M)
          if (gas_temperature.LT.REACTION_TMIN(N)) reaction_rates(N)=0.d0
          if (gas_temperature.GT.REACTION_TMAX(N)) reaction_rates(N)=0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then

          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=A(N)*((REACTION_TMIN(N)/300.D0)**B(N))*EXP(-C(N)/REACTION_TMIN(N))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=A(N)*((REACTION_TMAX(N)/300.D0)**B(N))*EXP(-C(N)/REACTION_TMAX(N))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  IONPOL1 FORMULA
    if (RATE_FORMULA(J).EQ.4) then
      reaction_rates(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/gas_temperature)**0.5))

      ! Check for temperature bounderies
      if (gas_temperature.LT.REACTION_TMIN(J)) reaction_rates(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/REACTION_TMIN(J))**0.5))
      if (gas_temperature.GT.REACTION_TMAX(J)) reaction_rates(J)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/REACTION_TMAX(J))**0.5))

      ! Check for the presence of several rate coefficients present in the network for the
      ! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)

        do M=1,W
          N=INDICE(M)
          if (gas_temperature.LT.REACTION_TMIN(N)) reaction_rates(N)= 0.d0
          if (gas_temperature.GT.REACTION_TMAX(N)) reaction_rates(N)= 0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/REACTION_TMIN(N))**0.5))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/REACTION_TMAX(N))**0.5))
          endif
        endif

        W=1
        INDICE(:)=0
        distmin(:)=9999.d0
        distmax(:)=9999.d0
      endif
    endif

    !---------------  IONPOL2 FORMULA
    if (RATE_FORMULA(J).EQ.5) then
      
      ! Check for temperature boundaries and apply the formula in consequence
      if (gas_temperature.LT.REACTION_TMIN(J)) then
        reaction_rates(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/REACTION_TMIN(J))**0.5)+(C(J)**2*300.d0/(10.526*REACTION_TMIN(J))))
      else if (gas_temperature.GT.REACTION_TMAX(J)) then
        reaction_rates(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/REACTION_TMAX(J))**0.5)+(C(J)**2*300.d0/(10.526*REACTION_TMAX(J))))
      else
        reaction_rates(J)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/gas_temperature)**0.5)+(C(J)**2*300.d0/(10.526*gas_temperature)))
      endif

      ! Check for the presence of several rate coefficients present in the network for the
      !! the same reaction
      if (REACTION_ID(J+1).EQ.REACTION_ID(J)) then
        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)
        W = W + 1
      endif

      if ((REACTION_ID(J+1).NE.REACTION_ID(J)).AND.(W.NE.1)) then

        INDICE(W)=J
        distmin(w)=REACTION_TMIN(j)-gas_temperature
        distmax(w)=gas_temperature-REACTION_TMAX(j)

        do M=1,W
          N=INDICE(M)
          if (gas_temperature.LT.REACTION_TMIN(N)) reaction_rates(N)= 0.d0
          if (gas_temperature.GT.REACTION_TMAX(N)) reaction_rates(N)= 0.d0
        enddo

        if (maxval(reaction_rates(indice(1:w))).lt.1.d-99) then
          if (minval(abs(distmin)).lt.minval(abs(distmax))) then
            N=indice(minloc(abs(distmin),dim=1))
            reaction_rates(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/REACTION_TMIN(N))**0.5) + &
                              (C(N)**2*300.d0/(10.526*REACTION_TMIN(N))))
          else
            N=indice(minloc(abs(distmax),dim=1))
            reaction_rates(N)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/REACTION_TMAX(N))**0.5) + &
                              (C(N)**2*300.d0/(10.526*REACTION_TMAX(N))))
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
      THERMAL_HOPING_RATE(K)=CHF(K)*EXP(-DIFFUSION_BARRIER(K)/dust_temperature)/nb_sites_per_grain
      CR_HOPING_RATE(K)=CHF(K)*EXP(-DIFFUSION_BARRIER(K)/PEAK_GRAIN_TEMPERATURE)/nb_sites_per_grain
      EVAPORATION_RATES(K)=CHF(K)*EXP(-DESORPTION_ENERGY(K)/dust_temperature)
      EVAPORATION_RATESCR(K)=CHF(K)*EXP(-DESORPTION_ENERGY(K)/PEAK_GRAIN_TEMPERATURE)
    enddo

    ! ========= Rxn ITYPE 15 - thermal evaporation
    ! ITYPE 15: Thermal evaporation
    do J=type_id_start(15),type_id_stop(15)
      reaction_rates(J)=A(J)*branching_ratio(J)*EVAPORATION_RATES(reagent_1_idx(J))
    enddo

    ! ========= Rxn ITYPE 16
    ! ITYPE 16: Cosmic-ray evaporation
    do J=type_id_start(16),type_id_stop(16)
      reaction_rates(J)=A(J)*branching_ratio(J)*((CR_IONISATION_RATE+X_IONISATION_RATE)/1.3D-17)&
      *CHF(reagent_1_idx(J))*CRFE*PEAK_DURATION*EXP(-DESORPTION_ENERGY(reagent_1_idx(J))/PEAK_GRAIN_TEMPERATURE)
    enddo


    if (type_id_start(98).ne.0) then
      ! ========= Rxn ITYPE 98 test the storage of H2S under a refractory form
      ! ITYPE 98: storage of H2S under a refractory form
      do j=type_id_start(98),type_id_stop(98)
        reaction_rates(J)=A(J)*(T300**B(J))*EXP(-C(J)*TI)
      enddo
    endif


    ! ====== Rxn ITYPE 17
    ! ITYPE 17: Photodissociations by Cosmic rays on grain surfaces
    ! Add X-rays in case X_IONISATION_RATE is not 0
    do J=type_id_start(17),type_id_stop(17)
      reaction_rates(J)=A(J)*(CR_IONISATION_RATE + X_IONISATION_RATE)
    enddo

    ! ====== Rxn ITYPE 18
    ! ITYPE 18: Photodissociations by Cosmic rays on grain surfaces
    do J=type_id_start(18),type_id_stop(18)
      reaction_rates(J)=A(J)*(CR_IONISATION_RATE+X_IONISATION_RATE)
    enddo

  endif

  ! When dust is turned off, zero all dust rates==========================
  if ((IS_GRAIN_REACTIONS.EQ.0).AND.(timestep.EQ.1)) then
    do J=type_id_start(14),type_id_stop(99)
      reaction_rates(J)=0.d0
      branching_ratio(J)=0.d0
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
  real(double_precision) :: ACTIV,BARR,DIFF,ACTIVCR,BARRCR,DIFFCR
  real(double_precision) :: XNH2,XNCO
  real(double_precision) :: TETABIS,TETABIS1,TETABIS2,TETABIS3
  real(double_precision) :: T300, TI, TSQ
  real(double_precision) :: YMOD1, YMOD2
  integer IMOD1,IMOD2
  integer :: j, l
  REAL(double_precision) :: XNDTOT         ! Sum of all abundances on grain surfaces
                                           ! Used to compute the photodesorption by FUV photons
                                           ! proposed by Hassel,G and following Oberg mesurments
  REAL(double_precision) :: MLAY           ! Number of layers from which species can desorb
  REAL(double_precision) :: SUMLAY         ! Total number of layers on the grain surface
  REAL(double_precision) :: UVCR           ! Scaling factor for CR generated UV
                                           ! The reference used is 1.3x10^-17 s-1


  T300=gas_temperature/300.D+0
  TI=1.0d00/gas_temperature
  TSQ=SQRT(gas_temperature)

  XNH2=Y(indH2)
  XNCO=Y(indCO)
  
  XNDTOT = 0.0E+00
  DO J = nb_gaseous_species+1,nb_species
    IF(species_name(J)(:1).NE.'J          ') PRINT*, "Warning: sum of all the species present on ", &
                                                  & "grain surface include gas-phase species"
    XNDTOT = XNDTOT + Y(J)
  ENDDO

  MLAY = 5.0E+00
  SUMLAY = XNDTOT*GTODN/nb_sites_per_grain
  UVCR = 1.300D-17 / CR_IONISATION_RATE

  !PRINT*, MLAY, SUMLAY


  ! Density/Av-dependent==================================================
  !
  ! VW 02/07 the treatment for the H2 and CO photodissociation was corrected
  ! for larger and smaller values of H2, CO column densities and Av not in the 
  ! data tables. I also added the choice of using this approximation in gg_control.d
  !
  ! ====== Rxn ITYPE 13
  ! ITYPE 2: Gas phase photodissociations/ionisations by UV
  do J=type_id_start(3),type_id_stop(3)
    reaction_rates(J)=A(J)*EXP(-C(J)*visual_extinction)*UV_FLUX

    ! MODIFY THE H2 AND CO PHOTODISSOCIATION if IS_ABSORPTION EQ 1
    if (IS_ABSORPTION.EQ.1) then 

      ! ====== Compute the H2 self-shielding
      if (REACTION_SUBSTANCES_NAMES(1,J).EQ.YH2) then
        TETABIS=1.D0

        NH2 = visual_extinction/AV_NH_ratio * XNH2

        ! ======= Linear extrapolation of the shielding factors
        do L=1,NL1-1
          if ((N1H2(L).LE.NH2).AND.(N1H2(L+1).GE.NH2)) then
            TETABIS=T1H2(L)+(NH2-N1H2(L))*(T1H2(L+1)-T1H2(L))/(N1H2(L+1)-N1H2(L))
          endif
        enddo
        if (NH2.GT.N1H2(NL1)) TETABIS = T1H2(NL1)

        reaction_rates(J)=2.54D-11*TETABIS

        reaction_rates(J)=reaction_rates(J)*UV_FLUX
      endif

      ! ====== Compute the CO self-shielding
      if (REACTION_SUBSTANCES_NAMES(1,J).EQ.YCO) then

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

        reaction_rates(J)=1.03D-10*TETABIS1*TETABIS2*TETABIS3
        reaction_rates(J)=reaction_rates(J)*UV_FLUX

      endif

    endif

  enddo

  ! Continually time-dependent grain rates================================
  if (IS_GRAIN_REACTIONS.NE.0) then

    ! ====== Rxn ITYPE 99
    ! ITYPE 99: Adsorption on grains
    do J=type_id_start(99),type_id_stop(99)
      ! ========= Set accretion rates
      ACCRETION_RATES(reagent_1_idx(J))=CONDSP(reagent_1_idx(J))*TSQ*Y(reagent_1_idx(J))*H_number_density
      ACCRETION_RATES(reagent_2_idx(J))=ACCRETION_RATES(reagent_1_idx(J))
      reaction_rates(J)=A(J)*branching_ratio(J)*ACCRETION_RATES(reagent_1_idx(J))/Y(reagent_1_idx(J))/GTODN
    enddo

    ! ====== Rxn ITYPE 14
    ! ITYPE 14: Grain surface reactions
    do J=type_id_start(14),type_id_stop(14)
      IMOD1=0
      IMOD2=0
      BARR=1.0d0
      ! --------- Calculate activation energy barrier multiplier
      if (ACTIVATION_ENERGY(J).GE.1.0D-40) then
        ACTIV=ACTIVATION_ENERGY(J)/dust_temperature
        ! ------------ Choose fastest of classical or tunnelling
        if (ACTIV.GT.quantum_activation_energy(J)) ACTIV=quantum_activation_energy(J)
        BARR=EXP(-ACTIV)
      endif

      ! --------- Thermal hopping diffusion method
      THERMAL_DIFFUSION_RATE_1(J)=THERMAL_HOPING_RATE(reagent_1_idx(J))
      THERMAL_DIFFUSION_RATE_2(J)=THERMAL_HOPING_RATE(reagent_2_idx(J))

      ! --------- Check for JH,JH2
      if (REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH)  IMOD1=1
      if (REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH2) IMOD1=2
      if (REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH)  IMOD2=1
      if (REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH2) IMOD2=2

      ! --------- QM for JH,JH2 only - others are too heavy
      if (IMOD1+IMOD2.NE.0) then
        ! ------------ QM1 - Tunnelling (if it's faster than thermal)
        if (GRAIN_TUNNELING_DIFFUSION.EQ.1) then
          if ((IMOD1.NE.0).AND.(TUNNELING_RATE_TYPE_1(reagent_1_idx(J)).GT.THERMAL_DIFFUSION_RATE_1(J))) then
            THERMAL_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_1(reagent_1_idx(J))
          endif
          if ((IMOD2.NE.0).AND.(TUNNELING_RATE_TYPE_1(reagent_2_idx(J)).GT.THERMAL_DIFFUSION_RATE_2(J))) then
            THERMAL_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_1(reagent_2_idx(J))
          endif
        endif
        ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
        if (GRAIN_TUNNELING_DIFFUSION.EQ.2) then
          if ((IMOD1.NE.0).AND.(TUNNELING_RATE_TYPE_2(reagent_1_idx(J)).GT.THERMAL_DIFFUSION_RATE_1(J))) then
            THERMAL_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_2(reagent_1_idx(J))
          endif
          if ((IMOD2.NE.0).AND.(TUNNELING_RATE_TYPE_2(reagent_2_idx(J)).GT.THERMAL_DIFFUSION_RATE_2(J))) then
            THERMAL_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_2(reagent_2_idx(J))
          endif
        endif
        ! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
        if (GRAIN_TUNNELING_DIFFUSION.EQ.3) then
          if (IMOD1.NE.0) then
            if (TUNNELING_RATE_TYPE_1(reagent_1_idx(J)).GT.THERMAL_DIFFUSION_RATE_1(J)) then
              THERMAL_DIFFUSION_RATE_1(J) = TUNNELING_RATE_TYPE_1(reagent_1_idx(J))
            endif
            if (TUNNELING_RATE_TYPE_2(reagent_1_idx(J)).GT.THERMAL_DIFFUSION_RATE_1(J)) then
              THERMAL_DIFFUSION_RATE_1(J) = TUNNELING_RATE_TYPE_2(reagent_1_idx(J))
            endif
          endif
          if (IMOD2.NE.0) then
            if (TUNNELING_RATE_TYPE_1(reagent_2_idx(J)).GT.THERMAL_DIFFUSION_RATE_2(J))  then
              THERMAL_DIFFUSION_RATE_2(J) = TUNNELING_RATE_TYPE_1(reagent_2_idx(J))
            endif
            if (TUNNELING_RATE_TYPE_2(reagent_2_idx(J)).GT.THERMAL_DIFFUSION_RATE_2(J))  then
              THERMAL_DIFFUSION_RATE_2(J) = TUNNELING_RATE_TYPE_2(reagent_2_idx(J))
            endif
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
          if ((REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH).OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JHe        ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JC         ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JN         ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JO         ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JS         ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JSi        ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JFe        ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JNa        ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JMg        ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JP         ').OR.&
          (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JCl        ')) IMOD1=3
          if ((REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH).OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JHe        ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JC         ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JN         ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JO         ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JS         ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JSi        ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JFe        ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JNa        ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JMg        ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JP         ').OR.&
          (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JCl        ')) IMOD2=3
        endif

        ! ------------ Modify rates (THERMAL_DIFFUSION_RATE_1 & THERMAL_DIFFUSION_RATE_2) according to their own evap/acc rates
        YMOD1=Y(reagent_1_idx(J))
        YMOD2=Y(reagent_2_idx(J))
        call modify_specific_rates(J,IMOD1,IMOD2,BARR,YMOD1,YMOD2)
      endif

      DIFF=THERMAL_DIFFUSION_RATE_1(J)+THERMAL_DIFFUSION_RATE_2(J)

      reaction_rates(J)=A(J)*branching_ratio(J)*BARR*DIFF*GTODN/H_number_density
      ! reaction_rates(J)=0.D0
    enddo

    ! Grain surface reactions (Cosmic Rays) DTEMP=70K
    ! ====== Rxn ITYPE 21
    ! ITYPE 21: Grain surface reactions
    if (type_id_start(21).NE.0) then
       do J=type_id_start(21),type_id_stop(21)
            IMOD1=0
            IMOD2=0
            BARRCR=1.0D+0
            ! --------- Calculate activation energy barrier multiplier
            IF (ACTIVATION_ENERGY(J).GE.1.0D-40) THEN
               ACTIVCR=ACTIVATION_ENERGY(J)/PEAK_GRAIN_TEMPERATURE
               ! ------------ Choose fastest of classical or tunnelling
               IF (ACTIVCR.GT.quantum_activation_energy(J)) ACTIVCR=quantum_activation_energy(J)
               BARRCR=EXP(-ACTIVCR)
            ENDIF

            ! --------- Thermal hopping diffusion method
            CR_DIFFUSION_RATE_1(J)=CR_HOPING_RATE(reagent_1_idx(J))
            CR_DIFFUSION_RATE_2(J)=CR_HOPING_RATE(reagent_2_idx(J))

            ! --------- Check for JH,JH2
            IF (REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH)  IMOD1=1
            IF (REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH2) IMOD1=2
            IF (REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH)  IMOD2=1
            IF (REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH2) IMOD2=2
            !IF (REACTION_SUBSTANCES_NAMES(1,J).EQ.YJO)  IMOD1=4
            !IF (REACTION_SUBSTANCES_NAMES(2,J).EQ.YJO)  IMOD2=4

            ! --------- QM for JH,JH2 only - others are too heavy
            IF (IMOD1+IMOD2.NE.0) THEN
              ! ------------ QM1 - Tunnelling (if it's faster than thermal)
              IF (grain_tunneling_diffusion.EQ.1) THEN
                if ((IMOD1.NE.0).AND.(TUNNELING_RATE_TYPE_1(reagent_1_idx(J)).GT.CR_DIFFUSION_RATE_1(J))) then
                  CR_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_1(reagent_1_idx(J))
                endif
                if ((IMOD2.NE.0).AND.(TUNNELING_RATE_TYPE_1(reagent_2_idx(J)).GT.CR_DIFFUSION_RATE_2(J))) then
                  CR_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_1(reagent_2_idx(J))
                endif
              ENDIF
              ! ------------ QM2 - Tunnelling: use estimated width of lowest energy band (if it's faster than thermal)
              IF (grain_tunneling_diffusion.EQ.2) THEN
                if ((IMOD1.NE.0).AND.(TUNNELING_RATE_TYPE_2(reagent_1_idx(J)).GT.CR_DIFFUSION_RATE_1(J))) then
                  CR_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_2(reagent_1_idx(J))
                endif
                if ((IMOD2.NE.0).AND.(TUNNELING_RATE_TYPE_2(reagent_2_idx(J)).GT.CR_DIFFUSION_RATE_2(J))) then
                  CR_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_2(reagent_2_idx(J))
                endif
              ENDIF
              ! ------------ QM3 - Fastest out of thermal, QM1, QM2 rates
              IF (grain_tunneling_diffusion.EQ.3) THEN
                IF (IMOD1.NE.0) THEN
                  if (TUNNELING_RATE_TYPE_1(reagent_1_idx(J)).GT.CR_DIFFUSION_RATE_1(J)) then
                    CR_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_1(reagent_1_idx(J))
                  endif
                  if (TUNNELING_RATE_TYPE_2(reagent_1_idx(J)).GT.CR_DIFFUSION_RATE_1(J)) then
                    CR_DIFFUSION_RATE_1(J)=TUNNELING_RATE_TYPE_2(reagent_1_idx(J))
                  endif
                ENDIF
                IF (IMOD2.NE.0) THEN
                  if (TUNNELING_RATE_TYPE_1(reagent_2_idx(J)).GT.CR_DIFFUSION_RATE_2(J)) then
                    CR_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_1(reagent_2_idx(J))
                  endif
                  if (TUNNELING_RATE_TYPE_2(reagent_2_idx(J)).GT.CR_DIFFUSION_RATE_2(J)) then
                    CR_DIFFUSION_RATE_2(J)=TUNNELING_RATE_TYPE_2(reagent_2_idx(J))
                  endif
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
                  IF ((REACTION_SUBSTANCES_NAMES(1,J).EQ.YJH).OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JHe        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JC         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JN         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JO         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JS         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JSi        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JFe        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JNa        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JMg        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JP         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(1,J).EQ.'JCl        ')) IMOD1=3
                  IF ((REACTION_SUBSTANCES_NAMES(2,J).EQ.YJH).OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JHe        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JC         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JN         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JO         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JS         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JSi        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JFe        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JNa        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JMg        ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JP         ').OR.&
                      (REACTION_SUBSTANCES_NAMES(2,J).EQ.'JCl        ')) IMOD2=3
               ENDIF

                  ! ------------ Modify rates (THERMAL_DIFFUSION_RATE_1 & THERMAL_DIFFUSION_RATE_2) according to their own evap/acc rates
                  YMOD1=Y(reagent_1_idx(J))
                  YMOD2=Y(reagent_2_idx(J))
                  CALL modify_specific_rates_cr(J,IMOD1,IMOD2,BARRCR,YMOD1,YMOD2)
            ENDIF

            DIFFCR=CR_DIFFUSION_RATE_1(J)+CR_DIFFUSION_RATE_2(J)

            reaction_rates(J)=(CR_IONISATION_RATE/1.3D-17)*CRFE*PEAK_DURATION*A(J)*branching_ratio(J)* &
                               BARRCR*DIFFCR*GTODN/H_number_density

       enddo
    endif

    ! ====== Rxn ITYPE 19 - 20
    ! ITYPE 19: Photodissociations by UV photons on grain surfaces
    ! ITYPE 20: Photodissociations by UV photons on grain surfaces
    do J=type_id_start(19),type_id_stop(20)
      reaction_rates(J)=A(J)*EXP(-C(J)*visual_extinction)*UV_FLUX
    enddo

! ============================================================
! Photodesorption, when used appears through ITYPES 66 and 67
! ============================================================
! ========= Rxn ITYPE 66
! ITYPE 66: Photodesorption by external UV
! 1.d8 is I_ISRF-FUV from Oberg et al. 2007, ApJ, 662, 23
    IF (type_id_start(66).NE.0) THEN
        DO J = type_id_start(66),type_id_stop(66)
!---- Used for all species
           reaction_rates(J)=A(J)/SITE_DENSITY*UV_FLUX*1.d8*EXP(-2.*visual_extinction)
!---- Specific cases
           CALL photodesorption_special_cases(J,SUMLAY)
!---- If there is more than MLAY on the grain surface, then we take into account that only
!     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
           IF(SUMLAY.GE.MLAY) reaction_rates(J) = reaction_rates(J) * MLAY / SUMLAY
        ENDDO
    ENDIF

! ========= Rxn ITYPE 67
! ITYPE 67: Photodesorption by CR generated UV
    IF (type_id_start(67).NE.0) THEN
        DO J = type_id_start(67),type_id_stop(67)
!---- Used for all species
           reaction_rates(J)=A(J)/SITE_DENSITY*1.d4*UVCR
           CALL photodesorption_special_cases(J,SUMLAY)
!---- If there is more than MLAY on the grain surface, then we take into account that only
!     the upper layers can photodesorb: this is done by assigning a reducing factor to the rate coefficient
           IF(SUMLAY.GE.MLAY) reaction_rates(J) = reaction_rates(J) * MLAY / SUMLAY
        ENDDO
    ENDIF

    ! Useful for testing
    ! To disable some reaction types
    !do j=type_id_start(14),type_id_stop(15)
    !reaction_rates(j)=0.
    !enddo

    !do j=type_id_start(17),type_id_stop(20)
    !reaction_rates(j)=0.
    !enddo
  endif

  ! Continually time-dependent gas phase rates============================
  ! H2 formation
  ! branching_ratio(1) and branching_ratio(2) are zero if IS_GRAIN_REACTIONS=1
  ! cf GRAINRATE
  ! VW Fev 2012 - this process has been removed
  !      do j=type_id_start(0),type_id_stop(0)
  !      if ((REACTION_SUBSTANCES_NAMES(1,J).eq.YH).and.(REACTION_SUBSTANCES_NAMES(2,j).eq.YH)) then
  !      reaction_rates(j)=branching_ratio(j)*A(j)*(T300**B(j))*GTODN/H_number_density/Y(reagent_1_idx(j))
  !      endif
  !      enddo

  ! if rate acoefficients are too small, put them to 0
  where (reaction_rates.lt.MINIMUM_RATE_COEFFICIENT) reaction_rates=0.d0

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

  EVAP_OVER_ACC_RATIO_1(J)=0.d0
  EVAP_OVER_ACC_RATIO_2(J)=0.d0

  ! --- Check value of x = t_acc/t_evap
  ! EVAPORATION_RATES = 1/t_evap
  ! ACCRETION_RATES = 1/t_acc
  if (ACCRETION_RATES(reagent_1_idx(J)).GT.0.d0) then
    EVAP_OVER_ACC_RATIO_1(J) = EVAPORATION_RATES(reagent_1_idx(J)) / ACCRETION_RATES(reagent_1_idx(J))
  endif
  if (ACCRETION_RATES(reagent_2_idx(J)).GT.0.d0) then
    EVAP_OVER_ACC_RATIO_2(J)=EVAPORATION_RATES(reagent_2_idx(J))/ACCRETION_RATES(reagent_2_idx(J))
  endif
  ! Hence x = 0 if t_evap or t_acc = 0

  ! --- Assign max rates

  if (BARR.EQ.1.0d0) then
    if (IMOD1.NE.0) then
      if (EVAP_OVER_ACC_RATIO_1(J).LT.1.0d0) then ! accretion dominates
        if (THERMAL_DIFFUSION_RATE_1(J).GT.ACCRETION_RATES(reagent_1_idx(J))) then
          THERMAL_DIFFUSION_RATE_1(J) = ACCRETION_RATES(reagent_1_idx(J))
        endif
      else ! evaporation dominates
        if (THERMAL_DIFFUSION_RATE_1(J).GT.EVAPORATION_RATES(reagent_1_idx(J))) then
          THERMAL_DIFFUSION_RATE_1(J) = EVAPORATION_RATES(reagent_1_idx(J))
        endif
      endif
    endif

    if (IMOD2.NE.0) then
      if (EVAP_OVER_ACC_RATIO_2(J).LT.1.0d0) then ! accretion dominates
        if (THERMAL_DIFFUSION_RATE_2(J).GT.ACCRETION_RATES(reagent_2_idx(J))) then
          THERMAL_DIFFUSION_RATE_2(J) = ACCRETION_RATES(reagent_2_idx(J))
        endif
      else ! evaporation dominates
        if (THERMAL_DIFFUSION_RATE_2(J).GT.EVAPORATION_RATES(reagent_2_idx(J))) then
          THERMAL_DIFFUSION_RATE_2(J) = EVAPORATION_RATES(reagent_2_idx(J))
        endif
      endif
    endif
  endif

  ! --- Species rate to compare chosen by fastest diffusion rate
  if (BARR.NE.1.0d0) then
    PICK=0.d0

    TESTREF1=ACCRETION_RATES(reagent_1_idx(J))
    if (EVAP_OVER_ACC_RATIO_1(J).GE.1.0d0) TESTREF1=EVAPORATION_RATES(reagent_1_idx(J))
    TESTREF2=ACCRETION_RATES(reagent_2_idx(J))
    if (EVAP_OVER_ACC_RATIO_2(J).GE.1.0d0) TESTREF2=EVAPORATION_RATES(reagent_2_idx(J))

    if (THERMAL_DIFFUSION_RATE_1(J).GE.THERMAL_DIFFUSION_RATE_2(J)) then
      TESTNUM=(THERMAL_DIFFUSION_RATE_1(J)+THERMAL_DIFFUSION_RATE_2(J))*BARR*YMOD2*GTODN
      if (YMOD2*GTODN.LT.1.0d0) TESTNUM=(THERMAL_DIFFUSION_RATE_1(J)+THERMAL_DIFFUSION_RATE_2(J))*BARR
      if (TESTNUM.GT.TESTREF1) PICK=1.d0
    endif
    if (THERMAL_DIFFUSION_RATE_2(J).GT.THERMAL_DIFFUSION_RATE_1(J)) then
      TESTNUM=(THERMAL_DIFFUSION_RATE_1(J)+THERMAL_DIFFUSION_RATE_2(J))*BARR*YMOD1*GTODN
      if (YMOD1*GTODN.LT.1.0d0) TESTNUM=(THERMAL_DIFFUSION_RATE_1(J)+THERMAL_DIFFUSION_RATE_2(J))*BARR
      if (TESTNUM.GT.TESTREF2) PICK=2.d0
    endif

    if (PICK.EQ.1) then
      THERMAL_DIFFUSION_RATE_1(J)=TESTREF1/BARR/YMOD2/GTODN
      if (YMOD2*GTODN.LT.1.0d0) THERMAL_DIFFUSION_RATE_1(J)=TESTREF1/BARR
      THERMAL_DIFFUSION_RATE_2(J)=0.d0
    endif

    if (PICK.EQ.2) then
      THERMAL_DIFFUSION_RATE_2(J)=TESTREF2/BARR/YMOD1/GTODN
      if (YMOD1*GTODN.LT.1.0d0) THERMAL_DIFFUSION_RATE_2(J)=TESTREF2/BARR
      THERMAL_DIFFUSION_RATE_1(J)=0.d0
    endif

  endif

  return
  end subroutine modify_specific_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Laura Reboussin
!
!> @date 2014
!
! DESCRIPTION:
!> @brief Modify some reaction rates on the grain surface (itype=21) in
!! various conditions. Test to estimates the fastest process and replace them.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  SUBROUTINE modify_specific_rates_cr(J,IMOD1,IMOD2,BARRCR,YMOD1,YMOD2)

  use global_variables

  implicit none

  integer :: J,IMOD1,IMOD2
  real(kind=8) :: BARRCR,YMOD1,YMOD2, PICK, TESTREF1, TESTREF2, TESTNUM

  EVAP_OVER_ACC_RATIO_1(J)=0.0D+0
  EVAP_OVER_ACC_RATIO_2(J)=0.0D+0

  ! --- Check value of x = t_acc/t_evap
  ! EVAPORATION_RATES = 1/t_evap
  ! ACCRETION_RATES = 1/t_acc
  IF (ACCRETION_RATES(reagent_1_idx(J)).GT.0.0D+0) THEN
     EVAP_OVER_ACC_RATIO_1(J) = EVAPORATION_RATESCR(reagent_1_idx(J)) / ACCRETION_RATES(reagent_1_idx(J))
  ENDIF
  IF (ACCRETION_RATES(reagent_2_idx(J)).GT.0.0D+0) THEN
     EVAP_OVER_ACC_RATIO_2(J)=EVAPORATION_RATESCR(reagent_2_idx(J))/ACCRETION_RATES(reagent_2_idx(J))
  ENDIF
  ! Hence x = 0 if t_evap or t_acc = 0

  ! --- Assign max rates

  IF (BARRCR.EQ.1.0D+0) THEN
     IF (IMOD1.NE.0) THEN
        if (EVAP_OVER_ACC_RATIO_1(J).LT.1.0D+0) THEN ! accretion dominates
          if (CR_DIFFUSION_RATE_1(J).GT.ACCRETION_RATES(reagent_1_idx(J))) then
            CR_DIFFUSION_RATE_1(J) = ACCRETION_RATES(reagent_1_idx(J))
          endif
        else ! evaporation dominates
          if (CR_DIFFUSION_RATE_1(J).GT.EVAPORATION_RATESCR(reagent_1_idx(J))) then
            CR_DIFFUSION_RATE_1(J) = EVAPORATION_RATESCR(reagent_1_idx(J))
          endif
        endif
     ENDIF

     IF (IMOD2.NE.0) THEN
        if (EVAP_OVER_ACC_RATIO_2(J).LT.1.0D+0) THEN ! accretion dominates
          if (CR_DIFFUSION_RATE_2(J).GT.ACCRETION_RATES(reagent_2_idx(J))) then
            CR_DIFFUSION_RATE_2(J) = ACCRETION_RATES(reagent_2_idx(J))
          endif
        else ! evaporation dominates
          if (CR_DIFFUSION_RATE_2(J).GT.EVAPORATION_RATESCR(reagent_2_idx(J))) then
            CR_DIFFUSION_RATE_2(J) = EVAPORATION_RATESCR(reagent_2_idx(J))
          endif
        endif
     ENDIF
  ENDIF

  ! --- Species rate to compare chosen by fastest diffusion rate
  IF (BARRCR.NE.1.0D+0) THEN
     PICK=0

     TESTREF1=ACCRETION_RATES(reagent_1_idx(J))
     IF (EVAP_OVER_ACC_RATIO_1(J).GE.1.0D+0) TESTREF1=EVAPORATION_RATESCR(reagent_1_idx(J))
     TESTREF2=ACCRETION_RATES(reagent_2_idx(J))
     IF (EVAP_OVER_ACC_RATIO_2(J).GE.1.0D+0) TESTREF2=EVAPORATION_RATESCR(reagent_2_idx(J))

     IF (CR_DIFFUSION_RATE_1(J).GE.CR_DIFFUSION_RATE_2(J)) THEN
        TESTNUM=(CR_DIFFUSION_RATE_1(J)+CR_DIFFUSION_RATE_2(J))*BARRCR*YMOD2*GTODN
        IF (YMOD2*GTODN.LT.1.0D+0) TESTNUM=(CR_DIFFUSION_RATE_1(J)+CR_DIFFUSION_RATE_2(J))*BARRCR
        IF (TESTNUM.GT.TESTREF1) PICK=1
     ENDIF
     IF (CR_DIFFUSION_RATE_2(J).GT.CR_DIFFUSION_RATE_1(J)) THEN
        TESTNUM=(CR_DIFFUSION_RATE_1(J)+CR_DIFFUSION_RATE_2(J))*BARRCR*YMOD1*GTODN
        IF (YMOD1*GTODN.LT.1.0D+0) TESTNUM=(CR_DIFFUSION_RATE_1(J)+CR_DIFFUSION_RATE_2(J))*BARRCR
        IF (TESTNUM.GT.TESTREF2) PICK=2
     ENDIF

     IF (PICK.EQ.1) THEN
        CR_DIFFUSION_RATE_1(J)=TESTREF1/BARRCR/YMOD2/GTODN
        IF (YMOD2*GTODN.LT.1.0D+0) CR_DIFFUSION_RATE_1(J)=TESTREF1/BARRCR
        CR_DIFFUSION_RATE_2(J)=0.0D+0
     ENDIF

     IF (PICK.EQ.2) THEN
        CR_DIFFUSION_RATE_2(J)=TESTREF2/BARRCR/YMOD1/GTODN
        IF (YMOD1*GTODN.LT.1.0D+0) CR_DIFFUSION_RATE_2(J)=TESTREF2/BARRCR
        CR_DIFFUSION_RATE_1(J)=0.0D+0
     ENDIF

  ENDIF

RETURN
END SUBROUTINE modify_specific_rates_cr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Treat special cases of photodesorption such as CO2, CO, H2O, CH3OH and N2
!!\n Data from Oberg et al. 2009:
!!\n      a - A&A, 496, 281-293
!!\n      b - ApJ, 693, 1209-1218
!!\n      c - A&A, 504, 891-913
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SUBROUTINE photodesorption_special_cases(J,SUMLAY)
! Treat special cases of photodesorption such as CO2, CO, H2O, CH3OH and N2
! Data from Oberg et al. 2009:
!      a - A&A, 496, 281-293
!      b - ApJ, 693, 1209-1218
!      c - A&A, 504, 891-913

use global_variables

IMPLICIT NONE

INTEGER, INTENT(IN) :: J
REAL(double_precision), INTENT(IN) :: SUMLAY
REAL(double_precision) :: LTD
REAL(double_precision) :: fH2O

!------ Photodesorption of CO2: photodesorbs as either CO2 or CO
IF (dust_temperature.LE.3.5E+01) THEN
    IF (REACTION_SUBSTANCES_NAMES(4,J) == 'CO2        ') THEN
        reaction_rates(J) = reaction_rates(J) * 1.2E-03 * ( 1.0E+00 - EXP(-SUMLAY/2.9E+00) ) / A(J)
    ENDIF
    IF((REACTION_SUBSTANCES_NAMES(4,J) == 'CO         ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'O          ') .OR. &
       (REACTION_SUBSTANCES_NAMES(4,J) == 'O          ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'CO         ')) THEN
        reaction_rates(J) = reaction_rates(J) * 1.1E-03 * ( 1.0E+00 - EXP(-SUMLAY/4.6E+00) ) / A(J)
    ENDIF
ELSEIF(dust_temperature.GT.3.5E+01) THEN
    IF (REACTION_SUBSTANCES_NAMES(4,J) == 'CO2        ') THEN
        reaction_rates(J) = reaction_rates(J) * 2.2E-03 * ( 1.0E+00 - EXP(-SUMLAY/5.8E+00) ) / A(J)
    ENDIF
    IF((REACTION_SUBSTANCES_NAMES(4,J) == 'CO         ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'O          ') .OR. &
       (REACTION_SUBSTANCES_NAMES(4,J) == 'O          ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'CO         ')) THEN
       reaction_rates(J) = reaction_rates(J) * 2.2E-04 * SUMLAY / A(J)
    ENDIF
ENDIF
!------ Photodesorption of CO
IF(REACTION_SUBSTANCES_NAMES(1,J) == 'JCO        ' .AND. REACTION_SUBSTANCES_NAMES(4,J) == 'CO         ') THEN
   reaction_rates(J) = reaction_rates(J) * ( 2.7E-03 - 1.7E-04 * (dust_temperature - 15E+00)) / A(J)
ENDIF
!------ Photodesorption of N2
IF(REACTION_SUBSTANCES_NAMES(1,J) == 'JN2        ' .AND. REACTION_SUBSTANCES_NAMES(4,J) == 'N2         ') THEN
   reaction_rates(J) = reaction_rates(J) * 4.0E-04 / A(J)
ENDIF
!------ Photodesorption of CH3OH
IF(REACTION_SUBSTANCES_NAMES(1,J) == 'JCH3OH     ' .AND. REACTION_SUBSTANCES_NAMES(4,J) == 'CH3OH      ') THEN
   reaction_rates(J) = reaction_rates(J) * 2.1E-03 / A(J)
ENDIF
!------ Photodesorption of H2O
IF(REACTION_SUBSTANCES_NAMES(1,J) == 'JH2O       ') THEN
   LTD = 6.0E-01 + 2.4E-02 * dust_temperature
   fH2O = 4.2E-01 + 2.0E-03 * dust_temperature
   reaction_rates(J) = reaction_rates(J) * 1.0E-03 * (1.3E+00 + 3.2E-02*dust_temperature) &
           & * ( 1.0E+00 - EXP(-SUMLAY/LTD) ) * fH2O / A(J)
   IF((REACTION_SUBSTANCES_NAMES(4,J) == 'OH         ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'H          ') .OR. &
      (REACTION_SUBSTANCES_NAMES(4,J) == 'H          ' .AND. REACTION_SUBSTANCES_NAMES(5,J) == 'OH         ')) THEN
      reaction_rates(J) = reaction_rates(J) * (1.0E+00 - fH2O)/fH2O
   ENDIF
ENDIF

END SUBROUTINE photodesorption_special_cases

! ======================================================================

end module ode_solver