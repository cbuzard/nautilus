!******************************************************************************
! MODULE: dynamic
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that solve an Ordinary Differential Equation using the 
!! Livermore solver. All routines comes from the ODEPACK package, but correspond
!! only to the DLSODES solver. 
!
!> @warning implicite none is not present so far. Error occurs during the 
!! compilation because some routine need integer and get REAL. I don't know how
!! to solve this.
!
!******************************************************************************

module ode_solve


!~   implicit none
  
!~   private
!~   
!~   public :: dlsodes
  
  contains

SUBROUTINE DLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
EXTERNAL F, JAC
INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
DIMENSION Y(*), ATOL(*), RWORK(LRW), IWORK(LIW)
!-----------------------------------------------------------------------
! This is the 12 November 2003 version of
! DLSODES: Livermore Solver for Ordinary Differential Equations
!          with general Sparse Jacobian matrix.
!
! This version is in double precision.
!
! DLSODES solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DLSODES is a variant of the DLSODE package, and is intended for
! problems in which the Jacobian matrix df/dy has an arbitrary
! sparse structure (when the problem is stiff).
!
! Authors:       Alan C. Hindmarsh
!                Center for Applied Scientific Computing, L-561
!                Lawrence Livermore National Laboratory
!                Livermore, CA 94551
! and
!                Andrew H. Sherman
!                J. S. Nolen and Associates
!                Houston, TX 77084
!-----------------------------------------------------------------------
! References:
! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
!     North-Holland, Amsterdam, 1983, pp. 55-64.
!
! 2.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
!     Yale Sparse Matrix Package: I. The Symmetric Codes,
!     Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
!
! 3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
!     Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
!     Research Report No. 114, Dept. of Computer Sciences, Yale
!     University, 1977.
!-----------------------------------------------------------------------
! Summary of Usage.
!
! Communication between the user and the DLSODES package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form:
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are two standard
! choices for the method flag, MF = 121 and MF = 222.  In both cases,
! DLSODES requires the Jacobian matrix in some form, and it treats this
! matrix in general sparse form, with sparsity structure determined
! internally.  (For options where the user supplies the sparsity
! structure, see the full description of MF below.)
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 121), but if this is not feasible, DLSODES will
! compute it internally by difference quotients (MF = 222).
! If you are supplying the Jacobian, provide a subroutine of the form:
!               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
!               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
! Here NEQ, T, Y, and J are input arguments, and the JAC routine is to
! load the array PDJ (of length NEQ) with the J-th column of df/dy.
! I.e., load PDJ(i) with df(i)/dy(J) for all relevant values of i.
! The arguments IAN and JAN should be ignored for normal situations.
! DLSODES will call the JAC routine with J = 1,2,...,NEQ.
! Only nonzero elements need be loaded.  Usually, a crude approximation
! to df/dy, possibly with fewer nonzero elements, will suffice.
!
! D. Write a main program which calls Subroutine DLSODES once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages by
! DLSODES.  On the first call to DLSODES, supply arguments as follows:
! F      = name of subroutine for right-hand side vector f.
!          This name must be declared External in calling program.
! NEQ    = number of first order ODEs.
! Y      = array of initial values, of length NEQ.
! T      = the initial value of the independent variable t.
! TOUT   = first point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = relative tolerance parameter (scalar).
! ATOL   = absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional inputs used.
! RWORK  = real work array of length at least:
!             20 + 16*NEQ            for MF = 10,
!             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
!                                    for MF = 121 or 222,
!          where:
!          NNZ    = the number of nonzero elements in the sparse
!                   Jacobian (if this is unknown, use an estimate), and
!          LENRAT = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          In any case, the required size of RWORK cannot generally
!          be predicted in advance if MF = 121 or 222, and the value
!          above is a rough estimate of a crude lower bound.  Some
!          experimentation with this size may be necessary.
!          (When known, the correct required length is an optional
!          output, available in IWORK(17).)
! LRW    = declared length of RWORK (in user dimension).
! IWORK  = integer work array of length at least 30.
! LIW    = declared length of IWORK (in user dimension).
! JAC    = name of subroutine for Jacobian matrix (MF = 121).
!          If used, this name must be declared External in calling
!          program.  If not used, pass a dummy name.
! MF     = method flag.  Standard values are:
!          10  for nonstiff (Adams) method, no Jacobian used
!          121 for stiff (BDF) method, user-supplied sparse Jacobian
!          222 for stiff method, internally generated sparse Jacobian
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL.
!
! E. The output from the first call (or any call) is:
!      Y = array of computed values of y(t) vector.
!      T = corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DLSODES was successful, negative otherwise.
!          -1 means excess work done on this call (perhaps wrong MF).
!          -2 means excess accuracy requested (tolerances too small).
!          -3 means illegal input detected (see printed message).
!          -4 means repeated error test failures (check all inputs).
!          -5 means repeated convergence failures (perhaps bad Jacobian
!             supplied or wrong choice of MF or tolerances).
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!          -7 means a fatal error return flag came from sparse solver
!             CDRV by way of DPRJS or DSOLSS.  Should never happen.
!          A return with ISTATE = -1, -4, or -5 may result from using
!          an inappropriate sparsity structure, one that is quite
!          different from the initial structure.  Consider calling
!          DLSODES again with ISTATE = 3 to force the structure to be
!          reevaluated.  See the full description of ISTATE below.
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call DLSODES again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! Example Problem.
!
! The following is a simple example problem, with the coding
! needed for its solution by DLSODES.  The problem is from chemical
! kinetics, and consists of the following 12 rate equations:
!    dy1/dt  = -rk1*y1
!    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
!                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
!    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
!                + rk11*rk14*y4 + rk12*rk14*y6
!    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
!    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
!    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
!    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
!    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
!    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
!    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
!                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
!                - rk6*y10 - rk9*y10
!    dy11/dt = rk10*y8
!    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
!                - rk15*y2*y12 - rk17*y10*y12
!
! with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
!      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
!      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
!      rk15 = rk17 = 100.0.
!
! The t interval is from 0 to 1000, and the initial conditions
! are y1 = 1, y2 = y3 = ... = y12 = 0.  The problem is stiff.
!
! The following coding solves this problem with DLSODES, using MF = 121
! and printing results at t = .1, 1., 10., 100., 1000.  It uses
! ITOL = 1 and mixed relative/absolute tolerance controls.
! During the run and at the end, statistical quantities of interest
! are printed (see optional outputs in the full description below).
!
!     EXTERNAL FEX, JEX
!     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
!     DIMENSION Y(12), RWORK(500), IWORK(30)
!     DATA LRW/500/, LIW/30/
!     NEQ = 12
!     DO 10 I = 1,NEQ
! 10    Y(I) = 0.0D0
!     Y(1) = 1.0D0
!     T = 0.0D0
!     TOUT = 0.1D0
!     ITOL = 1
!     RTOL = 1.0D-4
!     ATOL = 1.0D-6
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     MF = 121
!     DO 40 IOUT = 1,5
!       CALL DLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
!    1     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
!       WRITE(6,30)T,IWORK(11),RWORK(11),(Y(I),I=1,NEQ)
! 30    FORMAT(//' At t =',D11.3,4X,
!    1    ' No. steps =',I5,4X,' Last step =',D11.3/
!    2    '  Y array =  ',4D14.5/13X,4D14.5/13X,4D14.5)
!       IF (ISTATE .LT. 0) GOTO 80
!       TOUT = TOUT*10.0D0
! 40    CONTINUE
!     LENRW = IWORK(17)
!     LENIW = IWORK(18)
!     NST = IWORK(11)
!     NFE = IWORK(12)
!     NJE = IWORK(13)
!     NLU = IWORK(21)
!     NNZ = IWORK(19)
!     NNZLU = IWORK(25) + IWORK(26) + NEQ
!     WRITE (6,70) LENRW,LENIW,NST,NFE,NJE,NLU,NNZ,NNZLU
! 70  FORMAT(//' Required RWORK size =',I4,'   IWORK size =',I4/
!    1   ' No. steps =',I4,'   No. f-s =',I4,'   No. J-s =',I4,
!    2   '   No. LU-s =',I4/' No. of nonzeros in J =',I5,
!    3   '   No. of nonzeros in LU =',I5)
!     STOP
! 80  WRITE(6,90)ISTATE
! 90  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)
!     DOUBLE PRECISION T, Y, YDOT
!     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
!    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
!     DIMENSION Y(12), YDOT(12)
!     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
!    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
!    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
!    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
!    4   RK19/50.0D0/, RK20/50.0D0/
!     YDOT(1)  = -RK1*Y(1)
!     YDOT(2)  = RK1*Y(1) + RK11*RK14*Y(4) + RK19*RK14*Y(5)
!    1           - RK3*Y(2)*Y(3) - RK15*Y(2)*Y(12) - RK2*Y(2)
!     YDOT(3)  = RK2*Y(2) - RK5*Y(3) - RK3*Y(2)*Y(3) - RK7*Y(10)*Y(3)
!    1           + RK11*RK14*Y(4) + RK12*RK14*Y(6)
!     YDOT(4)  = RK3*Y(2)*Y(3) - RK11*RK14*Y(4) - RK4*Y(4)
!     YDOT(5)  = RK15*Y(2)*Y(12) - RK19*RK14*Y(5) - RK16*Y(5)
!     YDOT(6)  = RK7*Y(10)*Y(3) - RK12*RK14*Y(6) - RK8*Y(6)
!     YDOT(7)  = RK17*Y(10)*Y(12) - RK20*RK14*Y(7) - RK18*Y(7)
!     YDOT(8)  = RK9*Y(10) - RK13*RK14*Y(8) - RK10*Y(8)
!     YDOT(9)  = RK4*Y(4) + RK16*Y(5) + RK8*Y(6) + RK18*Y(7)
!     YDOT(10) = RK5*Y(3) + RK12*RK14*Y(6) + RK20*RK14*Y(7)
!    1           + RK13*RK14*Y(8) - RK7*Y(10)*Y(3) - RK17*Y(10)*Y(12)
!    2           - RK6*Y(10) - RK9*Y(10)
!     YDOT(11) = RK10*Y(8)
!     YDOT(12) = RK6*Y(10) + RK19*RK14*Y(5) + RK20*RK14*Y(7)
!    1           - RK15*Y(2)*Y(12) - RK17*Y(10)*Y(12)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, J, IA, JA, PDJ)
!     DOUBLE PRECISION T, Y, PDJ
!     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
!    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
!     DIMENSION Y(12), IA(*), JA(*), PDJ(12)
!     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
!    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
!    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
!    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
!    4   RK19/50.0D0/, RK20/50.0D0/
!     GOTO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), J
! 1   PDJ(1) = -RK1
!     PDJ(2) = RK1
!     RETURN
! 2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
!     PDJ(3) = RK2 - RK3*Y(3)
!     PDJ(4) = RK3*Y(3)
!     PDJ(5) = RK15*Y(12)
!     PDJ(12) = -RK15*Y(12)
!     RETURN
! 3   PDJ(2) = -RK3*Y(2)
!     PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
!     PDJ(4) = RK3*Y(2)
!     PDJ(6) = RK7*Y(10)
!     PDJ(10) = RK5 - RK7*Y(10)
!     RETURN
! 4   PDJ(2) = RK11*RK14
!     PDJ(3) = RK11*RK14
!     PDJ(4) = -RK11*RK14 - RK4
!     PDJ(9) = RK4
!     RETURN
! 5   PDJ(2) = RK19*RK14
!     PDJ(5) = -RK19*RK14 - RK16
!     PDJ(9) = RK16
!     PDJ(12) = RK19*RK14
!     RETURN
! 6   PDJ(3) = RK12*RK14
!     PDJ(6) = -RK12*RK14 - RK8
!     PDJ(9) = RK8
!     PDJ(10) = RK12*RK14
!     RETURN
! 7   PDJ(7) = -RK20*RK14 - RK18
!     PDJ(9) = RK18
!     PDJ(10) = RK20*RK14
!     PDJ(12) = RK20*RK14
!     RETURN
! 8   PDJ(8) = -RK13*RK14 - RK10
!     PDJ(10) = RK13*RK14
!     PDJ(11) = RK10
! 9   RETURN
! 10  PDJ(3) = -RK7*Y(3)
!     PDJ(6) = RK7*Y(3)
!     PDJ(7) = RK17*Y(12)
!     PDJ(8) = RK9
!     PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
!     PDJ(12) = RK6 - RK17*Y(12)
! 11  RETURN
! 12  PDJ(2) = -RK15*Y(2)
!     PDJ(5) = RK15*Y(2)
!     PDJ(7) = RK17*Y(10)
!     PDJ(10) = -RK17*Y(10)
!     PDJ(12) = -RK15*Y(2) - RK17*Y(10)
!     RETURN
!     END
!
! The output of this program (on a Cray-1 in single precision)
! is as follows:
!
!
! At t =  1.000e-01     No. steps =   12     Last step =  1.515e-02
!  Y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
!                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
!                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
!
!
! At t =  1.000e+00     No. steps =   33     Last step =  7.880e-02
!  Y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
!                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
!                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
!
!
! At t =  1.000e+01     No. steps =   48     Last step =  1.239e+00
!  Y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
!                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
!                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
!
!
! At t =  1.000e+02     No. steps =   91     Last step =  3.764e+00
!  Y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
!                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
!                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
!
!
! At t =  1.000e+03     No. steps =  111     Last step =  4.156e+02
!  Y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
!               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
!                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
!
!
! Required RWORK size = 442   IWORK size =  30
! No. steps = 111   No. f-s = 142   No. J-s =   2   No. LU-s =  20
! No. of nonzeros in J =   44   No. of nonzeros in LU =   50
!
!-----------------------------------------------------------------------
! Full Description of User Interface to DLSODES.
!
! The user interface to DLSODES consists of the following parts.
!
! 1.   The call sequence to Subroutine DLSODES, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is a description of
!      optional inputs available through the call sequence, and then
!      a description of optional outputs (in the work arrays).
!
! 2.   Descriptions of other routines in the DLSODES package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      Common, and obtain specified derivatives of the solution y(t).
!
! 3.   Descriptions of Common blocks to be declared in overlay
!      or similar environments, or to be saved when doing an interrupt
!      of the problem and continued solution later.
!
! 4.   Description of two routines in the DLSODES package, either of
!      which the user may replace with his/her own version, if desired.
!      These relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part 1.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional inputs and optional outputs.  (The term output here refers
! to the return from Subroutine DLSODES to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 on input.
!
! The descriptions of the call arguments are as follows.
!
! F      = the name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter y(1),...,y(NEQ).
!          F must be declared External in the calling program.
!
!          Subroutine F may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in F) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y below.
!
!          If quantities computed in the F routine are needed
!          externally to DLSODES, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DINTDY instead.
!
! NEQ    = the size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may be decreased, but not increased, during the problem.
!          If NEQ is decreased (with ISTATE = 3 on input), the
!          remaining components of Y should be left undisturbed, if
!          these are to be accessed in F and/or JAC.
!
!          Normally, NEQ is a scalar, and it is generally referred to
!          as a scalar in this user interface description.  However,
!          NEQ may be an array, with NEQ(1) set to the system size.
!          (The DLSODES package accesses only NEQ(1).)  In either case,
!          this parameter is passed as the NEQ argument in all calls
!          to F and JAC.  Hence, if it is an array, locations
!          NEQ(2),... may be used to store other integer data and pass
!          it to F and/or JAC.  Subroutines F and/or JAC must include
!          NEQ in a Dimension statement in that case.
!
! Y      = a real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          on the first call, Y must contain the vector of initial
!          values.  On output, Y contains the computed solution vector,
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.  Hence its length may exceed NEQ, and locations
!          Y(NEQ+1),... may be used to store other real data and
!          pass it to F and/or JAC.  (The DLSODES package accesses only
!          Y(1),...,Y(NEQ).)
!
! T      = the independent variable.  On input, T is used only on the
!          first call, as the initial point of the integration.
!          on output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = the next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal T interval, whose endpoints are
!          TCUR - HU and TCUR (see optional outputs, below, for
!          TCUR and HU).
!
! ITOL   = an indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = a relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = an absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!             The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector E = (E(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!          and the RMS-norm (root-mean-square norm) here is
!          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part 4 below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = an index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at t = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          On input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
!             the conditional inputs IA and JA,
!             and any of the optional inputs except H0.
!             In particular, if MITER = 1 or 2, a call with ISTATE = 3
!             will cause the sparsity structure of the problem to be
!             recomputed (or reread from IA and JA if MOSS = 0).
!          Note:  a preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful for the
!          purpose of outputting the initial conditions.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 on input.
!
!          On output, ISTATE has the following values and meanings.
!           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again
!              (the excess work step counter will be reset to 0).
!              In addition, the user may increase MXSTEP to avoid
!              this error return (see below on optional inputs).
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          -7  means a fatal error return flag came from the sparse
!              solver CDRV by way of DPRJS or DSOLSS (numerical
!              factorization or backsolve).  This should never happen.
!              The integration was successful as far as T.
!
!          Note: an error return with ISTATE = -1, -4, or -5 and with
!          MITER = 1 or 2 may mean that the sparsity structure of the
!          problem has changed significantly since it was last
!          determined (or input).  In that case, one can attempt to
!          complete the integration by setting ISTATE = 3 on the next
!          call, so that a new structure determination is done.
!
!          Note:  since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other inputs, before
!          calling the solver again.
!
! IOPT   = an integer flag to specify whether or not any optional
!          inputs are being used on this call.  Input only.
!          The optional inputs are listed separately below.
!          IOPT = 0 means no optional inputs are being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means one or more optional inputs are being used.
!
! RWORK  = a work array used for a mixture of real (double precision)
!          and integer work space.
!          The length of RWORK (in real words) must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = 0                                    if MITER = 0,
!          LWM = 2*NNZ + 2*NEQ + (NNZ+9*NEQ)/LENRAT   if MITER = 1,
!          LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT  if MITER = 2,
!          LWM = NEQ + 2                              if MITER = 3.
!          In the above formulas,
!          NNZ    = number of nonzero elements in the Jacobian matrix.
!          LENRAT = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          the minimum length of RWORK is:
!             20 + 16*NEQ        for MF = 10,
!             20 + 16*NEQ + LWM  for MF = 11, 111, 211, 12, 112, 212,
!             22 + 17*NEQ        for MF = 13,
!             20 +  9*NEQ        for MF = 20,
!             20 +  9*NEQ + LWM  for MF = 21, 121, 221, 22, 122, 222,
!             22 + 10*NEQ        for MF = 23.
!          If MITER = 1 or 2, the above formula for LWM is only a
!          crude lower bound.  The required length of RWORK cannot
!          be readily predicted in general, as it depends on the
!          sparsity structure of the problem.  Some experimentation
!          may be necessary.
!
!          The first 20 words of RWORK are reserved for conditional
!          and optional inputs and optional outputs.
!
!          The following word in RWORK is a conditional input:
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = the length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = an integer work array.  The length of IWORK must be at least
!             31 + NEQ + NNZ   if MOSS = 0 and MITER = 1 or 2, or
!             30               otherwise.
!          (NNZ is the number of nonzero elements in df/dy.)
!
!          In DLSODES, IWORK is used only for conditional and
!          optional inputs and optional outputs.
!
!          The following two blocks of words in IWORK are conditional
!          inputs, required if MOSS = 0 and MITER = 1 or 2, but not
!          otherwise (see the description of MF for MOSS).
!            IWORK(30+j) = IA(j)     (j=1,...,NEQ+1)
!            IWORK(31+NEQ+k) = JA(k) (k=1,...,NNZ)
!          The two arrays IA and JA describe the sparsity structure
!          to be assumed for the Jacobian matrix.  JA contains the row
!          indices where nonzero elements occur, reading in columnwise
!          order, and IA contains the starting locations in JA of the
!          descriptions of columns 1,...,NEQ, in that order, with
!          IA(1) = 1.  Thus, for each column index j = 1,...,NEQ, the
!          values of the row index i in column j where a nonzero
!          element may occur are given by
!            i = JA(k),  where   IA(j) .le. k .lt. IA(j+1).
!          If NNZ is the total number of nonzero locations assumed,
!          then the length of the JA array is NNZ, and IA(NEQ+1) must
!          be NNZ + 1.  Duplicate entries are not allowed.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note:  The work arrays must not be altered between calls to DLSODES
! for the same problem, except possibly for the conditional and
! optional inputs, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DLSODES between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = name of user-supplied routine (MITER = 1 or MOSS = 1) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
!               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
!          where NEQ, T, Y, J, IAN, and JAN are input, and the array
!          PDJ, of length NEQ, is to be loaded with column J
!          of the Jacobian on output.  Thus df(i)/dy(J) is to be
!          loaded into PDJ(i) for all relevant values of i.
!          Here T and Y have the same meaning as in Subroutine F,
!          and J is a column index (1 to NEQ).  IAN and JAN are
!          undefined in calls to JAC for structure determination
!          (MOSS = 1).  otherwise, IAN and JAN are structure
!          descriptors, as defined under optional outputs below, and
!          so can be used to determine the relevant row indices i, if
!          desired.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with greater sparsity) will do.
!               In any case, PDJ is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Calls to JAC are made with J = 1,...,NEQ, in that order, and
!          each such set of calls is preceded by a call to F with the
!          same arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user Common block by F and not recomputed by JAC,
!          if desired.  JAC must not alter its input arguments.
!          JAC must be declared External in the calling program.
!               Subroutine JAC may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y above.
!
! MF     = the method flag.  Used only for input.
!          MF has three decimal digits-- MOSS, METH, MITER--
!             MF = 100*MOSS + 10*METH + MITER.
!          MOSS indicates the method to be used to obtain the sparsity
!          structure of the Jacobian matrix if MITER = 1 or 2:
!            MOSS = 0 means the user has supplied IA and JA
!                     (see descriptions under IWORK above).
!            MOSS = 1 means the user has supplied JAC (see below)
!                     and the structure will be obtained from NEQ
!                     initial calls to JAC.
!            MOSS = 2 means the structure will be obtained from NEQ+1
!                     initial calls to F.
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on Backward
!                     Differentiation Formulas (BDFs).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      sparse Jacobian, given by Subroutine JAC.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) sparse Jacobian
!                      (using NGP extra calls to F per df/dy value,
!                      where NGP is an optional output described below.)
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!          If MITER = 1 or MOSS = 1, the user must supply a Subroutine
!          JAC (the name is arbitrary) as described above under JAC.
!          Otherwise, a dummy argument can be used.
!
!          The standard choices for MF are:
!            MF = 10  for a nonstiff problem,
!            MF = 21 or 22 for a stiff problem with IA/JA supplied
!                     (21 if JAC is supplied, 22 if not),
!            MF = 121 for a stiff problem with JAC supplied,
!                     but not IA/JA,
!            MF = 222 for a stiff problem with neither IA/JA nor
!                     JAC supplied.
!          The sparseness structure can be changed during the
!          problem by making a call to DLSODES with ISTATE = 3.
!-----------------------------------------------------------------------
! Optional Inputs.
!
! The following is a list of the optional inputs provided for in the
! call sequence.  (See also Part 2.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of these inputs requires IOPT = 1, and in that
! case all of these inputs are examined.  A value of zero for any
! of these optional inputs will cause the default value to be used.
! Thus to use a subset of the optional inputs, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! Name    Location      Meaning and Default Value
!
! H0      RWORK(5)  the step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  the maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  the minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! SETH    RWORK(8)  the element threshhold for sparsity determination
!                   when MOSS = 1 or 2.  If the absolute value of
!                   an estimated Jacobian element is .le. SETH, it
!                   will be assumed to be absent in the structure.
!                   The default value of SETH is 0.
!
! MAXORD  IWORK(5)  the maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!-----------------------------------------------------------------------
! Optional Outputs.
!
! As optional additional output from DLSODES, the variables listed
! below are quantities related to the performance of DLSODES
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of these outputs are defined
! on any successful return from DLSODES, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, outputs relevant to the error will be defined,
! as noted below.
!
! Name    Location      Meaning
!
! HU      RWORK(11) the step size in t last used (successfully).
!
! HCUR    RWORK(12) the step size to be attempted on the next step.
!
! TCUR    RWORK(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  On output, TCUR
!                   will always be at least as far as the argument
!                   T, but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) the number of steps taken for the problem so far.
!
! NFE     IWORK(12) the number of f evaluations for the problem so far,
!                   excluding those for structure determination
!                   (MOSS = 2).
!
! NJE     IWORK(13) the number of Jacobian evaluations for the problem
!                   so far, excluding those for structure determination
!                   (MOSS = 1).
!
! NQU     IWORK(14) the method order last used (successfully).
!
! NQCUR   IWORK(15) the order to be attempted on the next step.
!
! IMXER   IWORK(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( E(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) the length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) the length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NNZ     IWORK(19) the number of nonzero elements in the Jacobian
!                   matrix, including the diagonal (MITER = 1 or 2).
!                   (This may differ from that given by IA(NEQ+1)-1
!                   if MOSS = 0, because of added diagonal entries.)
!
! NGP     IWORK(20) the number of groups of column indices, used in
!                   difference quotient Jacobian aproximations if
!                   MITER = 2.  This is also the number of extra f
!                   evaluations needed for each Jacobian evaluation.
!
! NLU     IWORK(21) the number of sparse LU decompositions for the
!                   problem so far.
!
! LYH     IWORK(22) the base address in RWORK of the history array YH,
!                   described below in this list.
!
! IPIAN   IWORK(23) the base address of the structure descriptor array
!                   IAN, described below in this list.
!
! IPJAN   IWORK(24) the base address of the structure descriptor array
!                   JAN, described below in this list.
!
! NZL     IWORK(25) the number of nonzero elements in the strict lower
!                   triangle of the LU factorization used in the chord
!                   iteration (MITER = 1 or 2).
!
! NZU     IWORK(26) the number of nonzero elements in the strict upper
!                   triangle of the LU factorization used in the chord
!                   iteration (MITER = 1 or 2).
!                   The total number of nonzeros in the factorization
!                   is therefore NZL + NZU + NEQ.
!
! The following four arrays are segments of the RWORK array which
! may also be of interest to the user as optional outputs.
! For each array, the table below gives its internal name,
! its base address, and its description.
! For YH and ACOR, the base addresses are in RWORK (a real array).
! The integer arrays IAN and JAN are to be obtained by declaring an
! integer array IWK and identifying IWK(1) with RWORK(21), using either
! an equivalence statement or a subroutine call.  Then the base
! addresses IPIAN (of IAN) and IPJAN (of JAN) in IWK are to be obtained
! as optional outputs IWORK(23) and IWORK(24), respectively.
! Thus IAN(1) is IWK(IPIAN), etc.
!
! Name    Base Address      Description
!
! IAN    IPIAN (in IWK)  structure descriptor array of size NEQ + 1.
! JAN    IPJAN (in IWK)  structure descriptor array of size NNZ.
!         (see above)    IAN and JAN together describe the sparsity
!                        structure of the Jacobian matrix, as used by
!                        DLSODES when MITER = 1 or 2.
!                        JAN contains the row indices of the nonzero
!                        locations, reading in columnwise order, and
!                        IAN contains the starting locations in JAN of
!                        the descriptions of columns 1,...,NEQ, in
!                        that order, with IAN(1) = 1.  Thus for each
!                        j = 1,...,NEQ, the row indices i of the
!                        nonzero locations in column j are
!                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
!                        Note that IAN(NEQ+1) = NNZ + 1.
!                        (If MOSS = 0, IAN/JAN may differ from the
!                        input IA/JA because of a different ordering
!                        in each column, and added diagonal entries.)
!
! YH      LYH            the Nordsieck history array, of size NYH by
!          (optional     (NQCUR + 1), where NYH is the initial value
!           output)      of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the solution,
!                        evaluated at t = TCUR.  The base address LYH
!                        is another optional output, listed above.
!
! ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
!                        corrections on each step, scaled on output
!                        to represent the estimated local error in y
!                        on the last step.  This is the vector E  in
!                        the description of the error control.  It is
!                        defined only on a successful return from
!                        DLSODES.
!
!-----------------------------------------------------------------------
! Part 2.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DLSODES.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     Form of Call                  Function
!   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
!                             output of messages from DLSODES, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!   CALL XSETF(MFLAG)         Set a flag to control the printing of
!                             messages by DLSODES.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!   CALL DSRCMS(RSAV,ISAV,JOB) saves and restores the contents of
!                             the internal Common blocks used by
!                             DLSODES (see Part 3 below).
!                             RSAV must be a real array of length 224
!                             or more, and ISAV must be an integer
!                             array of length 71 or more.
!                             JOB=1 means save Common into RSAV/ISAV.
!                             JOB=2 means restore Common from RSAV/ISAV.
!                                DSRCMS is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with DLSODES.
!
!   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
!        (see below)          orders, at a specified point t, if
!                             desired.  It may be called only after
!                             a successful return from DLSODES.
!
! The detailed instructions for using DINTDY are as follows.
! The form of the call is:
!
!   LYH = IWORK(22)
!   CALL DINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = value of independent variable where answers are desired
!             (normally the same as the T last returned by DLSODES).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional outputs for TCUR and HU.)
! K         = integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (See optional outputs).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DLSODES directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with DINTDY.
! LYH       = the base address of the history array YH, obtained
!             as an optional output as shown above.
! NYH       = column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = a real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part 3.  Common Blocks.
!
! If DLSODES is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to DLSODES, and
!   (2) the two internal Common blocks
!         /DLS001/  of length  255  (218 double precision words
!                      followed by 37 integer words),
!         /DLSS01/  of length  40  (6 double precision words
!                      followed by 34 integer words),
!
! If DLSODES is used on a system in which the contents of internal
! Common blocks are not preserved between calls, the user should
! declare the above Common blocks in the calling program to insure
! that their contents are preserved.
!
! If the solution of a given problem by DLSODES is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more problems, the user should save,
! following the return from the last DLSODES call prior to the
! interruption, the contents of the call sequence variables and the
! internal Common blocks, and later restore these values before the
! next DLSODES call for that problem.  To save and restore the Common
! blocks, use Subroutine DSRCMS (see Part 2 above).
!
!-----------------------------------------------------------------------
! Part 4.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the DLSODES package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODES call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
! routine (see below), and also used by DLSODES in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! optional outputs.  In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of H**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in DEWSET the statements:
!     DOUBLE PRECISION RLS
!     COMMON /DLS001/ RLS(218),ILS(37)
!     NQ = ILS(33)
!     NST = ILS(34)
!     H = RLS(212)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) DVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM (N, V, W)
! where
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by Subroutine DEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of DVNORM suitable for use in the error control in DLSODES.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM routine might:
!   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of y.
!-----------------------------------------------------------------------
!
!***REVISION HISTORY  (YYYYMMDD)
! 19810120  DATE WRITTEN
! 19820315  Upgraded MDI in ODRV package: operates on M + M-transpose.
! 19820426  Numerous revisions in use of work arrays;
!           use wordlength ratio LENRAT; added IPISP & LRAT to Common;
!           added optional outputs IPIAN/IPJAN;
!           numerous corrections to comments.
! 19830503  Added routine CNTNZU; added NZL and NZU to /LSS001/;
!           changed ADJLR call logic; added optional outputs NZL & NZU;
!           revised counter initializations; revised PREP stmt. numbers;
!           corrections to comments throughout.
! 19870320  Corrected jump on test of umax in CDRV routine;
!           added ISTATE = -7 return.
! 19870330  Major update: corrected comments throughout;
!           removed TRET from Common; rewrote EWSET with 4 loops;
!           fixed t test in INTDY; added Cray directives in STODE;
!           in STODE, fixed DELP init. and logic around PJAC call;
!           combined routines to save/restore Common;
!           passed LEVEL = 0 in error message calls (except run abort).
! 20010425  Major update: convert source lines to upper case;
!           added *DECK lines; changed from 1 to * in dummy dimensions;
!           changed names R1MACH/D1MACH to RUMACH/DUMACH;
!           renamed routines for uniqueness across single/double prec.;
!           converted intrinsic names to generic form;
!           removed ILLIN and NTREP (data loaded) from Common;
!           removed all 'own' variables from Common;
!           changed error messages to quoted strings;
!           replaced XERRWV/XERRWD with 1993 revised version;
!           converted prologues, comments, error messages to mixed case;
!           converted arithmetic IF statements to logical IF statements;
!           numerous corrections to prologues and internal comments.
! 20010507  Converted single precision source to double precision.
! 20020502  Corrected declarations in descriptions of user routines.
! 20031105  Restored 'own' variables to Common blocks, to enable
!           interrupt/restart feature.
! 20031112  Added SAVE statements for data-loaded constants.
!
!-----------------------------------------------------------------------
! Other routines in the DLSODES package.
!
! In addition to Subroutine DLSODES, the DLSODES package includes the
! following subroutines and function routines:
!  DIPREP   acts as an iterface between DLSODES and DPREP, and also does
!           adjusting of work space pointers and work arrays.
!  DPREP    is called by DIPREP to compute sparsity and do sparse matrix
!           preprocessing if MITER = 1 or 2.
!  JGROUP   is called by DPREP to compute groups of Jacobian column
!           indices for use when MITER = 2.
!  ADJLR    adjusts the length of required sparse matrix work space.
!           It is called by DPREP.
!  CNTNZU   is called by DPREP and counts the nonzero elements in the
!           strict upper triangle of J + J-transpose, where J = df/dy.
!  DINTDY   computes an interpolated value of the y vector at t = TOUT.
!  DSTODE   is the core integrator, which does one step of the
!           integration and the associated error control.
!  DCFODE   sets all method coefficients and test constants.
!  DPRJS    computes and preprocesses the Jacobian matrix J = df/dy
!           and the Newton iteration matrix P = I - h*l0*J.
!  DSOLSS   manages solution of linear system in chord iteration.
!  DEWSET   sets the error weight vector EWT before each step.
!  DVNORM   computes the weighted RMS-norm of a vector.
!  DSRCMS   is a user-callable routine to save and restore
!           the contents of the internal Common blocks.
!  ODRV     constructs a reordering of the rows and columns of
!           a matrix by the minimum degree algorithm.  ODRV is a
!           driver routine which calls Subroutines MD, MDI, MDM,
!           MDP, MDU, and SRO.  See Ref. 2 for details.  (The ODRV
!           module has been modified since Ref. 2, however.)
!  CDRV     performs reordering, symbolic factorization, numerical
!           factorization, or linear system solution operations,
!           depending on a path argument ipath.  CDRV is a
!           driver routine which calls Subroutines NROC, NSFC,
!           NNFC, NNSC, and NNTC.  See Ref. 3 for details.
!           DLSODES uses CDRV to solve linear systems in which the
!           coefficient matrix is  P = I - con*J, where I is the
!           identity, con is a scalar, and J is an approximation to
!           the Jacobian df/dy.  Because CDRV deals with rowwise
!           sparsity descriptions, CDRV works with P-transpose, not P.
!  DUMACH   computes the unit roundoff in a machine-independent manner.
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
!           error messages and warnings.  XERRWD is machine-dependent.
! Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
! All the others are subroutines.
!
!-----------------------------------------------------------------------
INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT
INTEGER LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP
INTEGER IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, I1, I2, IFLAG, IMAX, IMUL, IMXER, IPFLAG, IPGO, IREM, J, KGO, LENRAT, LENYHT, LENIW, LENRW, LF0, LIA
INTEGER LJA, LRTEM, LWTEM, LYHD, LYHN, MF1, MORD, MXHNL0, MXSTP0, NCOLM
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
DIMENSION MORD(2)
LOGICAL IHIT
CHARACTER*60 MSG
SAVE LENRAT, MORD, MXSTP0, MXHNL0
!-----------------------------------------------------------------------
! The following two internal Common blocks contain
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine ("own" variables), and
! (b) variables which are communicated between subroutines.
! The block DLS001 is declared in subroutines DLSODES, DIPREP, DPREP,
! DINTDY, DSTODE, DPRJS, and DSOLSS.
! The block DLSS01 is declared in subroutines DLSODES, DIPREP, DPREP,
! DPRJS, and DSOLSS.
! Groups of variables are replaced by dummy arrays in the Common
! declarations in routines where those variables are not used.
!-----------------------------------------------------------------------
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, & 
IOWNS( 6), ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, & 
MSBP, MXNCF , N, NQ, NST, NFE, NJE, NQU
!
COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, &
IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, &
NLU, NNZ, NSP, NZL, NZU
!
DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
!-----------------------------------------------------------------------
! In the Data statement below, set LENRAT equal to the ratio of
! the wordlength for a real number to that for an integer.  Usually,
! LENRAT = 1 for single precision and 2 for double precision.  If the
! true ratio is not an integer, use the next smaller integer (.ge. 1).
!-----------------------------------------------------------------------
DATA LENRAT/2/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
IWORK(1:LWM) = 0
!-----------------------------------------------------------------------
IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GOTO 601
IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GOTO 602
IF (ISTATE .EQ. 1) GOTO 10
IF (INIT .EQ. 0) GOTO 603
IF (ISTATE .EQ. 2) GOTO 200
GOTO 20
10   INIT = 0
IF (TOUT .EQ. T) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all inputs and various initializations.
! If ISTATE = 1, the final setting of work space pointers, the matrix
! preprocessing, and other initializations are done in Block C.
!
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
20   IF (NEQ .LE. 0) GOTO 604
IF (ISTATE .EQ. 1) GOTO 25
IF (NEQ .GT. N) GOTO 605
25   N = NEQ
IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GOTO 606
IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GOTO 607
MOSS = MF/100
MF1 = MF - 100*MOSS
METH = MF1/10
MITER = MF1 - 10*METH
IF (MOSS .LT. 0 .OR. MOSS .GT. 2) GOTO 608
IF (METH .LT. 1 .OR. METH .GT. 2) GOTO 608
IF (MITER .LT. 0 .OR. MITER .GT. 3) GOTO 608
IF (MITER .EQ. 0 .OR. MITER .EQ. 3) MOSS = 0
! Next process and check the optional inputs. --------------------------
IF (IOPT .EQ. 1) GOTO 40
MAXORD = MORD(METH)
MXSTEP = MXSTP0
MXHNIL = MXHNL0
IF (ISTATE .EQ. 1) H0 = 0.0D0
HMXI = 0.0D0
HMIN = 0.0D0
SETH = 0.0D0
GOTO 60
40   MAXORD = IWORK(5)
IF (MAXORD .LT. 0) GOTO 611
IF (MAXORD .EQ. 0) MAXORD = 100
MAXORD = MIN(MAXORD,MORD(METH))
MXSTEP = IWORK(6)
IF (MXSTEP .LT. 0) GOTO 612
IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
MXHNIL = IWORK(7)
IF (MXHNIL .LT. 0) GOTO 613
IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
IF (ISTATE .NE. 1) GOTO 50
H0 = RWORK(5)
IF ((TOUT - T)*H0 .LT. 0.0D0) GOTO 614
50   HMAX = RWORK(6)
IF (HMAX .LT. 0.0D0) GOTO 615
HMXI = 0.0D0
IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
HMIN = RWORK(7)
IF (HMIN .LT. 0.0D0) GOTO 616
SETH = RWORK(8)
IF (SETH .LT. 0.0D0) GOTO 609
! Check RTOL and ATOL for legality. ------------------------------------
60   RTOLI = RTOL
ATOLI = ATOL(1)
DO 65 I = 1,N
IF (ITOL .GE. 3) RTOLI = RTOL
IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
IF (RTOLI .LT. 0.0D0) GOTO 619
IF (ATOLI .LT. 0.0D0) GOTO 620
65     CONTINUE
!-----------------------------------------------------------------------
! Compute required work array lengths, as far as possible, and test
! these against LRW and LIW.  Then set tentative pointers for work
! arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  WM, YH, SAVF, EWT, ACOR.
! If MITER = 1 or 2, the required length of the matrix work space WM
! is not yet known, and so a crude minimum value is used for the
! initial tests of LRW and LIW, and YH is temporarily stored as far
! to the right in RWORK as possible, to leave the maximum amount
! of space for WM for matrix preprocessing.  Thus if MITER = 1 or 2
! and MOSS .ne. 2, some of the segments of RWORK are temporarily
! omitted, as they are not needed in the preprocessing.  These
! omitted segments are: ACOR if ISTATE = 1, EWT and ACOR if ISTATE = 3
! and MOSS = 1, and SAVF, EWT, and ACOR if ISTATE = 3 and MOSS = 0.
!-----------------------------------------------------------------------
LRAT = LENRAT
IF (ISTATE .EQ. 1) NYH = N
LWMIN = 0
IF (MITER .EQ. 1) LWMIN = 4*N + 10*N/LRAT
IF (MITER .EQ. 2) LWMIN = 4*N + 11*N/LRAT
IF (MITER .EQ. 3) LWMIN = N + 2
LENYH = (MAXORD+1)*NYH
LREST = LENYH + 3*N
LENRW = 20 + LWMIN + LREST
IWORK(17) = LENRW
LENIW = 30
IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3) LENIW = LENIW + N + 1
IWORK(18) = LENIW
IF (LENRW .GT. LRW) GOTO 617
IF (LENIW .GT. LIW) GOTO 618
LIA = 31
IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3) LENIW = LENIW + IWORK(LIA+N) - 1
IWORK(18) = LENIW
IF (LENIW .GT. LIW) GOTO 618
LJA = LIA + N + 1
LIA = MIN(LIA,LIW)
LJA = MIN(LJA,LIW)
LWM = 21
IF (ISTATE .EQ. 1) NQ = 1
NCOLM = MIN(NQ+1,MAXORD+2)
LENYHM = NCOLM*NYH
LENYHT = LENYH
IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENYHT = LENYHM
IMUL = 2
IF (ISTATE .EQ. 3) IMUL = MOSS
IF (MOSS .EQ. 2) IMUL = 3
LRTEM = LENYHT + IMUL*N
LWTEM = LWMIN
IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LWTEM = LRW - 20 - LRTEM
LENWK = LWTEM
LYHN = LWM + LWTEM
LSAVF = LYHN + LENYHT
LEWT = LSAVF + N
LACOR = LEWT + N
ISTATC = ISTATE
IF (ISTATE .EQ. 1) GOTO 100
!-----------------------------------------------------------------------
! ISTATE = 3.  Move YH to its new location.
! Note that only the part of YH needed for the next step, namely
! MIN(NQ+1,MAXORD+2) columns, is actually moved.
! A temporary error weight array EWT is loaded if MOSS = 2.
! Sparse matrix processing is done in DIPREP/DPREP if MITER = 1 or 2.
! If MAXORD was reduced below NQ, then the pointers are finally set
! so that SAVF is identical to YH(*,MAXORD+2).
!-----------------------------------------------------------------------
LYHD = LYH - LYHN
IMAX = LYHN - 1 + LENYHM
! Move YH.  Move right if LYHD < 0; move left if LYHD > 0. -------------
IF (LYHD .LT. 0) THEN
DO 72 I = LYHN,IMAX
J = IMAX + LYHN - I
72       RWORK(J) = RWORK(J+LYHD)
ENDIF
IF (LYHD .GT. 0) THEN
DO 76 I = LYHN,IMAX
76       RWORK(I) = RWORK(I+LYHD)
ENDIF
80   LYH = LYHN
IWORK(22) = LYH
IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GOTO 92
IF (MOSS .NE. 2) GOTO 85
! Temporarily load EWT if MITER = 1 or 2 and MOSS = 2. -----------------
CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
DO 82 I = 1,N
IF (RWORK(I+LEWT-1) .LE. 0.0D0) GOTO 621
82     RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
85   CONTINUE
! DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
LSAVF = MIN(LSAVF,LRW)
LEWT = MIN(LEWT,LRW)
LACOR = MIN(LACOR,LRW)
CALL DIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC)
LENRW = LWM - 1 + LENWK + LREST
IWORK(17) = LENRW
IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
IPGO = -IPFLAG + 1
GOTO (90, 628, 629, 630, 631, 632, 633), IPGO
90   IWORK(22) = LYH
IF (LENRW .GT. LRW) GOTO 617
! Set flag to signal parameter changes to DSTODE. ----------------------
92   JSTART = -1
IF (N .EQ. NYH) GOTO 200
! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
I1 = LYH + L*NYH
I2 = LYH + (MAXORD + 1)*NYH - 1
IF (I1 .GT. I2) GOTO 200
DO 95 I = I1,I2
95     RWORK(I) = 0.0D0
GOTO 200
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! the sparse matrix preprocessing (MITER = 1 or 2), and the
! calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
100  CONTINUE
LYH = LYHN
IWORK(22) = LYH
TN = T
NST = 0
H = 1.0D0
NNZ = 0
NGP = 0
NZL = 0
NZU = 0
! Load the initial value vector in YH. ---------------------------------
DO 105 I = 1,N
105    RWORK(I+LYH-1) = Y(I)
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
LF0 = LYH + NYH
CALL F (NEQ, T, Y, RWORK(LF0))
NFE = 1
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
DO 110 I = 1,N
IF (RWORK(I+LEWT-1) .LE. 0.0D0) GOTO 621
110    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GOTO 120
! DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
LACOR = MIN(LACOR,LRW)
CALL DIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC)
LENRW = LWM - 1 + LENWK + LREST
IWORK(17) = LENRW
IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
IPGO = -IPFLAG + 1
GOTO (115, 628, 629, 630, 631, 632, 633), IPGO
115  IWORK(22) = LYH
IF (LENRW .GT. LRW) GOTO 617
! Check TCRIT for legality (ITASK = 4 or 5). ---------------------------
120  CONTINUE
IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GOTO 125
TCRIT = RWORK(1)
IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GOTO 625
IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0) H0 = TCRIT - T
! Initialize all remaining parameters. ---------------------------------
125  UROUND = DUMACH()
JSTART = 0
IF (MITER .NE. 0) RWORK(LWM) = SQRT(UROUND)
MSBJ = 50
NSLJ = 0
CCMXJ = 0.2D0
PSMALL = 1000.0D0*UROUND
RBIG = 0.01D0/PSMALL
NHNIL = 0
NJE = 0
NLU = 0
NSLAST = 0
HU = 0.0D0
NQU = 0
CCMAX = 0.3D0
MAXCOR = 3
MSBP = 20
MXNCF = 10
!-----------------------------------------------------------------------
! The coding below computes the step size, H0, to be attempted on the
! first step, unless the user has supplied a value for this.
! First check that TOUT - T differs significantly from zero.
! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
! so as to be between 100*UROUND and 1.0E-3.
! Then the computed value H0 is given by..
!                                      NEQ
!   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
!                                       1
! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
!         f(i)   = i-th component of initial value of f,
!         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
! The sign of H0 is inferred from the initial values of TOUT and T.
! ABS(H0) is made .le. ABS(TOUT-T) in any case.
!-----------------------------------------------------------------------
LF0 = LYH + NYH
IF (H0 .NE. 0.0D0) GOTO 180
TDIST = ABS(TOUT - T)
W0 = MAX(ABS(T),ABS(TOUT))
IF (TDIST .LT. 2.0D0*UROUND*W0) GOTO 622
TOL = RTOL
IF (ITOL .LE. 2) GOTO 140
DO 130 I = 1,N
130    TOL = MAX(TOL,RTOL)
140  IF (TOL .GT. 0.0D0) GOTO 160
ATOLI = ATOL(1)
DO 150 I = 1,N
IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
AYI = ABS(Y(I))
IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
150    CONTINUE
160  TOL = MAX(TOL,100.0D0*UROUND)
TOL = MIN(TOL,0.001D0)
SUM = DVNORM (N, RWORK(LF0), RWORK(LEWT))
SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
H0 = 1.0D0/SQRT(SUM)
H0 = MIN(H0,TDIST)
H0 = SIGN(H0,TOUT-T)
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180  RH = ABS(H0)*HMXI
IF (RH .GT. 1.0D0) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
H = H0
DO 190 I = 1,N
190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
GOTO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
200  NSLAST = NST
GOTO (210, 250, 220, 230, 240), ITASK
210  IF ((TN - TOUT)*H .LT. 0.0D0) GOTO 250
CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
IF (IFLAG .NE. 0) GOTO 627
T = TOUT
GOTO 420
220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
IF ((TP - TOUT)*H .GT. 0.0D0) GOTO 623
IF ((TN - TOUT)*H .LT. 0.0D0) GOTO 250
GOTO 400
230  TCRIT = RWORK(1)
IF ((TN - TCRIT)*H .GT. 0.0D0) GOTO 624
IF ((TCRIT - TOUT)*H .LT. 0.0D0) GOTO 625
IF ((TN - TOUT)*H .LT. 0.0D0) GOTO 245
CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
IF (IFLAG .NE. 0) GOTO 627
T = TOUT
GOTO 420
240  TCRIT = RWORK(1)
IF ((TN - TCRIT)*H .GT. 0.0D0) GOTO 624
245  HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
IF (IHIT) GOTO 400
TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GOTO 250
H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
IF (ISTATE .EQ. 2) JSTART = -2
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DSTODE.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
250  CONTINUE
IF ((NST-NSLAST) .GE. MXSTEP) GOTO 500
CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
DO 260 I = 1,N
IF (RWORK(I+LEWT-1) .LE. 0.0D0) GOTO 510
260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
IF (TOLSF .LE. 1.0D0) GOTO 280
TOLSF = TOLSF*2.0D0
IF (NST .EQ. 0) GOTO 626
GOTO 520
280  IF ((TN + H) .NE. TN) GOTO 290
NHNIL = NHNIL + 1
IF (NHNIL .GT. MXHNIL) GOTO 290
MSG = 'DLSODES- Warning..Internal T (=R1) and H (=R2) are'
CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='      such that in the machine, T + H = T on the next step  '
CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '     (H = step size). Solver will continue anyway.'
CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
IF (NHNIL .LT. MXHNIL) GOTO 290
MSG = 'DLSODES- Above warning has been issued I1 times.  '
CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '     It will not be issued again for this problem.'
CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
290  CONTINUE
!-----------------------------------------------------------------------
!    CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,WM,F,JAC,DPRJS,DSOLSS)
!-----------------------------------------------------------------------
CALL DSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT), RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), &
IWORK(LWM), F, JAC, DPRJS, DSOLSS)
KGO = 1 - KFLAG
GOTO (300, 530, 540, 550), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
300  INIT = 1
GOTO (310, 400, 330, 340, 350), ITASK
! ITASK = 1.  if TOUT has been reached, interpolate. -------------------
310  IF ((TN - TOUT)*H .LT. 0.0D0) GOTO 250
CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
T = TOUT
GOTO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330  IF ((TN - TOUT)*H .GE. 0.0D0) GOTO 400
GOTO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340  IF ((TN - TOUT)*H .LT. 0.0D0) GOTO 345
CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
T = TOUT
GOTO 420
345  HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
IF (IHIT) GOTO 400
TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GOTO 250
H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
JSTART = -2
GOTO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350  HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from DLSODES.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional outputs are loaded into the
! work arrays before returning.
!-----------------------------------------------------------------------
400  DO 410 I = 1,N
410    Y(I) = RWORK(I+LYH-1)
T = TN
IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GOTO 420
IF (IHIT) T = TCRIT
420  ISTATE = 2
RWORK(11) = HU
RWORK(12) = H
RWORK(13) = TN
IWORK(11) = NST
IWORK(12) = NFE
IWORK(13) = NJE
IWORK(14) = NQU
IWORK(15) = NQ
IWORK(19) = NNZ
IWORK(20) = NGP
IWORK(21) = NLU
IWORK(25) = NZL
IWORK(26) = NZU
RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! If there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH and T is set to TN.
! The optional outputs are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
500  MSG = 'DLSODES- At current T (=R1), MXSTEP (=I1) steps   '
CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      taken on this call before reaching TOUT     '
CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
ISTATE = -1
GOTO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510  EWTI = RWORK(LEWT+I-1)
MSG = 'DLSODES- At T (=R1), EWT(I1) has become R2 .le. 0.'
CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
ISTATE = -6
GOTO 580
! Too much accuracy requested for machine precision. -------------------
520  MSG = 'DLSODES- At T (=R1), too much accuracy requested  '
CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      for precision of machine..  See TOLSF (=R2) '
CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
RWORK(14) = TOLSF
ISTATE = -2
GOTO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530  MSG = 'DLSODES- At T(=R1) and step size H(=R2), the error'
CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      test failed repeatedly or with ABS(H) = HMIN'
CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
ISTATE = -4
GOTO 560
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540  MSG = 'DLSODES- At T (=R1) and step size H (=R2), the    '
CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      corrector convergence failed repeatedly     '
CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      or with ABS(H) = HMIN   '
CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
ISTATE = -5
GOTO 560
! KFLAG = -3.  Fatal error flag returned by DPRJS or DSOLSS (CDRV). ----
550  MSG = 'DLSODES- At T (=R1) and step size H (=R2), a fatal'
CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      error flag was returned by CDRV (by way of  '
CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG = '      Subroutine DPRJS or DSOLSS)       '
CALL XERRWD (MSG, 40, 207, 0, 0, 0, 0, 2, TN, H)
ISTATE = -7
GOTO 580
! Compute IMXER if relevant. -------------------------------------------
560  BIG = 0.0D0
IMXER = 1
DO 570 I = 1,N
SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
IF (BIG .GE. SIZE) GOTO 570
BIG = SIZE
IMXER = I
570    CONTINUE
IWORK(16) = IMXER
! Set Y vector, T, and optional outputs. -------------------------------
580  DO 590 I = 1,N
590    Y(I) = RWORK(I+LYH-1)
T = TN
RWORK(11) = HU
RWORK(12) = H
RWORK(13) = TN
IWORK(11) = NST
IWORK(12) = NFE
IWORK(13) = NJE
IWORK(14) = NQU
IWORK(15) = NQ
IWORK(19) = NNZ
IWORK(20) = NGP
IWORK(21) = NLU
IWORK(25) = NZL
IWORK(26) = NZU
RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.  If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
601  MSG = 'DLSODES- ISTATE (=I1) illegal.'
CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
IF (ISTATE .LT. 0) GOTO 800
GOTO 700
602  MSG = 'DLSODES- ITASK (=I1) illegal. '
CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
GOTO 700
603  MSG = 'DLSODES- ISTATE.gt.1 but DLSODES not initialized. '
CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
GOTO 700
604  MSG = 'DLSODES- NEQ (=I1) .lt. 1     '
CALL XERRWD (MSG, 30, 4, 0, 1, NEQ, 0, 0, 0.0D0, 0.0D0)
GOTO 700
605  MSG = 'DLSODES- ISTATE = 3 and NEQ increased (I1 to I2). '
CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ, 0, 0.0D0, 0.0D0)
GOTO 700
606  MSG = 'DLSODES- ITOL (=I1) illegal.  '
CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
GOTO 700
607  MSG = 'DLSODES- IOPT (=I1) illegal.  '
CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
GOTO 700
608  MSG = 'DLSODES- MF (=I1) illegal.    '
CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
GOTO 700
609  MSG = 'DLSODES- SETH (=R1) .lt. 0.0  '
CALL XERRWD (MSG, 30, 9, 0, 0, 0, 0, 1, SETH, 0.0D0)
GOTO 700
611  MSG = 'DLSODES- MAXORD (=I1) .lt. 0  '
CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
GOTO 700
612  MSG = 'DLSODES- MXSTEP (=I1) .lt. 0  '
CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
GOTO 700
613  MSG = 'DLSODES- MXHNIL (=I1) .lt. 0  '
CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
GOTO 700
614  MSG = 'DLSODES- TOUT (=R1) behind T (=R2)      '
CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
MSG = '      Integration direction is given by H0 (=R1)  '
CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
GOTO 700
615  MSG = 'DLSODES- HMAX (=R1) .lt. 0.0  '
CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
GOTO 700
616  MSG = 'DLSODES- HMIN (=R1) .lt. 0.0  '
CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
GOTO 700
617  MSG = 'DLSODES- RWORK length is insufficient to proceed. '
CALL XERRWD (MSG, 50, 17, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
GOTO 700
618  MSG = 'DLSODES- IWORK length is insufficient to proceed. '
CALL XERRWD (MSG, 50, 18, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENIW (=I1), exceeds LIW (=I2)'
CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
GOTO 700
619  MSG = 'DLSODES- RTOL(I1) is R1 .lt. 0.0        '
CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
GOTO 700
620  MSG = 'DLSODES- ATOL(I1) is R1 .lt. 0.0        '
CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
GOTO 700
621  EWTI = RWORK(LEWT+I-1)
MSG = 'DLSODES- EWT(I1) is R1 .le. 0.0         '
CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
GOTO 700
622  MSG='DLSODES- TOUT(=R1) too close to T(=R2) to start integration.'
CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
GOTO 700
623  MSG='DLSODES- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
GOTO 700
624  MSG='DLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
GOTO 700
625  MSG='DLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
GOTO 700
626  MSG = 'DLSODES- At start of problem, too much accuracy   '
CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='      requested for precision of machine..  See TOLSF (=R1) '
CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
RWORK(14) = TOLSF
GOTO 700
627  MSG = 'DLSODES- Trouble in DINTDY.  ITASK = I1, TOUT = R1'
CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
GOTO 700
628  MSG='DLSODES- RWORK length insufficient (for Subroutine DPREP).  '
CALL XERRWD (MSG, 60, 28, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 28, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
GOTO 700
629  MSG='DLSODES- RWORK length insufficient (for Subroutine JGROUP). '
CALL XERRWD (MSG, 60, 29, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 29, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
GOTO 700
630  MSG='DLSODES- RWORK length insufficient (for Subroutine ODRV).   '
CALL XERRWD (MSG, 60, 30, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 30, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
GOTO 700
631  MSG='DLSODES- Error from ODRV in Yale Sparse Matrix Package.     '
CALL XERRWD (MSG, 60, 31, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
IMUL = (IYS - 1)/N
IREM = IYS - IMUL*N
MSG='      At T (=R1), ODRV returned error flag = I1*NEQ + I2.   '
CALL XERRWD (MSG, 60, 31, 0, 2, IMUL, IREM, 1, TN, 0.0D0)
GOTO 700
632  MSG='DLSODES- RWORK length insufficient (for Subroutine CDRV).   '
CALL XERRWD (MSG, 60, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 32, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
GOTO 700
633  MSG='DLSODES- Error from CDRV in Yale Sparse Matrix Package.     '
CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
IMUL = (IYS - 1)/N
IREM = IYS - IMUL*N
MSG='      At T (=R1), CDRV returned error flag = I1*NEQ + I2.   '
CALL XERRWD (MSG, 60, 33, 0, 2, IMUL, IREM, 1, TN, 0.0D0)
IF (IMUL .EQ. 2) THEN
MSG='        Duplicate entry in sparsity structure descriptors.  '
CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
ENDIF
IF (IMUL .EQ. 3 .OR. IMUL .EQ. 6) THEN
MSG='        Insufficient storage for NSFC (called by CDRV).     '
CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
ENDIF
!
700  ISTATE = -3
RETURN
!
800  MSG = 'DLSODES- Run aborted.. apparent infinite loop.    '
CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
RETURN
!----------------------- End of Subroutine DLSODES ---------------------
END subroutine DLSODES

SUBROUTINE DSOLSS (WK, IWK, X, TEM)
INTEGER IWK
DOUBLE PRECISION WK, X, TEM
DIMENSION WK(*), IWK(*), X(*), TEM(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I
DOUBLE PRECISION DI, HL0, PHL0, R
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls CDRV to accomplish this.
! If MITER = 3 it updates the coefficient H*EL0 in the diagonal
! matrix, and then computes the solution.
! communication with DSOLSS uses the following variables:
! WK    = real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WK(3).
!         WK also contains the following matrix-related data:
!         WK(1) = SQRT(UROUND) (not used here),
!         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
! IWK   = integer work space for matrix-related data, assumed to
!         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!         are assumed to have identical locations.
! X     = the right-hand side vector on input, and the solution vector
!         on output, of length N.
! TEM   = vector of work space of length N, not used in this version.
! IERSL = output flag (in Common).
!         IERSL = 0  if no trouble occurred.
!         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
!                    This should never occur and is considered fatal.
!         IERSL = 1  if a singular matrix arose with MITER = 3.
! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
IERSL = 0
GOTO (100, 100, 300), MITER
100  CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),WK(IPA),X,X,NSP,IWK(IPISP),WK(IPRSP),IESP,4,IERSL)
IF (IERSL .NE. 0) IERSL = -1
RETURN
!
300  PHL0 = WK(2)
HL0 = H*EL0
WK(2) = HL0
IF (HL0 .EQ. PHL0) GOTO 330
R = HL0/PHL0
DO 320 I = 1,N
DI = 1.0D0 - R*(1.0D0 - 1.0D0/WK(I+2))
IF (ABS(DI) .EQ. 0.0D0) GOTO 390
320    WK(I+2) = 1.0D0/DI
330  DO 340 I = 1,N
340    X(I) = WK(I+2)*X(I)
RETURN
390  IERSL = 1
RETURN
!
!----------------------- End of Subroutine DSOLSS ----------------------
END subroutine DSOLSS

subroutine cdrv(n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
!*** subroutine cdrv
!*** driver for subroutines for solving sparse nonsymmetric systems of
!       linear equations (compressed pointer storage)
!
!
!    parameters
!    class abbreviations are--
!       n - integer variable
!       f - real variable
!       v - supplies a value to the driver
!       r - returns a result from the driver
!       i - used internally by the driver
!       a - array
!
! class - parameter
! ------+----------
!       -
!         the nonzero entries of the coefficient matrix m are stored
!    row-by-row in the array a.  to identify the individual nonzero
!    entries in each row, we need to know in which column each entry
!    lies.  the column indices which correspond to the nonzero entries
!    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  in addition, we need to know where each row starts and
!    how long it is.  the index positions in ja and a where the rows of
!    m begin are stored in the array ia.  i.e., if m(i,j) is the first
!    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
!    ia(i) = k.  moreover, the index in ja and a of the first location
!    following the last element in the last row is stored in ia(n+1).
!    thus, the number of entries in the i-th row is given by
!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!    consecutively in
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!    and the corresponding column indices are stored consecutively in
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!    for example, the 5 by 5 matrix
!                ( 1. 0. 2. 0. 0.)
!                ( 0. 3. 0. 0. 0.)
!            m = ( 0. 4. 5. 6. 0.)
!                ( 0. 0. 0. 7. 0.)
!                ( 0. 0. 0. 8. 9.)
!    would be stored as
!               - 1  2  3  4  5  6  7  8  9
!            ---+--------------------------
!            ia - 1  3  4  7  8 10
!            ja - 1  3  2  2  3  4  4  4  5
!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
! nv    - n     - number of variables/equations.
! fva   - a     - nonzero entries of the coefficient matrix m, stored
!       -           by rows.
!       -           size = number of nonzero entries in m.
! nva   - ia    - pointers to delimit the rows in a.
!       -           size = n+1.
! nva   - ja    - column numbers corresponding to the elements of a.
!       -           size = size of a.
! fva   - b     - right-hand side b.  b and z can the same array.
!       -           size = n.
! fra   - z     - solution x.  b and z can be the same array.
!       -           size = n.
!
!         the rows and columns of the original matrix m can be
!    reordered (e.g., to reduce fillin or ensure numerical stability)
!    before calling the driver.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
!    in the original order.
!         if the columns have been reordered (i.e.,  c(i).ne.i  for some
!    i), then the driver will call a subroutine (nroc) which rearranges
!    each row of ja and a, leaving the rows in the original order, but
!    placing the elements of each row in increasing order with respect
!    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
!    been called already.
!
! nva   - r     - ordering of the rows of m.
!       -           size = n.
! nva   - c     - ordering of the columns of m.
!       -           size = n.
! nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
!       -           ic(c(i)) = i  for i=1,...,n.
!       -           size = n.
!
!         the solution of the system of linear equations is divided into
!    three stages --
!      nsfc -- the matrix m is processed symbolically to determine where
!               fillin will occur during the numeric factorization.
!      nnfc -- the matrix m is factored numerically into the product ldu
!               of a unit lower triangular matrix l, a diagonal matrix
!               d, and a unit upper triangular matrix u, and the system
!               mx = b  is solved.
!      nnsc -- the linear system  mx = b  is solved using the ldu
!  or           factorization from nnfc.
!      nntc -- the transposed linear system  mt x = b  is solved using
!               the ldu factorization from nnf.
!    for several systems whose coefficient matrices have the same
!    nonzero structure, nsfc need be done only once (for the first
!    system).  then nnfc is done once for each additional system.  for
!    several systems with the same coefficient matrix, nsfc and nnfc
!    need be done only once (for the first system).  then nnsc or nntc
!    is done once for each additional right-hand side.
!
! nv    - path  - path specification.  values and their meanings are --
!       -           1  perform nroc, nsfc, and nnfc.
!       -           2  perform nnfc only  (nsfc is assumed to have been
!       -               done in a manner compatible with the storage
!       -               allocation used in the driver).
!       -           3  perform nnsc only  (nsfc and nnfc are assumed to
!       -               have been done in a manner compatible with the
!       -               storage allocation used in the driver).
!       -           4  perform nntc only  (nsfc and nnfc are assumed to
!       -               have been done in a manner compatible with the
!       -               storage allocation used in the driver).
!       -           5  perform nroc and nsfc.
!
!         various errors are detected by the driver and the individual
!    subroutines.
!
! nr    - flag  - error flag.  values and their meanings are --
!       -             0     no errors detected
!       -             n+k   null row in a  --  row = k
!       -            2n+k   duplicate entry in a  --  row = k
!       -            3n+k   insufficient storage in nsfc  --  row = k
!       -            4n+1   insufficient storage in nnfc
!       -            5n+k   null pivot  --  row = k
!       -            6n+k   insufficient storage in nsfc  --  row = k
!       -            7n+1   insufficient storage in nnfc
!       -            8n+k   zero pivot  --  row = k
!       -           10n+1   insufficient storage in cdrv
!       -           11n+1   illegal path specification
!
!         working storage is needed for the factored form of the matrix
!    m plus various temporary vectors.  the arrays isp and rsp should be
!    equivalenced.  integer storage is allocated from the beginning of
!    isp and real storage from the end of rsp.
!
! nv    - nsp   - declared dimension of rsp.  nsp generally must
!       -           be larger than  8n+2 + 2k  (where  k = (number of
!       -           nonzero entries in m)).
! nvira - isp   - integer working storage divided up into various arrays
!       -           needed by the subroutines.  isp and rsp should be
!       -           equivalenced.
!       -           size = lratio*nsp.
! fvira - rsp   - real working storage divided up into various arrays
!       -           needed by the subroutines.  isp and rsp should be
!       -           equivalenced.
!       -           size = nsp.
! nr    - esp   - if sufficient storage was available to perform the
!       -           symbolic factorization (nsfc), then esp is set to
!       -           the amount of excess storage provided (negative if
!       -           insufficient storage was available to perform the
!       -           numeric factorization (nnfc)).
!
!
!  conversion to double precision
!
!    to convert these routines for double precision arrays..
!    (1) use the double precision declarations in place of the real
!    declarations in each subprogram, as given in comment cards.
!    (2) change the data-loaded value of the integer  lratio
!    in subroutine cdrv, as indicated below.
!    (3) change e0 to d0 in the constants in statement number 10
!    in subroutine nnfc and the line following that.
!
integer  r(*), c(*), ic(*),  ia(*), ja(*),  isp(*), esp,  path,flag,  d, u, q, row, tmp, ar,  umax
!     real  a(*), b(*), z(*), rsp(*)
double precision  a(*), b(*), z(*), rsp(*)
!
!  set lratio equal to the ratio between the length of floating point
!  and integer array data.  e. g., lratio = 1 for (real, integer),
!  lratio = 2 for (double precision, integer)
!
data lratio/2/
!
if (path.lt.1 .or. 5.lt.path)  GOTO 111
!******initialize and divide up temporary storage  *******************
il   = 1
ijl  = il  + (n+1)
iu   = ijl +   n
iju  = iu  + (n+1)
irl  = iju +   n
jrl  = irl +   n
jl   = jrl +   n
!
!  ******  reorder a if necessary, call nsfc if flag is set  ***********
if ((path-1) * (path-5) .ne. 0)  GOTO 5
max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
jlmax = max/2
q     = jl   + jlmax
ira   = q    + (n+1)
jra   = ira  +   n
irac  = jra  +   n
iru   = irac +   n
jru   = iru  +   n
jutmp = jru  +   n
jumax = lratio*nsp  + 1 - jutmp
esp = max/lratio
if (jlmax.le.0 .or. jumax.le.0)  GOTO 110
!
do 1 i=1,n
if (c(i).ne.i)  GOTO 2
1      continue
GOTO 3
2    ar = nsp + 1 - n
call  nroc(n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
if (flag.ne.0)  GOTO 100
!
3    call  nsfc(n, r, ic, ia,ja,jlmax, isp(il), isp(jl), isp(ijl),jumax, isp(iu), isp(jutmp), isp(iju),isp(q), isp(ira), &
isp(jra), isp(irac),isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
if(flag .ne. 0)  GOTO 100
!  ******  move ju next to jl  *****************************************
jlmax = isp(ijl+n-1)
ju    = jl + jlmax
jumax = isp(iju+n-1)
if (jumax.le.0)  GOTO 5
do 4 j=1,jumax
4      isp(ju+j-1) = isp(jutmp+j-1)
!
!  ******  call remaining subroutines  *********************************
5  jlmax = isp(ijl+n-1)
ju    = jl  + jlmax
jumax = isp(iju+n-1)
l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
lmax  = isp(il+n) - 1
d     = l   + lmax
u     = d   + n
row   = nsp + 1 - n
tmp   = row - n
umax  = tmp - u
esp   = umax - (isp(iu+n) - 1)
!
if ((path-1) * (path-2) .ne. 0)  GOTO 6
if (umax.lt.0)  GOTO 110
call nnfc(n,  r, c, ic,  ia, ja, a, z, b,lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),umax, isp(iu), &
isp(ju), isp(iju), rsp(u),rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
if(flag .ne. 0)  GOTO 100
!
6  if ((path-3) .ne. 0)  GOTO 7
call nnsc(n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),z, b,  rsp(tmp))
!
7  if ((path-4) .ne. 0)  GOTO 8
call nntc(n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),z, b,  rsp(tmp))
8  return
!
! ** error.. error detected in nroc, nsfc, nnfc, or nnsc
100  return
! ** error.. insufficient storage
110  flag = 10*n + 1
return
! ** error.. illegal path specification
111  flag = 11*n + 1
return
end subroutine cdrv

subroutine nnsc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
!*** subroutine nnsc
!*** numerical solution of sparse nonsymmetric system of linear
!      equations given ldu-factorization (compressed pointer storage)
!
!
!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!       output variables.. z
!
!       parameters used internally..
! fia   - tmp   - temporary vector which gets result of solving  ly = b.
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row of
!      u or l  to be used.
!
integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
!     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
double precision  l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
!
!  ******  set tmp to reordered b  *************************************
do 1 k=1,n
1    tmp(k) = b(r(k))
!  ******  solve  ly = b  by forward substitution  *********************
do 3 k=1,n
jmin = il(k)
jmax = il(k+1) - 1
tmpk = -d(k) * tmp(k)
tmp(k) = -tmpk
if (jmin .gt. jmax) GOTO 3
ml = ijl(k) - jmin
do 2 j=jmin,jmax
2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
3    continue
!  ******  solve  ux = y  by back substitution  ************************
k = n
do 6 i=1,n
sum = -tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax) GOTO 5
mu = iju(k) - jmin
do 4 j=jmin,jmax
4      sum = sum + u(j) * tmp(ju(mu+j))
5    tmp(k) = -sum
z(c(k)) = -sum
k = k - 1
6    continue
return
end subroutine nnsc

subroutine nntc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
!*** subroutine nntc
!*** numeric solution of the transpose of a sparse nonsymmetric system
!      of linear equations given lu-factorization (compressed pointer
!      storage)
!
!
!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!       output variables.. z
!
!       parameters used internally..
! fia   - tmp   - temporary vector which gets result of solving ut y = b
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row of
!      u or l  to be used.
!
integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
!     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
double precision l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
!
!  ******  set tmp to reordered b  *************************************
do 1 k=1,n
1    tmp(k) = b(c(k))
!  ******  solve  ut y = b  by forward substitution  *******************
do 3 k=1,n
jmin = iu(k)
jmax = iu(k+1) - 1
tmpk = -tmp(k)
if (jmin .gt. jmax) GOTO 3
mu = iju(k) - jmin
do 2 j=jmin,jmax
2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
3    continue
!  ******  solve  lt x = y  by back substitution  **********************
k = n
do 6 i=1,n
sum = -tmp(k)
jmin = il(k)
jmax = il(k+1) - 1
if (jmin .gt. jmax) GOTO 5
ml = ijl(k) - jmin
do 4 j=jmin,jmax
4      sum = sum + l(j) * tmp(jl(ml+j))
5    tmp(k) = -sum * d(k)
z(r(k)) = tmp(k)
k = k - 1
6    continue
return
end subroutine nntc

subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
!
!       ----------------------------------------------------------------
!
!               yale sparse matrix package - nonsymmetric codes
!                    solving the system of equations mx = b
!
!    i.   calling sequences
!         the coefficient matrix can be processed by an ordering routine
!    (e.g., to reduce fillin or ensure numerical stability) before using
!    the remaining subroutines.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
!    is used, then nroc should be used to reorder the coefficient matrix
!    the calling sequence is --
!        (       (matrix ordering))
!        (nroc   (matrix reordering))
!         nsfc   (symbolic factorization to determine where fillin will
!                  occur during numeric factorization)
!         nnfc   (numeric factorization into product ldu of unit lower
!                  triangular matrix l, diagonal matrix d, and unit
!                  upper triangular matrix u, and solution of linear
!                  system)
!         nnsc   (solution of linear system for additional right-hand
!                  side using ldu factorization from nnfc)
!    (if only one system of equations is to be solved, then the
!    subroutine trk should be used.)
!
!    ii.  storage of sparse matrices
!         the nonzero entries of the coefficient matrix m are stored
!    row-by-row in the array a.  to identify the individual nonzero
!    entries in each row, we need to know in which column each entry
!    lies.  the column indices which correspond to the nonzero entries
!    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  in addition, we need to know where each row starts and
!    how long it is.  the index positions in ja and a where the rows of
!    m begin are stored in the array ia.  i.e., if m(i,j) is the first
!    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
!    ia(i) = k.  moreover, the index in ja and a of the first location
!    following the last element in the last row is stored in ia(n+1).
!    thus, the number of entries in the i-th row is given by
!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!    consecutively in
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!    and the corresponding column indices are stored consecutively in
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!    for example, the 5 by 5 matrix
!                ( 1. 0. 2. 0. 0.)
!                ( 0. 3. 0. 0. 0.)
!            m = ( 0. 4. 5. 6. 0.)
!                ( 0. 0. 0. 7. 0.)
!                ( 0. 0. 0. 8. 9.)
!    would be stored as
!               - 1  2  3  4  5  6  7  8  9
!            ---+--------------------------
!            ia - 1  3  4  7  8 10
!            ja - 1  3  2  2  3  4  4  4  5
!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
!         the strict upper (lower) triangular portion of the matrix
!    u (l) is stored in a similar fashion using the arrays  iu, ju, u
!    (il, jl, l)  except that an additional array iju (ijl) is used to
!    compress storage of ju (jl) by allowing some sequences of column
!    (row) indices to used for more than one row (column)  (n.b., l is
!    stored by columns).  iju(k) (ijl(k)) points to the starting
!    location in ju (jl) of entries for the kth row (column).
!    compression in ju (jl) occurs in two ways.  first, if a row
!    (column) i was merged into the current row (column) k, and the
!    number of elements merged in from (the tail portion of) row
!    (column) i is the same as the final length of row (column) k, then
!    the kth row (column) and the tail of row (column) i are identical
!    and iju(k) (ijl(k)) points to the start of the tail.  second, if
!    some tail portion of the (k-1)st row (column) is identical to the
!    head of the kth row (column), then iju(k) (ijl(k)) points to the
!    start of that tail portion.  for example, the nonzero structure of
!    the strict upper triangular part of the matrix
!            d 0 x x x
!            0 d 0 x x
!            0 0 d x 0
!            0 0 0 d x
!            0 0 0 0 d
!    would be represented as
!                - 1 2 3 4 5 6
!            ----+------------
!             iu - 1 4 6 7 8 8
!             ju - 3 4 5 4
!            iju - 1 2 4 3           .
!    the diagonal entries of l and u are assumed to be equal to one and
!    are not stored.  the array d contains the reciprocals of the
!    diagonal entries of the matrix d.
!
!    iii. additional storage savings
!         in nsfc, r and ic can be the same array in the calling
!    sequence if no reordering of the coefficient matrix has been done.
!         in nnfc, r, c, and ic can all be the same array if no
!    reordering has been done.  if only the rows have been reordered,
!    then c and ic can be the same array.  if the row and column
!    orderings are the same, then r and c can be the same array.  z and
!    row can be the same array.
!         in nnsc or nntc, r and c can be the same array if no
!    reordering has been done or if the row and column orderings are the
!    same.  z and b can be the same array.  however, then b will be
!    destroyed.
!
!    iv.  parameters
!         following is a list of parameters to the programs.  names are
!    uniform among the various subroutines.  class abbreviations are --
!       n - integer variable
!       f - real variable
!       v - supplies a value to a subroutine
!       r - returns a result from a subroutine
!       i - used internally by a subroutine
!       a - array
!
! class - parameter
! ------+----------
! fva   - a     - nonzero entries of the coefficient matrix m, stored
!       -           by rows.
!       -           size = number of nonzero entries in m.
! fva   - b     - right-hand side b.
!       -           size = n.
! nva   - c     - ordering of the columns of m.
!       -           size = n.
! fvra  - d     - reciprocals of the diagonal entries of the matrix d.
!       -           size = n.
! nr    - flag  - error flag.  values and their meanings are --
!       -            0     no errors detected
!       -            n+k   null row in a  --  row = k
!       -           2n+k   duplicate entry in a  --  row = k
!       -           3n+k   insufficient storage for jl  --  row = k
!       -           4n+1   insufficient storage for l
!       -           5n+k   null pivot  --  row = k
!       -           6n+k   insufficient storage for ju  --  row = k
!       -           7n+1   insufficient storage for u
!       -           8n+k   zero pivot  --  row = k
! nva   - ia    - pointers to delimit the rows of a.
!       -           size = n+1.
! nvra  - ijl   - pointers to the first element in each column in jl,
!       -           used to compress storage in jl.
!       -           size = n.
! nvra  - iju   - pointers to the first element in each row in ju, used
!       -           to compress storage in ju.
!       -           size = n.
! nvra  - il    - pointers to delimit the columns of l.
!       -           size = n+1.
! nvra  - iu    - pointers to delimit the rows of u.
!       -           size = n+1.
! nva   - ja    - column numbers corresponding to the elements of a.
!       -           size = size of a.
! nvra  - jl    - row numbers corresponding to the elements of l.
!       -           size = jlmax.
! nv    - jlmax - declared dimension of jl.  jlmax must be larger than
!       -           the number of nonzeros in the strict lower triangle
!       -           of m plus fillin minus compression.
! nvra  - ju    - column numbers corresponding to the elements of u.
!       -           size = jumax.
! nv    - jumax - declared dimension of ju.  jumax must be larger than
!       -           the number of nonzeros in the strict upper triangle
!       -           of m plus fillin minus compression.
! fvra  - l     - nonzero entries in the strict lower triangular portion
!       -           of the matrix l, stored by columns.
!       -           size = lmax.
! nv    - lmax  - declared dimension of l.  lmax must be larger than
!       -           the number of nonzeros in the strict lower triangle
!       -           of m plus fillin  (il(n+1)-1 after nsfc).
! nv    - n     - number of variables/equations.
! nva   - r     - ordering of the rows of m.
!       -           size = n.
! fvra  - u     - nonzero entries in the strict upper triangular portion
!       -           of the matrix u, stored by rows.
!       -           size = umax.
! nv    - umax  - declared dimension of u.  umax must be larger than
!       -           the number of nonzeros in the strict upper triangle
!       -           of m plus fillin  (iu(n+1)-1 after nsfc).
! fra   - z     - solution x.
!       -           size = n.
!
!       ----------------------------------------------------------------
!
!*** subroutine nroc
!*** reorders rows of a, leaving row order unchanged
!
!
!       input parameters.. n, ic, ia, ja, a
!       output parameters.. ja, a, flag
!
!       parameters used internally..
! nia   - p     - at the kth step, p is a linked list of the reordered
!       -           column indices of the kth row of a.  p(n+1) points
!       -           to the first entry in the list.
!       -           size = n+1.
! nia   - jar   - at the kth step,jar contains the elements of the
!       -           reordered column indices of a.
!       -           size = n.
! fia   - ar    - at the kth step, ar contains the elements of the
!       -           reordered row of a.
!       -           size = n.
!
integer  ic(*), ia(*), ja(*), jar(*), p(*), flag
!     real  a(*), ar(*)
double precision  a(*), ar(*)
!
!  ******  for each nonempty row  *******************************
do 5 k=1,n
jmin = ia(k)
jmax = ia(k+1) - 1
if(jmin .gt. jmax) GOTO 5
p(n+1) = n + 1
!  ******  insert each element in the list  *********************
do 3 j=jmin,jmax
newj = ic(ja(j))
i = n + 1
1      if(p(i) .ge. newj) GOTO 2
i = p(i)
GOTO 1
2      if(p(i) .eq. newj) GOTO 102
p(newj) = p(i)
p(i) = newj
jar(newj) = ja(j)
ar(newj) = a(j)
3      continue
!  ******  replace old row in ja and a  *************************
i = n + 1
do 4 j=jmin,jmax
i = p(i)
ja(j) = jar(i)
4      a(j) = ar(i)
5    continue
flag = 0
return
!
! ** error.. duplicate entry in a
102  flag = n + k
return
end subroutine nroc

subroutine nsfc(n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,q, ira,jra, irac, irl,jrl, iru,jru, flag)
!*** subroutine nsfc
!*** symbolic ldu-factorization of nonsymmetric sparse matrix
!      (compressed pointer storage)
!
!
!       input variables.. n, r, ic, ia, ja, jlmax, jumax.
!       output variables.. il, jl, ijl, iu, ju, iju, flag.
!
!       parameters used internally..
! nia   - q     - suppose  m*  is the result of reordering  m.  if
!       -           processing of the ith row of  m*  (hence the ith
!       -           row of  u) is being done,  q(j)  is initially
!       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
!       -           values need not be stored, each entry points to the
!       -           next nonzero and  q(n+1)  points to the first.  n+1
!       -           indicates the end of the list.  for example, if n=9
!       -           and the 5th row of  m*  is
!       -              0 x x 0 x 0 0 x 0
!       -           then  q  will initially be
!       -              a a a a 8 a a 10 5           (a - arbitrary).
!       -           as the algorithm proceeds, other elements of  q
!       -           are inserted in the list because of fillin.
!       -           q  is used in an analogous manner to compute the
!       -           ith column of  l.
!       -           size = n+1.
! nia   - ira,  - vectors used to find the columns of  m.  at the kth
! nia   - jra,      step of the factorization,  irac(k)  points to the
! nia   - irac      head of a linked list in  jra  of row indices i
!       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
!       -           indicates the end of the list.  ira(i)  (i.ge.k)
!       -           points to the smallest j such that j .ge. k and
!       -           m(i,j)  is nonzero.
!       -           size of each = n.
! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
! nia   - jrl       of the factorization,  jrl(k)  points to the head
!       -           of a linked list in  jrl  of column indices j
!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!       -           points to the smallest i such that i .ge. k and
!       -           l(i,j)  is nonzero.
!       -           size of each = n.
! nia   - iru,  - vectors used in a manner analogous to  irl and jrl
! nia   - jru       to find the columns of  u.
!       -           size of each = n.
!
!  internal variables..
!    jlptr - points to the last position used in  jl.
!    juptr - points to the last position used in  ju.
!    jmin,jmax - are the indices in  a or u  of the first and last
!                elements to be examined in a given row.
!                for example,  jmin=ia(k), jmax=ia(k+1)-1.
!
integer cend, qm, rend, rk, vj
integer ia(*), ja(*), ira(*), jra(*), il(*), jl(*), ijl(*)
integer iu(*), ju(*), iju(*), irl(*), jrl(*), iru(*), jru(*)
integer r(*), ic(*), q(*), irac(*), flag
!
!  ******  initialize pointers  ****************************************
np1 = n + 1
jlmin = 1
jlptr = 0
il(1) = 1
jumin = 1
juptr = 0
iu(1) = 1
do 1 k=1,n
irac(k) = 0
jra(k) = 0
jrl(k) = 0
1    jru(k) = 0
!  ******  initialize column pointers for a  ***************************
do 2 k=1,n
rk = r(k)
iak = ia(rk)
if (iak .ge. ia(rk+1))  GOTO 101
jaiak = ic(ja(iak))
if (jaiak .gt. k)  GOTO 105
jra(k) = irac(jaiak)
irac(jaiak) = k
2    ira(k) = iak
!
!  ******  for each column of l and row of u  **************************
do 41 k=1,n
!
!  ******  initialize q for computing kth column of l  *****************
q(np1) = np1
luk = -1
!  ******  by filling in kth column of a  ******************************
vj = irac(k)
if (vj .eq. 0)  GOTO 5
3      qm = np1
4      m = qm
qm =  q(m)
if (qm .lt. vj)  GOTO 4
if (qm .eq. vj)  GOTO 102
luk = luk + 1
q(m) = vj
q(vj) = qm
vj = jra(vj)
if (vj .ne. 0)  GOTO 3
!  ******  link through jru  *******************************************
5    lastid = 0
lasti = 0
ijl(k) = jlptr
i = k
6      i = jru(i)
if (i .eq. 0)  GOTO 10
qm = np1
jmin = irl(i)
jmax = ijl(i) + il(i+1) - il(i) - 1
long = jmax - jmin
if (long .lt. 0)  GOTO 6
jtmp = jl(jmin)
if (jtmp .ne. k)  long = long + 1
if (jtmp .eq. k)  r(i) = -r(i)
if (lastid .ge. long)  GOTO 7
lasti = i
lastid = long
!  ******  and merge the corresponding columns into the kth column  ****
7      do 9 j=jmin,jmax
vj = jl(j)
8        m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 8
if (qm .eq. vj)  GOTO 9
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
9        continue
GOTO 6
!  ******  lasti is the longest column merged into the kth  ************
!  ******  see if it equals the entire kth column  *********************
10    qm = q(np1)
if (qm .ne. k)  GOTO 105
if (luk .eq. 0)  GOTO 17
if (lastid .ne. luk)  GOTO 11
!  ******  if so, jl can be compressed  ********************************
irll = irl(lasti)
ijl(k) = irll + 1
if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
GOTO 17
!  ******  if not, see if kth column can overlap the previous one  *****
11    if (jlmin .gt. jlptr)  GOTO 15
qm = q(qm)
do 12 j=jlmin,jlptr
if (jl(j) - qm)  12, 13, 15
12      continue
GOTO 15
13    ijl(k) = j
do 14 i=j,jlptr
if (jl(i) .ne. qm)  GOTO 15
qm = q(qm)
if (qm .gt. n)  GOTO 17
14      continue
jlptr = j - 1
!  ******  move column indices from q to jl, update vectors  ***********
15    jlmin = jlptr + 1
ijl(k) = jlmin
if (luk .eq. 0)  GOTO 17
jlptr = jlptr + luk
if (jlptr .gt. jlmax)  GOTO 103
qm = q(np1)
do 16 j=jlmin,jlptr
qm = q(qm)
16        jl(j) = qm
17    irl(k) = ijl(k)
il(k+1) = il(k) + luk
!
!  ******  initialize q for computing kth row of u  ********************
q(np1) = np1
luk = -1
!  ******  by filling in kth row of reordered a  ***********************
rk = r(k)
jmin = ira(k)
jmax = ia(rk+1) - 1
if (jmin .gt. jmax)  GOTO 20
do 19 j=jmin,jmax
vj = ic(ja(j))
qm = np1
18      m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 18
if (qm .eq. vj)  GOTO 102
luk = luk + 1
q(m) = vj
q(vj) = qm
19      continue
!  ******  link through jrl,  ******************************************
20    lastid = 0
lasti = 0
iju(k) = juptr
i = k
i1 = jrl(k)
21      i = i1
if (i .eq. 0)  GOTO 26
i1 = jrl(i)
qm = np1
jmin = iru(i)
jmax = iju(i) + iu(i+1) - iu(i) - 1
long = jmax - jmin
if (long .lt. 0)  GOTO 21
jtmp = ju(jmin)
if (jtmp .eq. k)  GOTO 22
!  ******  update irl and jrl, *****************************************
long = long + 1
cend = ijl(i) + il(i+1) - il(i)
irl(i) = irl(i) + 1
if (irl(i) .ge. cend)  GOTO 22
j = jl(irl(i))
jrl(i) = jrl(j)
jrl(j) = i
22      if (lastid .ge. long)  GOTO 23
lasti = i
lastid = long
!  ******  and merge the corresponding rows into the kth row  **********
23      do 25 j=jmin,jmax
vj = ju(j)
24        m = qm
qm = q(m)
if (qm .lt. vj)  GOTO 24
if (qm .eq. vj)  GOTO 25
luk = luk + 1
q(m) = vj
q(vj) = qm
qm = vj
25        continue
GOTO 21
!  ******  update jrl(k) and irl(k)  ***********************************
26    if (il(k+1) .le. il(k))  GOTO 27
j = jl(irl(k))
jrl(k) = jrl(j)
jrl(j) = k
!  ******  lasti is the longest row merged into the kth  ***************
!  ******  see if it equals the entire kth row  ************************
27    qm = q(np1)
if (qm .ne. k)  GOTO 105
if (luk .eq. 0)  GOTO 34
if (lastid .ne. luk)  GOTO 28
!  ******  if so, ju can be compressed  ********************************
irul = iru(lasti)
iju(k) = irul + 1
if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
GOTO 34
!  ******  if not, see if kth row can overlap the previous one  ********
28    if (jumin .gt. juptr)  GOTO 32
qm = q(qm)
do 29 j=jumin,juptr
if (ju(j) - qm)  29, 30, 32
29      continue
GOTO 32
30    iju(k) = j
do 31 i=j,juptr
if (ju(i) .ne. qm)  GOTO 32
qm = q(qm)
if (qm .gt. n)  GOTO 34
31      continue
juptr = j - 1
!  ******  move row indices from q to ju, update vectors  **************
32    jumin = juptr + 1
iju(k) = jumin
if (luk .eq. 0)  GOTO 34
juptr = juptr + luk
if (juptr .gt. jumax)  GOTO 106
qm = q(np1)
do 33 j=jumin,juptr
qm = q(qm)
33        ju(j) = qm
34    iru(k) = iju(k)
iu(k+1) = iu(k) + luk
!
!  ******  update iru, jru  ********************************************
i = k
35      i1 = jru(i)
if (r(i) .lt. 0)  GOTO 36
rend = iju(i) + iu(i+1) - iu(i)
if (iru(i) .ge. rend)  GOTO 37
j = ju(iru(i))
jru(i) = jru(j)
jru(j) = i
GOTO 37
36      r(i) = -r(i)
37      i = i1
if (i .eq. 0)  GOTO 38
iru(i) = iru(i) + 1
GOTO 35
!
!  ******  update ira, jra, irac  **************************************
38    i = irac(k)
if (i .eq. 0)  GOTO 41
39      i1 = jra(i)
ira(i) = ira(i) + 1
if (ira(i) .ge. ia(r(i)+1))  GOTO 40
irai = ira(i)
jairai = ic(ja(irai))
if (jairai .gt. i)  GOTO 40
jra(i) = irac(jairai)
irac(jairai) = i
40      i = i1
if (i .ne. 0)  GOTO 39
41    continue
!
ijl(n) = jlptr
iju(n) = juptr
flag = 0
return
!
! ** error.. null row in a
101  flag = n + rk
return
! ** error.. duplicate entry in a
102  flag = 2*n + rk
return
! ** error.. insufficient storage for jl
103  flag = 3*n + k
return
! ** error.. null pivot
105  flag = 5*n + k
return
! ** error.. insufficient storage for ju
106  flag = 6*n + k
return
end subroutine nsfc

subroutine nnfc(n, r,c,ic, ia,ja,a, z, b,lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,row, tmp, irl,jrl, flag)
!*** subroutine nnfc
!*** numerical ldu-factorization of sparse nonsymmetric matrix and
!      solution of system of linear equations (compressed pointer
!      storage)
!
!
!       input variables..  n, r, c, ic, ia, ja, a, b,
!                          il, jl, ijl, lmax, iu, ju, iju, umax
!       output variables.. z, l, d, u, flag
!
!       parameters used internally..
! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
! nia   - jrl       of the factorization,  jrl(k)  points to the head
!       -           of a linked list in  jrl  of column indices j
!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!       -           points to the smallest i such that i .ge. k and
!       -           l(i,j)  is nonzero.
!       -           size of each = n.
! fia   - row   - holds intermediate values in calculation of  u and l.
!       -           size = n.
! fia   - tmp   - holds new right-hand side  b*  for solution of the
!       -           equation ux = b*.
!       -           size = n.
!
!  internal variables..
!    jmin, jmax - indices of the first and last positions in a row to
!      be examined.
!    sum - used in calculating  tmp.
!
integer rk,umax
integer  r(*), c(*), ic(*), ia(*), ja(*), il(*), jl(*), ijl(*)
integer  iu(*), ju(*), iju(*), irl(*), jrl(*), flag
!     real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
!     real tmp(*), lki, sum, dk
double precision  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
double precision  tmp(*), lki, sum, dk
!
!  ******  initialize pointers and test storage  ***********************
if(il(n+1)-1 .gt. lmax) GOTO 104
if(iu(n+1)-1 .gt. umax) GOTO 107
do 1 k=1,n
irl(k) = il(k)
jrl(k) = 0
1    continue
!
!  ******  for each row  ***********************************************
do 19 k=1,n
!  ******  reverse jrl and zero row where kth row of l will fill in  ***
row(k) = 0
i1 = 0
if (jrl(k) .eq. 0) GOTO 3
i = jrl(k)
2    i2 = jrl(i)
jrl(i) = i1
i1 = i
row(i) = 0
i = i2
if (i .ne. 0) GOTO 2
!  ******  set row to zero where u will fill in  ***********************
3    jmin = iju(k)
jmax = jmin + iu(k+1) - iu(k) - 1
if (jmin .gt. jmax) GOTO 5
do 4 j=jmin,jmax
4      row(ju(j)) = 0
!  ******  place kth row of a in row  **********************************
5    rk = r(k)
jmin = ia(rk)
jmax = ia(rk+1) - 1
do 6 j=jmin,jmax
row(ic(ja(j))) = a(j)
6      continue
!  ******  initialize sum, and link through jrl  ***********************
sum = b(rk)
i = i1
if (i .eq. 0) GOTO 10
!  ******  assign the kth row of l and adjust row, sum  ****************
7      lki = -row(i)
!  ******  if l is not required, then comment out the following line  **
l(irl(i)) = -lki
sum = sum + lki * tmp(i)
jmin = iu(i)
jmax = iu(i+1) - 1
if (jmin .gt. jmax) GOTO 9
mu = iju(i) - jmin
do 8 j=jmin,jmax
8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
9      i = jrl(i)
if (i .ne. 0) GOTO 7
!
!  ******  assign kth row of u and diagonal d, set tmp(k)  *************
10    if (row(k) .eq. 0.0d0) GOTO 108
dk = 1.0d0 / row(k)
d(k) = dk
tmp(k) = sum * dk
if (k .eq. n) GOTO 19
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax)  GOTO 12
mu = iju(k) - jmin
do 11 j=jmin,jmax
11      u(j) = row(ju(mu+j)) * dk
12    continue
!
!  ******  update irl and jrl, keeping jrl in decreasing order  ********
i = i1
if (i .eq. 0) GOTO 18
14    irl(i) = irl(i) + 1
i1 = jrl(i)
if (irl(i) .ge. il(i+1)) GOTO 17
ijlb = irl(i) - il(i) + ijl(i)
j = jl(ijlb)
15    if (i .gt. jrl(j)) GOTO 16
j = jrl(j)
GOTO 15
16    jrl(i) = jrl(j)
jrl(j) = i
17    i = i1
if (i .ne. 0) GOTO 14
18    if (irl(k) .ge. il(k+1)) GOTO 19
j = jl(ijl(k))
jrl(k) = jrl(j)
jrl(j) = k
19    continue
!
!  ******  solve  ux = tmp  by back substitution  **********************
k = n
do 22 i=1,n
sum =  tmp(k)
jmin = iu(k)
jmax = iu(k+1) - 1
if (jmin .gt. jmax)  GOTO 21
mu = iju(k) - jmin
do 20 j=jmin,jmax
20      sum = sum - u(j) * tmp(ju(mu+j))
21    tmp(k) =  sum
z(c(k)) =  sum
22    k = k-1
flag = 0
return
!
! ** error.. insufficient storage for l
104  flag = 4*n + 1
return
! ** error.. insufficient storage for u
107  flag = 7*n + 1
return
! ** error.. zero pivot
108  flag = 8*n + k
return
end subroutine nnfc

SUBROUTINE DPRJS (NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WK,IWK,F,JAC)
EXTERNAL F,JAC
INTEGER NEQ, NYH, IWK
DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WK
DIMENSION Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, &
MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN,&
IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, &
NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMUL, J, JJ, JOK, JMAX, JMIN, K, KMAX, KMIN, NG
DOUBLE PRECISION CON, DI, FAC, HL0, PIJ, R, R0, RCON, RCONT, SRUR
!-----------------------------------------------------------------------
! DPRJS is called to compute and process the matrix
! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
! J is computed by columns, either by the user-supplied routine JAC
! if MITER = 1, or by finite differencing if MITER = 2.
! if MITER = 3, a diagonal approximation to J is used.
! if MITER = 1 or 2, and if the existing value of the Jacobian
! (as contained in P) is considered acceptable, then a new value of
! P is reconstructed from the old value.  In any case, when MITER
! is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
! P and its LU decomposition are stored (separately) in WK.
!
! In addition to variables described previously, communication
! with DPRJS uses the following:
! Y     = array containing predicted values on entry.
! FTEM  = work array of length N (ACOR in DSTODE).
! SAVF  = array containing f evaluated at predicted y.
! WK    = real work space for matrices.  On output it contains the
!         inverse diagonal matrix if MITER = 3, and P and its sparse
!         LU decomposition if MITER is 1 or 2.
!         Storage of matrix elements starts at WK(3).
!         WK also contains the following matrix-related data:
!         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
!         WK(2) = H*EL0, saved for later use if MITER = 3.
! IWK   = integer work space for matrix-related data, assumed to
!         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!         are assumed to have identical locations.
! EL0   = EL(1) (input).
! IERPJ = output error flag (in Common).
!       = 0 if no error.
!       = 1  if zero pivot found in CDRV.
!       = 2  if a singular matrix arose with MITER = 3.
!       = -1 if insufficient storage for CDRV (should not occur here).
!       = -2 if other error found in CDRV (should not occur here).
! JCUR  = output flag showing status of (approximate) Jacobian matrix:
!          = 1 to indicate that the Jacobian is now current, or
!          = 0 to indicate that a saved value was used.
! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
HL0 = H*EL0
CON = -HL0
IF (MITER .EQ. 3) GOTO 300
! See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
JOK = 1
IF (NST .EQ. 0 .OR. NST .GE. NSLJ+MSBJ) JOK = 0
IF (ICF .EQ. 1 .AND. ABS(RC - 1.0D0) .LT. CCMXJ) JOK = 0
IF (ICF .EQ. 2) JOK = 0
IF (JOK .EQ. 1) GOTO 250
!
! MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
20   JCUR = 1
NJE = NJE + 1
NSLJ = NST
IPLOST = 0
CONMIN = ABS(CON)
GOTO (100, 200), MITER
!
! If MITER = 1, call JAC, multiply by scalar, and add identity. --------
100  CONTINUE
KMIN = IWK(IPIAN)
DO 130 J = 1, N
KMAX = IWK(IPIAN+J) - 1
DO 110 I = 1,N
110      FTEM(I) = 0.0D0
CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), FTEM)
DO 120 K = KMIN, KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = FTEM(I)*CON
IF (I .EQ. J) WK(IBA+K) = WK(IBA+K) + 1.0D0
120      CONTINUE
KMIN = KMAX + 1
130    CONTINUE
GOTO 290
!
! If MITER = 2, make NGP calls to F to approximate J and P. ------------
200  CONTINUE
FAC = DVNORM(N, SAVF, EWT)
R0 = 1000.0D0 * ABS(H) * UROUND * N * FAC
IF (R0 .EQ. 0.0D0) R0 = 1.0D0
SRUR = WK(1)
JMIN = IWK(IPIGP)
DO 240 NG = 1,NGP
JMAX = IWK(IPIGP+NG) - 1
DO 210 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
210      Y(JJ) = Y(JJ) + R
CALL F (NEQ, TN, Y, FTEM)
DO 230 J = JMIN,JMAX
JJ = IWK(IBJGP+J)
Y(JJ) = YH(JJ,1)
R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
FAC = -HL0/R
KMIN =IWK(IBIAN+JJ)
KMAX =IWK(IBIAN+JJ+1) - 1
DO 220 K = KMIN,KMAX
I = IWK(IBJAN+K)
WK(IBA+K) = (FTEM(I) - SAVF(I))*FAC
IF (I .EQ. JJ) WK(IBA+K) = WK(IBA+K) + 1.0D0
220        CONTINUE
230      CONTINUE
JMIN = JMAX + 1
240    CONTINUE
NFE = NFE + NGP
GOTO 290
!
! If JOK = 1, reconstruct new P from old P. ----------------------------
250  JCUR = 0
RCON = CON/CON0
RCONT = ABS(CON)/CONMIN
IF (RCONT .GT. RBIG .AND. IPLOST .EQ. 1) GOTO 20
KMIN = IWK(IPIAN)
DO 275 J = 1,N
KMAX = IWK(IPIAN+J) - 1
DO 270 K = KMIN,KMAX
I = IWK(IBJAN+K)
PIJ = WK(IBA+K)
IF (I .NE. J) GOTO 260
PIJ = PIJ - 1.0D0
IF (ABS(PIJ) .GE. PSMALL) GOTO 260
IPLOST = 1
CONMIN = MIN(ABS(CON0),CONMIN)
260      PIJ = PIJ*RCON
IF (I .EQ. J) PIJ = PIJ + 1.0D0
WK(IBA+K) = PIJ
270      CONTINUE
KMIN = KMAX + 1
275    CONTINUE
!
! Do numerical factorization of P matrix. ------------------------------
290  NLU = NLU + 1
CON0 = CON
IERPJ = 0
DO 295 I = 1,N
295    FTEM(I) = 0.0D0
CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),FTEM,FTEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
IF (IYS .EQ. 0) RETURN
IMUL = (IYS - 1)/N
IERPJ = -2
IF (IMUL .EQ. 8) IERPJ = 1
IF (IMUL .EQ. 10) IERPJ = -1
RETURN
!
! If MITER = 3, construct a diagonal approximation to J and P. ---------
300  CONTINUE
JCUR = 1
NJE = NJE + 1
WK(2) = HL0
IERPJ = 0
R = EL0*0.1D0
DO 310 I = 1,N
310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
CALL F (NEQ, TN, Y, WK(3))
NFE = NFE + 1
DO 320 I = 1,N
R0 = H*SAVF(I) - YH(I,2)
DI = 0.1D0*R0 - H*(WK(I+2) - SAVF(I))
WK(I+2) = 1.0D0
IF (ABS(R0) .LT. UROUND/EWT(I)) GOTO 320
IF (ABS(DI) .EQ. 0.0D0) GOTO 330
WK(I+2) = 0.1D0*R0/DI
320    CONTINUE
RETURN
330  IERPJ = 2
RETURN
!----------------------- End of Subroutine DPRJS -----------------------
END subroutine dprjs

DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
!***BEGIN PROLOGUE  DVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the array V, with weights
!  contained in the array W of length N:
!    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DVNORM
!**End
INTEGER N,   I
DOUBLE PRECISION V, W,   SUM
DIMENSION V(N), W(N)
!
!***FIRST EXECUTABLE STATEMENT  DVNORM
SUM = 0.0D0
DO 10 I = 1,N
10     SUM = SUM + (V(I)*W(I))**2
DVNORM = SQRT(SUM/N)
RETURN
!----------------------- END OF FUNCTION DVNORM ------------------------
END function dvnorm

SUBROUTINE DSTODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM, F, JAC, PJAC, SLVS)
!***BEGIN PROLOGUE  DSTODE
!***SUBSIDIARY
!***PURPOSE  Performs one step of an ODEPACK integration.
!***TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DSTODE performs one step of the integration of an initial value
!  problem for a system of ordinary differential equations.
!  Note:  DSTODE is independent of the value of the iteration method
!  indicator MITER, when this is .ne. 0, and hence is independent
!  of the type of chord method used, or the Jacobian structure.
!  Communication with DSTODE is done with the following variables:
!
!  NEQ    = integer array containing problem size in NEQ(1), and
!           passed as the NEQ argument in all calls to F and JAC.
!  Y      = an array of length .ge. N used as the Y argument in
!           all calls to F and JAC.
!  YH     = an NYH by LMAX array containing the dependent variables
!           and their approximate scaled derivatives, where
!           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!           j-th derivative of y(i), scaled by h**j/factorial(j)
!           (j = 0,1,...,NQ).  on entry for the first step, the first
!           two columns of YH must be set from the initial values.
!  NYH    = a constant integer .ge. N, the first dimension of YH.
!  YH1    = a one-dimensional array occupying the same space as YH.
!  EWT    = an array of length N containing multiplicative weights
!           for local error measurements.  Local errors in Y(i) are
!           compared to 1.0/EWT(i) in various error tests.
!  SAVF   = an array of working storage, of length N.
!           Also used for input of YH(*,MAXORD+2) when JSTART = -1
!           and MAXORD .lt. the current order NQ.
!  ACOR   = a work array of length N, used for the accumulated
!           corrections.  On a successful return, ACOR(i) contains
!           the estimated one-step local error in Y(i).
!  WM,IWM = real and integer work arrays associated with matrix
!           operations in chord iteration (MITER .ne. 0).
!  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
!           and P = I - h*el0*JAC, if a chord method is being used.
!  SLVS   = name of routine to solve linear system in chord iteration.
!  CCMAX  = maximum relative change in h*el0 before PJAC is called.
!  H      = the step size to be attempted on the next step.
!           H is altered by the error control algorithm during the
!           problem.  H can be either positive or negative, but its
!           sign must remain constant throughout the problem.
!  HMIN   = the minimum absolute value of the step size h to be used.
!  HMXI   = inverse of the maximum absolute value of h to be used.
!           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
!           HMIN and HMXI may be changed at any time, but will not
!           take effect until the next change of h is considered.
!  TN     = the independent variable. TN is updated on each step taken.
!  JSTART = an integer used for input only, with the following
!           values and meanings:
!                0  perform the first step.
!            .gt.0  take a new step continuing from the last.
!               -1  take the next step with a new value of H, MAXORD,
!                     N, METH, MITER, and/or matrix parameters.
!               -2  take the next step with a new value of H,
!                     but with other inputs unchanged.
!           On return, JSTART is set to 1 to facilitate continuation.
!  KFLAG  = a completion code with the following meanings:
!                0  the step was succesful.
!               -1  the requested error could not be achieved.
!               -2  corrector convergence could not be achieved.
!               -3  fatal error in PJAC or SLVS.
!           A return with KFLAG = -1 or -2 means either
!           abs(H) = HMIN or 10 consecutive failures occurred.
!           On a return with KFLAG negative, the values of TN and
!           the YH array are as of the beginning of the last
!           step, and H is the last step size attempted.
!  MAXORD = the maximum order of integration method to be allowed.
!  MAXCOR = the maximum number of corrector iterations allowed.
!  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
!  MXNCF  = maximum number of convergence failures allowed.
!  METH/MITER = the method flags.  See description in driver.
!  N      = the number of first-order differential equations.
!  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
!  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DCFODE, DVNORM
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSTODE
!**End
EXTERNAL F, JAC, PJAC, SLVS
INTEGER NEQ, NYH, IWM
DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
DIMENSION Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),ACOR(*), WM(*), IWM(*)
INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L
INTEGER LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD

COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, &
LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU

!
!***FIRST EXECUTABLE STATEMENT  DSTODE
KFLAG = 0
TOLD = TN
NCF = 0
IERPJ = 0
IERSL = 0
JCUR = 0
ICF = 0
DELP = 0.0D0
IF (JSTART .GT. 0) GOTO 200
IF (JSTART .EQ. -1) GOTO 100
IF (JSTART .EQ. -2) GOTO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set to 2
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
IALTH = 2
RMAX = 10000.0D0
RC = 0.0D0
EL0 = 1.0D0
CRATE = 0.7D0
HOLD = H
MEO = METH
NSLP = 0
IPUP = MITER
IRET = 3
GOTO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100  IPUP = MITER
LMAX = MAXORD + 1
IF (IALTH .EQ. 1) IALTH = 2
IF (METH .EQ. MEO) GOTO 110
CALL DCFODE (METH, ELCO, TESCO)
MEO = METH
IF (NQ .GT. MAXORD) GOTO 120
IALTH = L
IRET = 1
GOTO 150
110  IF (NQ .LE. MAXORD) GOTO 160
120  NQ = MAXORD
L = LMAX
DO 125 I = 1,L
125    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
EXDN = 1.0D0/L
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
RH = MIN(RHDN,1.0D0)
IREDO = 3
IF (H .EQ. HOLD) GOTO 170
RH = MIN(RH,ABS(H/HOLD))
H = HOLD
GOTO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140  CALL DCFODE (METH, ELCO, TESCO)
150  DO 155 I = 1,L
155    EL(I) = ELCO(I,NQ)
NQNYH = NQ*NYH
RC = RC*EL(1)/EL0
EL0 = EL(1)
CONIT = 0.5D0/(NQ+2)
GOTO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160  IF (H .EQ. HOLD) GOTO 200
RH = H/HOLD
H = HOLD
IREDO = 3
GOTO 175
170  RH = MAX(RH,HMIN/ABS(H))
175  RH = MIN(RH,RMAX)
RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
R = 1.0D0
DO 180 J = 2,L
R = R*RH
DO 180 I = 1,N
180      YH(I,J) = YH(I,J)*R
H = H*RH
RC = RC*RH
IALTH = L
IF (IREDO .EQ. 0) GOTO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal Triangle matrix.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called, if a Jacobian is involved.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
IF (NST .GE. NSLP+MSBP) IPUP = MITER
TN = TN + H
I1 = NQNYH + 1
DO 215 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 210 I = I1,NQNYH
210      YH1(I) = YH1(I) + YH1(I+NYH)
215    CONTINUE
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the R.M.S. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220  M = 0
DO 230 I = 1,N
230    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
IF (IPUP .LE. 0) GOTO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*el(1)*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
IPUP = 0
RC = 1.0D0
NSLP = NST
CRATE = 0.7D0
IF (IERPJ .NE. 0) GOTO 430
250  DO 260 I = 1,N
260    ACOR(I) = 0.0D0
270  IF (MITER .NE. 0) GOTO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
DO 290 I = 1,N
SAVF(I) = H*SAVF(I) - YH(I,2)
290    Y(I) = SAVF(I) - ACOR(I)
DEL = DVNORM (N, Y, EWT)
DO 300 I = 1,N
Y(I) = YH(I,1) + EL(1)*SAVF(I)
300    ACOR(I) = SAVF(I)
GOTO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
350  DO 360 I = 1,N
360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
CALL SLVS (WM, IWM, Y, SAVF)
IF (IERSL .LT. 0) GOTO 430
IF (IERSL .GT. 0) GOTO 410
DEL = DVNORM (N, Y, EWT)
DO 380 I = 1,N
ACOR(I) = ACOR(I) + Y(I)
380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M.gt.0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
IF (DCON .LE. 1.0D0) GOTO 450
M = M + 1
IF (M .EQ. MAXCOR) GOTO 410
IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GOTO 410
DELP = DEL
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
GOTO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GOTO 430
ICF = 1
IPUP = MITER
GOTO 220
430  ICF = 2
NCF = NCF + 1
RMAX = 2.0D0
TN = TOLD
I1 = NQNYH + 1
DO 445 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 440 I = I1,NQNYH
440      YH1(I) = YH1(I) - YH1(I+NYH)
445    CONTINUE
IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GOTO 680
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 670
IF (NCF .EQ. MXNCF) GOTO 670
RH = 0.25D0
IPUP = MITER
IREDO = 1
GOTO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450  JCUR = 0
IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
IF (DSM .GT. 1.0D0) GOTO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
KFLAG = 0
IREDO = 0
NST = NST + 1
HU = H
NQU = NQ
DO 470 J = 1,L
DO 470 I = 1,N
470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
IALTH = IALTH - 1
IF (IALTH .EQ. 0) GOTO 520
IF (IALTH .GT. 1) GOTO 700
IF (L .EQ. LMAX) GOTO 700
DO 490 I = 1,N
490    YH(I,LMAX) = ACOR(I)
GOTO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500  KFLAG = KFLAG - 1
TN = TOLD
I1 = NQNYH + 1
DO 515 JB = 1,NQ
I1 = I1 - NYH
!dir$ ivdep
DO 510 I = I1,NQNYH
510      YH1(I) = YH1(I) - YH1(I+NYH)
515    CONTINUE
RMAX = 2.0D0
IF (ABS(H) .LE. HMIN*1.00001D0) GOTO 660
IF (KFLAG .LE. -3) GOTO 640
IREDO = 2
RHUP = 0.0D0
GOTO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520  RHUP = 0.0D0
IF (L .EQ. LMAX) GOTO 540
DO 530 I = 1,N
530    SAVF(I) = ACOR(I) - YH(I,LMAX)
DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
EXUP = 1.0D0/(L+1)
RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
540  EXSM = 1.0D0/L
RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
RHDN = 0.0D0
IF (NQ .EQ. 1) GOTO 560
DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
EXDN = 1.0D0/NQ
RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
560  IF (RHSM .GE. RHUP) GOTO 570
IF (RHUP .GT. RHDN) GOTO 590
GOTO 580
570  IF (RHSM .LT. RHDN) GOTO 580
NEWQ = NQ
RH = RHSM
GOTO 620
580  NEWQ = NQ - 1
RH = RHDN
IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
GOTO 620
590  NEWQ = L
RH = RHUP
IF (RH .LT. 1.1D0) GOTO 610
R = EL(L)/L
DO 600 I = 1,N
600    YH(I,NEWQ+1) = ACOR(I)*R
GOTO 630
610  IALTH = 3
GOTO 700
620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GOTO 610
IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, l, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
IF (NEWQ .EQ. NQ) GOTO 170
630  NQ = NEWQ
L = NQ + 1
IRET = 2
GOTO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640  IF (KFLAG .EQ. -10) GOTO 660
RH = 0.1D0
RH = MAX(HMIN/ABS(H),RH)
H = H*RH
DO 645 I = 1,N
645    Y(I) = YH(I,1)
CALL F (NEQ, TN, Y, SAVF)
NFE = NFE + 1
DO 650 I = 1,N
650    YH(I,2) = H*SAVF(I)
IPUP = MITER
IALTH = 5
IF (NQ .EQ. 1) GOTO 200
NQ = 1
L = 2
IRET = 3
GOTO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660  KFLAG = -1
GOTO 720
670  KFLAG = -2
GOTO 720
680  KFLAG = -3
GOTO 720
690  RMAX = 10.0D0
700  R = 1.0D0/TESCO(2,NQU)
DO 710 I = 1,N
710    ACOR(I) = ACOR(I)*R
720  HOLD = H
JSTART = 1
RETURN
!----------------------- END OF SUBROUTINE DSTODE ----------------------
END subroutine dstode

SUBROUTINE DCFODE (METH, ELCO, TESCO)
!***BEGIN PROLOGUE  DCFODE
!***SUBSIDIARY
!***PURPOSE  Set ODE integrator coefficients.
!***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DCFODE is called by the integrator routine to set coefficients
!  needed there.  The coefficients for the current method, as
!  given by the value of METH, are set for all orders and saved.
!  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
!  (A smaller value of the maximum order is also allowed.)
!  DCFODE is called once at the beginning of the problem,
!  and is not called again unless and until METH is changed.
!
!  The ELCO array contains the basic method coefficients.
!  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
!  order nq are stored in ELCO(i,nq).  They are given by a genetrating
!  polynomial, i.e.,
!      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
!  For the implicit Adams methods, l(x) is given by
!      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
!  For the BDF methods, l(x) is given by
!      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
!  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
!  The TESCO array contains test constants used for the
!  local error test and the selection of step size and/or order.
!  At order nq, TESCO(k,nq) is used for the selection of step
!  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
!  nq + 1 if k = 3.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DCFODE
!**End
INTEGER METH
INTEGER I, IB, NQ, NQM1, NQP1
DOUBLE PRECISION ELCO, TESCO
DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,RQFAC, RQ1FAC, TSIGN, XPIN
DIMENSION ELCO(13,12), TESCO(3,12)
DIMENSION PC(12)
!
!***FIRST EXECUTABLE STATEMENT  DCFODE
GOTO (100, 200), METH
!
100  ELCO(1,1) = 1.0D0
ELCO(2,1) = 1.0D0
TESCO(1,1) = 0.0D0
TESCO(2,1) = 2.0D0
TESCO(1,2) = 1.0D0
TESCO(3,12) = 0.0D0
PC(1) = 1.0D0
RQFAC = 1.0D0
DO 140 NQ = 2,12
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
RQ1FAC = RQFAC
RQFAC = RQFAC/NQ
NQM1 = NQ - 1
FNQM1 = NQM1
NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq-1). ----------------------------------
PC(NQ) = 0.0D0
DO 110 IB = 1,NQM1
I = NQP1 - IB
110      PC(I) = PC(I-1) + FNQM1*PC(I)
PC(1) = FNQM1*PC(1)
! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
PINT = PC(1)
XPIN = PC(1)/2.0D0
TSIGN = 1.0D0
DO 120 I = 2,NQ
TSIGN = -TSIGN
PINT = PINT + TSIGN*PC(I)/I
120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
! Store coefficients in ELCO and TESCO. --------------------------------
ELCO(1,NQ) = PINT*RQ1FAC
ELCO(2,NQ) = 1.0D0
DO 130 I = 2,NQ
130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
AGAMQ = RQFAC*XPIN
RAGQ = 1.0D0/AGAMQ
TESCO(2,NQ) = RAGQ
IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
TESCO(3,NQM1) = RAGQ
140    CONTINUE
RETURN
!
200  PC(1) = 1.0D0
RQ1FAC = 1.0D0
DO 230 NQ = 1,5
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
FNQ = NQ
NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq). ------------------------------------
PC(NQP1) = 0.0D0
DO 210 IB = 1,NQ
I = NQ + 2 - IB
210      PC(I) = PC(I-1) + FNQ*PC(I)
PC(1) = FNQ*PC(1)
! Store coefficients in ELCO and TESCO. --------------------------------
DO 220 I = 1,NQP1
220      ELCO(I,NQ) = PC(I)/PC(2)
ELCO(2,NQ) = 1.0D0
TESCO(1,NQ) = RQ1FAC
TESCO(2,NQ) = NQP1/ELCO(1,NQ)
TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
RQ1FAC = RQ1FAC/FNQ
230    CONTINUE
RETURN
!----------------------- END OF SUBROUTINE DCFODE ----------------------
END subroutine dcfode

SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
!***BEGIN PROLOGUE  DEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DEWSET
!**End
INTEGER N, ITOL
INTEGER I
DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
DIMENSION ATOL(*), YCUR(N), EWT(N)
!
!***FIRST EXECUTABLE STATEMENT  DEWSET
GOTO (10, 20, 30, 40), ITOL
10   CONTINUE
DO 15 I = 1,N
15     EWT(I) = RTOL*ABS(YCUR(I)) + ATOL(1)
RETURN
20   CONTINUE
DO 25 I = 1,N
25     EWT(I) = RTOL*ABS(YCUR(I)) + ATOL(I)
RETURN
30   CONTINUE
DO 35 I = 1,N
35     EWT(I) = RTOL*ABS(YCUR(I)) + ATOL(1)
RETURN
40   CONTINUE
DO 45 I = 1,N
45     EWT(I) = RTOL*ABS(YCUR(I)) + ATOL(I)
RETURN
!----------------------- END OF SUBROUTINE DEWSET ----------------------
END subroutine dewset

SUBROUTINE DIPREP (NEQ, Y, RWORK, IA, JA, IPFLAG, F, JAC)
EXTERNAL F, JAC
INTEGER NEQ, IA, JA, IPFLAG, IWORK
DOUBLE PRECISION Y, RWORK
DIMENSION Y(*), RWORK(*), IA(*), JA(*), IWORK(LWM)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION RLSS
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ RLSS(6), IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, &
IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
INTEGER I, IMAX, LEWTN, LYHD, LYHN
!-----------------------------------------------------------------------
! This routine serves as an interface between the driver and
! Subroutine DPREP.  It is called only if MITER is 1 or 2.
! Tasks performed here are:
!  * call DPREP,
!  * reset the required WM segment length LENWK,
!  * move YH back to its final location (following WM in RWORK),
!  * reset pointers for YH, SAVF, EWT, and ACOR, and
!  * move EWT to its new position if ISTATE = 1.
! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
! no trouble, and IPFLAG is the value of the DPREP error flag IPPER
! if there was trouble in Subroutine DPREP.
!-----------------------------------------------------------------------
IPFLAG = 0
IWORK(1:LWM) = 0
! Call DPREP to do matrix preprocessing operations. --------------------
CALL DPREP (NEQ, Y, RWORK(LYH), RWORK(LSAVF), RWORK(LEWT), RWORK(LACOR), IA, JA, RWORK(LWM), IWORK(LWM), IPFLAG, F, JAC)
LENWK = MAX(LREQ,LWMIN)
IF (IPFLAG .LT. 0) RETURN
! If DPREP was successful, move YH to end of required space for WM. ----
LYHN = LWM + LENWK
IF (LYHN .GT. LYH) RETURN
LYHD = LYH - LYHN
IF (LYHD .EQ. 0) GOTO 20
IMAX = LYHN - 1 + LENYHM
DO 10 I = LYHN,IMAX
10     RWORK(I) = RWORK(I+LYHD)
LYH = LYHN
! Reset pointers for SAVF, EWT, and ACOR. ------------------------------
20   LSAVF = LYH + LENYH
LEWTN = LSAVF + N
LACOR = LEWTN + N
IF (ISTATC .EQ. 3) GOTO 40
! If ISTATE = 1, move EWT (left) to its new position. ------------------
IF (LEWTN .GT. LEWT) RETURN
DO 30 I = 1,N
30     RWORK(I+LEWTN-1) = RWORK(I+LEWT-1)
40   LEWT = LEWTN
RETURN
!----------------------- End of Subroutine DIPREP ----------------------
END subroutine diprep

SUBROUTINE DPREP (NEQ, Y, YH, SAVF, EWT, FTEM, IA, JA, WK, IWK, IPPER, F, JAC)
EXTERNAL F,JAC
INTEGER NEQ, IA, JA, IWK, IPPER
DOUBLE PRECISION Y, YH, SAVF, EWT, FTEM, WK
DIMENSION Y(*), YH(*), SAVF(*), EWT(*), FTEM(*), IA(*), JA(*), WK(*), IWK(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD
INTEGER MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA
INTEGER LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
DOUBLE PRECISION ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, &
JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, IPIAN, IPJAN,&
IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA, LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ, NSLJ, NGP, &
NLU, NNZ, NSP, NZL, NZU
INTEGER I, IBR, IER, IPIL, IPIU, IPTT1, IPTT2, J, JFOUND, K, KNEW, KMAX, KMIN, LDIF, LENIGP, LIWK, MAXG, NP1, NZSUT
DOUBLE PRECISION DQ, DYJ, ERWT, FAC, YJ
!-----------------------------------------------------------------------
! This routine performs preprocessing related to the sparse linear
! systems that must be solved if MITER = 1 or 2.
! The operations that are performed here are:
!  * compute sparseness structure of Jacobian according to MOSS,
!  * compute grouping of column indices (MITER = 2),
!  * compute a new ordering of rows and columns of the matrix,
!  * reorder JA corresponding to the new ordering,
!  * perform a symbolic LU factorization of the matrix, and
!  * set pointers for segments of the IWK/WK array.
! In addition to variables described previously, DPREP uses the
! following for communication:
! YH     = the history array.  Only the first column, containing the
!          current Y vector, is used.  Used only if MOSS .ne. 0.
! SAVF   = a work array of length NEQ, used only if MOSS .ne. 0.
! EWT    = array of length NEQ containing (inverted) error weights.
!          Used only if MOSS = 2 or if ISTATE = MOSS = 1.
! FTEM   = a work array of length NEQ, identical to ACOR in the driver,
!          used only if MOSS = 2.
! WK     = a real work array of length LENWK, identical to WM in
!          the driver.
! IWK    = integer work array, assumed to occupy the same space as WK.
! LENWK  = the length of the work arrays WK and IWK.
! ISTATC = a copy of the driver input argument ISTATE (= 1 on the
!          first call, = 3 on a continuation call).
! IYS    = flag value from ODRV or CDRV.
! IPPER  = output error flag with the following values and meanings:
!          0  no error.
!         -1  insufficient storage for internal structure pointers.
!         -2  insufficient storage for JGROUP.
!         -3  insufficient storage for ODRV.
!         -4  other error flag from ODRV (should never occur).
!         -5  insufficient storage for CDRV.
!         -6  other error flag from CDRV.
!-----------------------------------------------------------------------
IBIAN = LRAT*2
IPIAN = IBIAN + 1
NP1 = N + 1
IPJAN = IPIAN + NP1
IBJAN = IPJAN - 1
LIWK = LENWK*LRAT
IF (IPJAN+N-1 .GT. LIWK) GOTO 210
IF (MOSS .EQ. 0) GOTO 30
!
IF (ISTATC .EQ. 3) GOTO 20
! ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
DO 10 I = 1,N
ERWT = 1.0D0/EWT(I)
FAC = 1.0D0 + 1.0D0/(I + 1.0D0)
Y(I) = Y(I) + FAC*SIGN(ERWT,Y(I))
10     CONTINUE
GOTO (70, 100), MOSS
!
20   CONTINUE
! ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
DO 25 I = 1,N
25     Y(I) = YH(I)
GOTO (70, 100), MOSS
!
! MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
30   KNEW = IPJAN
KMIN = IA(1)
IWK(IPIAN) = 1
DO 60 J = 1,N
JFOUND = 0
KMAX = IA(J+1) - 1
IF (KMIN .GT. KMAX) GOTO 45
DO 40 K = KMIN,KMAX
I = JA(K)
IF (I .EQ. J) JFOUND = 1
IF (KNEW .GT. LIWK) GOTO 210
IWK(KNEW) = I
KNEW = KNEW + 1
40       CONTINUE
IF (JFOUND .EQ. 1) GOTO 50
45     IF (KNEW .GT. LIWK) GOTO 210
IWK(KNEW) = J
KNEW = KNEW + 1
50     IWK(IPIAN+J) = KNEW + 1 - IPJAN
KMIN = KMAX + 1
60     CONTINUE
GOTO 140
!
! MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
70   CONTINUE
! A dummy call to F allows user to create temporaries for use in JAC. --
CALL F (NEQ, TN, Y, SAVF)
K = IPJAN
IWK(IPIAN) = 1
DO 90 J = 1,N
IF (K .GT. LIWK) GOTO 210
IWK(K) = J
K = K + 1
DO 75 I = 1,N
75       SAVF(I) = 0.0D0
CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), SAVF)
DO 80 I = 1,N
IF (ABS(SAVF(I)) .LE. SETH) GOTO 80
IF (I .EQ. J) GOTO 80
IF (K .GT. LIWK) GOTO 210
IWK(K) = I
K = K + 1
80       CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
90     CONTINUE
GOTO 140
!
! MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
100  K = IPJAN
IWK(IPIAN) = 1
CALL F (NEQ, TN, Y, SAVF)
DO 120 J = 1,N
IF (K .GT. LIWK) GOTO 210
IWK(K) = J
K = K + 1
YJ = Y(J)
ERWT = 1.0D0/EWT(J)
DYJ = SIGN(ERWT,YJ)
Y(J) = YJ + DYJ
CALL F (NEQ, TN, Y, FTEM)
Y(J) = YJ
DO 110 I = 1,N
DQ = (FTEM(I) - SAVF(I))/DYJ
IF (ABS(DQ) .LE. SETH) GOTO 110
IF (I .EQ. J) GOTO 110
IF (K .GT. LIWK) GOTO 210
IWK(K) = I
K = K + 1
110      CONTINUE
IWK(IPIAN+J) = K + 1 - IPJAN
120    CONTINUE
!
140  CONTINUE
IF (MOSS .EQ. 0 .OR. ISTATC .NE. 1) GOTO 150
! If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
DO 145 I = 1,N
145    Y(I) = YH(I)
150  NNZ = IWK(IPIAN+N) - 1
LENIGP = 0
IPIGP = IPJAN + NNZ
IF (MITER .NE. 2) GOTO 160
!
! Compute grouping of column indices (MITER = 2). ----------------------
MAXG = NP1
IPJGP = IPJAN + NNZ
IBJGP = IPJGP - 1
IPIGP = IPJGP + N
IPTT1 = IPIGP + NP1
IPTT2 = IPTT1 + N
LREQ = IPTT2 + N - 1
IF (LREQ .GT. LIWK) GOTO 220
CALL JGROUP (N, IWK(IPIAN), IWK(IPJAN), MAXG, NGP, IWK(IPIGP), IWK(IPJGP), IWK(IPTT1), IWK(IPTT2), IER)
IF (IER .NE. 0) GOTO 220
LENIGP = NGP + 1
!
! Compute new ordering of rows/columns of Jacobian. --------------------
160  IPR = IPIGP + LENIGP
IPC = IPR
IPIC = IPC + N
IPISP = IPIC + N
IPRSP = (IPISP - 2)/LRAT + 2
IESP = LENWK + 1 - IPRSP
IF (IESP .LT. 0) GOTO 230
IBR = IPR - 1
DO 170 I = 1,N
170    IWK(IBR+I) = I
NSP = LIWK + 1 - IPISP
CALL ODRV (N, IWK(IPIAN), IWK(IPJAN), WK, IWK(IPR), IWK(IPIC), NSP, IWK(IPISP), 1, IYS)
IF (IYS .EQ. 11*N+1) GOTO 240
IF (IYS .NE. 0) GOTO 230
!
! Reorder JAN and do symbolic LU factorization of matrix. --------------
IPA = LENWK + 1 - NNZ
NSP = IPA - IPRSP
LREQ = MAX(12*N/LRAT, 6*N/LRAT+2*N+NNZ) + 3
LREQ = LREQ + IPRSP - 1 + NNZ
IF (LREQ .GT. LENWK) GOTO 250
IBA = IPA - 1
DO 180 I = 1,NNZ
180    WK(IBA+I) = 0.0D0
IPISP = LRAT*(IPRSP - 1) + 1
CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN), WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)
LREQ = LENWK - IESP
IF (IYS .EQ. 10*N+1) GOTO 250
IF (IYS .NE. 0) GOTO 260
IPIL = IPISP
IPIU = IPIL + 2*N + 1
NZU = IWK(IPIL+N) - IWK(IPIL)
NZL = IWK(IPIU+N) - IWK(IPIU)
IF (LRAT .GT. 1) GOTO 190
CALL ADJLR (N, IWK(IPISP), LDIF)
LREQ = LREQ + LDIF
190  CONTINUE
IF (LRAT .EQ. 2 .AND. NNZ .EQ. N) LREQ = LREQ + 1
NSP = NSP + LREQ - LENWK
IPA = LREQ + 1 - NNZ
IBA = IPA - 1
IPPER = 0
RETURN
!
210  IPPER = -1
LREQ = 2 + (2*N + 1)/LRAT
LREQ = MAX(LENWK+1,LREQ)
RETURN
!
220  IPPER = -2
LREQ = (LREQ - 1)/LRAT + 1
RETURN
!
230  IPPER = -3
CALL CNTNZU (N, IWK(IPIAN), IWK(IPJAN), NZSUT)
LREQ = LENWK - IESP + (3*N + 4*NZSUT - 1)/LRAT + 1
RETURN
!
240  IPPER = -4
RETURN
!
250  IPPER = -5
RETURN
!
260  IPPER = -6
LREQ = LENWK
RETURN
!----------------------- End of Subroutine DPREP -----------------------
END subroutine dprep

SUBROUTINE JGROUP (N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
INTEGER N, IA, JA, MAXG, NGRP, IGP, JGP, INCL, JDONE, IER
DIMENSION IA(*), JA(*), IGP(*), JGP(*), INCL(*), JDONE(*)
!-----------------------------------------------------------------------
! This subroutine constructs groupings of the column indices of
! the Jacobian matrix, used in the numerical evaluation of the
! Jacobian by finite differences.
!
! Input:
! N      = the order of the matrix.
! IA,JA  = sparse structure descriptors of the matrix by rows.
! MAXG   = length of available storage in the IGP array.
!
! Output:
! NGRP   = number of groups.
! JGP    = array of length N containing the column indices by groups.
! IGP    = pointer array of length NGRP + 1 to the locations in JGP
!          of the beginning of each group.
! IER    = error indicator.  IER = 0 if no error occurred, or 1 if
!          MAXG was insufficient.
!
! INCL and JDONE are working arrays of length N.
!-----------------------------------------------------------------------
INTEGER I, J, K, KMIN, KMAX, NCOL, NG
!
IER = 0
DO 10 J = 1,N
10     JDONE(J) = 0
NCOL = 1
DO 60 NG = 1,MAXG
IGP(NG) = NCOL
DO 20 I = 1,N
20       INCL(I) = 0
DO 50 J = 1,N
! Reject column J if it is already in a group.--------------------------
IF (JDONE(J) .EQ. 1) GOTO 50
KMIN = IA(J)
KMAX = IA(J+1) - 1
DO 30 K = KMIN,KMAX
! Reject column J if it overlaps any column already in this group.------
I = JA(K)
IF (INCL(I) .EQ. 1) GOTO 50
30         CONTINUE
! Accept column J into group NG.----------------------------------------
JGP(NCOL) = J
NCOL = NCOL + 1
JDONE(J) = 1
DO 40 K = KMIN,KMAX
I = JA(K)
40         INCL(I) = 1
50       CONTINUE
! Stop if this group is empty (grouping is complete).-------------------
IF (NCOL .EQ. IGP(NG)) GOTO 70
60     CONTINUE
! Error return if not all columns were chosen (MAXG too small).---------
IF (NCOL .LE. N) GOTO 80
NG = MAXG
70   NGRP = NG - 1
RETURN
80   IER = 1
RETURN
!----------------------- End of Subroutine JGROUP ----------------------
END subroutine jgroup

subroutine odrv(n, ia,ja,a, p,ip, nsp,isp, path, flag)
!                                                                 5/2/83
!***********************************************************************
!  odrv -- driver for sparse matrix reordering routines
!***********************************************************************
!
!  description
!
!    odrv finds a minimum degree ordering of the rows and columns
!    of a matrix m stored in (ia,ja,a) format (see below).  for the
!    reordered matrix, the work and storage required to perform
!    gaussian elimination is (usually) significantly less.
!
!    note.. odrv and its subordinate routines have been modified to
!    compute orderings for general matrices, not necessarily having any
!    symmetry.  the miminum degree ordering is computed for the
!    structure of the symmetric matrix  m + m-transpose.
!    modifications to the original odrv module have been made in
!    the coding in subroutine mdi, and in the initial comments in
!    subroutines odrv and md.
!
!    if only the nonzero entries in the upper triangle of m are being
!    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
!    with the diagonal entries placed first in each row.  this is to
!    ensure that if m(i,j) will be in the upper triangle of m with
!    respect to the new ordering, then m(i,j) is stored in row i (and
!    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
!    strict lower triangle of m, then m(j,i) is stored in row j (and
!    thus m(i,j) is not stored).
!
!
!  storage of sparse matrices
!
!    the nonzero entries of the matrix m are stored row-by-row in the
!    array a.  to identify the individual nonzero entries in each row,
!    we need to know in which column each entry lies.  these column
!    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  to identify the individual rows, we need to know where
!    each row starts.  these row pointers are stored in the array ia.
!    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
!    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
!    the first location following the last element in the last row.
!    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
!    the nonzero entries in the i-th row are stored consecutively in
!
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!
!    and the corresponding column indices are stored consecutively in
!
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!
!    when the coefficient matrix is symmetric, only the nonzero entries
!    in the upper triangle need be stored.  for example, the matrix
!
!             ( 1  0  2  3  0 )
!             ( 0  4  0  0  0 )
!         m = ( 2  0  5  6  0 )
!             ( 3  0  6  7  8 )
!             ( 0  0  0  8  9 )
!
!    could be stored as
!
!            - 1  2  3  4  5  6  7  8  9 10 11 12 13
!         ---+--------------------------------------
!         ia - 1  4  5  8 12 14
!         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
!          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
!
!    or (symmetrically) as
!
!            - 1  2  3  4  5  6  7  8  9
!         ---+--------------------------
!         ia - 1  4  5  7  9 10
!         ja - 1  3  4  2  3  4  4  5  5
!          a - 1  2  3  4  5  6  7  8  9          .
!
!
!  parameters
!
!    n    - order of the matrix
!
!    ia   - integer one-dimensional array containing pointers to delimit
!           rows in ja and a.  dimension = n+1
!
!    ja   - integer one-dimensional array containing the column indices
!           corresponding to the elements of a.  dimension = number of
!           nonzero entries in (the upper triangle of) m
!
!    a    - real one-dimensional array containing the nonzero entries in
!           (the upper triangle of) m, stored by rows.  dimension =
!           number of nonzero entries in (the upper triangle of) m
!
!    p    - integer one-dimensional array used to return the permutation
!           of the rows and columns of m corresponding to the minimum
!           degree ordering.  dimension = n
!
!    ip   - integer one-dimensional array used to return the inverse of
!           the permutation returned in p.  dimension = n
!
!    nsp  - declared dimension of the one-dimensional array isp.  nsp
!           must be at least  3n+4k,  where k is the number of nonzeroes
!           in the strict upper triangle of m
!
!    isp  - integer one-dimensional array used for working storage.
!           dimension = nsp
!
!    path - integer path specification.  values and their meanings are -
!             1  find minimum degree ordering only
!             2  find minimum degree ordering and reorder symmetrically
!                  stored matrix (used when only the nonzero entries in
!                  the upper triangle of m are being stored)
!             3  reorder symmetrically stored matrix as specified by
!                  input permutation (used when an ordering has already
!                  been determined and only the nonzero entries in the
!                  upper triangle of m are being stored)
!             4  same as 2 but put diagonal entries at start of each row
!             5  same as 3 but put diagonal entries at start of each row
!
!    flag - integer error flag.  values and their meanings are -
!               0    no errors detected
!              9n+k  insufficient storage in md
!             10n+1  insufficient storage in odrv
!             11n+1  illegal path specification
!
!
!  conversion from real to double precision
!
!    change the real declarations in odrv and sro to double precision
!    declarations.
!
!-----------------------------------------------------------------------
!
integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,v, l, head,  tmp, q
integer n
!...  real  a(*)
double precision  a(*)
logical  dflag
!
!----initialize error flag and validate path specification
flag = 0
if (path.lt.1 .or. 5.lt.path)  GOTO 111
!
!----allocate storage and find minimum degree ordering
if ((path-1) * (path-2) * (path-4) .ne. 0)  GOTO 1
max = (nsp-n)/2
v    = 1
l    = v     +  max
head = l     +  max
next = head  +  n
if (max.lt.n)  GOTO 110
!
call  md(n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
if (flag.ne.0)  GOTO 100
!
!----allocate storage and symmetrically reorder matrix
1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  GOTO 2
tmp = (nsp+1) -      n
q   = tmp     - (ia(n+1)-1)
if (q.lt.1)  GOTO 110
!
dflag = path.eq.4 .or. path.eq.5
call sro(n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
!
2  return
!
! ** error -- error detected in md
100  return
! ** error -- insufficient storage
110  flag = 10*n + 1
return
! ** error -- illegal path specified
111  flag = 11*n + 1
return
end subroutine odrv

subroutine sro(n, ip, ia,ja,a, q, r, dflag)
!***********************************************************************
!  sro -- symmetric reordering of sparse symmetric matrix
!***********************************************************************
!
!  description
!
!    the nonzero entries of the matrix m are assumed to be stored
!    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
!    are stored if i ne j).
!
!    sro does not rearrange the order of the rows, but does move
!    nonzeroes from one row to another to ensure that if m(i,j) will be
!    in the upper triangle of m with respect to the new ordering, then
!    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
!    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
!    stored in row j (and thus m(i,j) is not stored).
!
!
!  additional parameters
!
!    q     - integer one-dimensional work array.  dimension = n
!
!    r     - integer one-dimensional work array.  dimension = number of
!            nonzero entries in the upper triangle of m
!
!    dflag - logical variable.  if dflag = .true., then store nonzero
!            diagonal elements at the beginning of the row
!
!-----------------------------------------------------------------------
!
integer  ip(*),  ia(*), ja(*),  q(*), r(*)
!...  real  a(*),  ak
double precision  a(*),  ak
logical  dflag
!
!
!--phase 1 -- find row in which to store each nonzero
!----initialize count of nonzeroes to be stored in each row
do 1 i=1,n
1     q(i) = 0
!
!----for each nonzero element a(j)
do 3 i=1,n
jmin = ia(i)
jmax = ia(i+1) - 1
if (jmin.gt.jmax)  GOTO 3
do 2 j=jmin,jmax
!
!--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
k = ja(j)
if (ip(k).lt.ip(i))  ja(j) = i
if (ip(k).ge.ip(i))  k = i
r(j) = k
!
!--------... and increment count of nonzeroes (=q(r(j)) in that row
2       q(k) = q(k) + 1
3     continue
!
!
!--phase 2 -- find new ia and permutation to apply to (ja,a)
!----determine pointers to delimit rows in permuted (ja,a)
do 4 i=1,n
ia(i+1) = ia(i) + q(i)
4     q(i) = ia(i+1)
!
!----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
!----for each nonzero element (in reverse order)
ilast = 0
jmin = ia(1)
jmax = ia(n+1) - 1
j = jmax
do 6 jdummy=jmin,jmax
i = r(j)
if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  GOTO 5
!
!------if dflag, then put diagonal nonzero at beginning of row
r(j) = ia(i)
ilast = i
GOTO 6
!
!------put (off-diagonal) nonzero in last unused location in row
5       q(i) = q(i) - 1
r(j) = q(i)
!
6     j = j-1
!
!
!--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
do 8 j=jmin,jmax
7     if (r(j).eq.j)  GOTO 8
k = r(j)
r(j) = r(k)
r(k) = k
jak = ja(k)
ja(k) = ja(j)
ja(j) = jak
ak = a(k)
a(k) = a(j)
a(j) = ak
GOTO 7
8     continue
!
return
end subroutine sro

subroutine md(n, ia,ja, max, v,l, head,last,next, mark, flag)
!***********************************************************************
!  md -- minimum degree algorithm (based on element model)
!***********************************************************************
!
!  description
!
!    md finds a minimum degree ordering of the rows and columns of a
!    general sparse matrix m stored in (ia,ja,a) format.
!    when the structure of m is nonsymmetric, the ordering is that
!    obtained for the symmetric matrix  m + m-transpose.
!
!
!  additional parameters
!
!    max  - declared dimension of the one-dimensional arrays v and l.
!           max must be at least  n+2k,  where k is the number of
!           nonzeroes in the strict upper triangle of m + m-transpose
!
!    v    - integer one-dimensional work array.  dimension = max
!
!    l    - integer one-dimensional work array.  dimension = max
!
!    head - integer one-dimensional work array.  dimension = n
!
!    last - integer one-dimensional array used to return the permutation
!           of the rows and columns of m corresponding to the minimum
!           degree ordering.  dimension = n
!
!    next - integer one-dimensional array used to return the inverse of
!           the permutation returned in last.  dimension = n
!
!    mark - integer one-dimensional work array (may be the same as v).
!           dimension = n
!
!    flag - integer error flag.  values and their meanings are -
!             0     no errors detected
!             9n+k  insufficient storage in md
!
!
!  definitions of internal parameters
!
!    ---------+---------------------------------------------------------
!    v(s)     - value field of list entry
!    ---------+---------------------------------------------------------
!    l(s)     - link field of list entry  (0 =) end of list)
!    ---------+---------------------------------------------------------
!    l(vi)    - pointer to element list of uneliminated vertex vi
!    ---------+---------------------------------------------------------
!    l(ej)    - pointer to boundary list of active element ej
!    ---------+---------------------------------------------------------
!    head(d)  - vj =) vj head of d-list d
!             -  0 =) no vertex in d-list d
!
!
!             -                  vi uneliminated vertex
!             -          vi in ek           -       vi not in ek
!    ---------+-----------------------------+---------------------------
!    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
!             -                             -  0 =) vi tail of d-list
!    ---------+-----------------------------+---------------------------
!    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
!             --vk =) compute degree        - vj =) vj last in d-list
!             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
!             -  0 =) do not compute degree -
!    ---------+-----------------------------+---------------------------
!    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
!
!
!             -                   vi eliminated vertex
!             -      ei active element      -           otherwise
!    ---------+-----------------------------+---------------------------
!    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
!             -       to be eliminated      -       to be eliminated
!    ---------+-----------------------------+---------------------------
!    last(vi) -  m =) size of ei = m        - undefined
!    ---------+-----------------------------+---------------------------
!    mark(vi) - -m =) overlap count of ei   - undefined
!             -       with ek = m           -
!             - otherwise nonnegative tag   -
!             -       .lt. mark(vk)         -
!
!-----------------------------------------------------------------------
!
integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),mark(*),  flag,  tag, dmin, vk,ek, tail
integer n
equivalence  (vk,ek)
!
!----initialization
tag = 0
call  mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
if (flag.ne.0)  return
!
k = 0
dmin = 1
!
!----while  k .lt. n  do
1  if (k.ge.n)  GOTO 4
!
!------search for vertex of minimum degree
2    if (head(dmin).gt.0)  GOTO 3
dmin = dmin + 1
GOTO 2
!
!------remove vertex vk of minimum degree from degree list
3    vk = head(dmin)
head(dmin) = next(vk)
if (head(dmin).gt.0)  last(head(dmin)) = -dmin
!
!------number vertex vk, adjust tag, and tag vk
k = k+1
next(vk) = -k
last(ek) = dmin - 1
tag = tag + last(ek)
mark(vk) = tag
!
!------form element ek from uneliminated neighbors of vk
call  mdm(vk,tail, v,l, last,next, mark)
!
!------purge inactive elements and do mass elimination
call  mdp(k,ek,tail, v,l, head,last,next, mark)
!
!------update degrees of uneliminated vertices in ek
call  mdu(ek,dmin, v,l, head,last,next, mark)
!
GOTO 1
!
!----generate inverse permutation from permutation
4  do 5 k=1,n
next(k) = -next(k)
5    last(next(k)) = k
!
return
end subroutine md

subroutine mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
!***********************************************************************
!  mdi -- initialization
!***********************************************************************
integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*), mark(*), tag,  flag,  sfs, vi,dvi, vj
integer n
!
!----initialize degrees, element lists, and degree lists
do 1 vi=1,n
mark(vi) = 1
l(vi) = 0
1    head(vi) = 0
sfs = n+1
!
!----create nonzero structure
!----for each nonzero entry a(vi,vj)
do 6 vi=1,n
jmin = ia(vi)
jmax = ia(vi+1) - 1
if (jmin.gt.jmax)  GOTO 6
do 5 j=jmin,jmax
vj = ja(j)
if (vj-vi) 2, 5, 4
!
!------if a(vi,vj) is in strict lower triangle
!------check for previous occurrence of a(vj,vi)
2      lvk = vi
kmax = mark(vi) - 1
if (kmax .eq. 0) GOTO 4
do 3 k=1,kmax
lvk = l(lvk)
if (v(lvk).eq.vj) GOTO 5
3        continue
!----for unentered entries a(vi,vj)
4        if (sfs.ge.max)  GOTO 101
!
!------enter vj in element list for vi
mark(vi) = mark(vi) + 1
v(sfs) = vj
l(sfs) = l(vi)
l(vi) = sfs
sfs = sfs+1
!
!------enter vi in element list for vj
mark(vj) = mark(vj) + 1
v(sfs) = vi
l(sfs) = l(vj)
l(vj) = sfs
sfs = sfs+1
5      continue
6    continue
!
!----create degree lists and initialize mark vector
do 7 vi=1,n
dvi = mark(vi)
next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
nextvi = next(vi)
if (nextvi.gt.0)  last(nextvi) = vi
7    mark(vi) = tag
!
return
!
! ** error-  insufficient storage
101  flag = 9*n + vi
return
end subroutine mdi

subroutine mdm(vk,tail, v,l, last,next, mark)
!***********************************************************************
!  mdm -- form element from uneliminated neighbors of vk
!***********************************************************************
integer  vk, tail, v(*), l(*), last(*), next(*), mark(*),tag, s,ls,vs,es, b,lb,vb, blp,blpmax
equivalence  (vs, es)
!
!----initialize tag and list of uneliminated neighbors
tag = mark(vk)
tail = vk
!
!----for each vertex/element vs/es in element list of vk
ls = l(vk)
1  s = ls
if (s.eq.0)  GOTO 5
ls = l(s)
vs = v(s)
if (next(vs).lt.0)  GOTO 2
!
!------if vs is uneliminated vertex, then tag and append to list of
!------uneliminated neighbors
mark(vs) = tag
l(tail) = s
tail = s
GOTO 4
!
!------if es is active element, then ...
!--------for each vertex vb in boundary list of element es
2      lb = l(es)
blpmax = last(es)
do 3 blp=1,blpmax
b = lb
lb = l(b)
vb = v(b)
!
!----------if vb is untagged vertex, then tag and append to list of
!----------uneliminated neighbors
if (mark(vb).ge.tag)  GOTO 3
mark(vb) = tag
l(tail) = b
tail = b
3        continue
!
!--------mark es inactive
mark(es) = tag
!
4    GOTO 1
!
!----terminate list of uneliminated neighbors
5  l(tail) = 0
!
return
end subroutine mdm

subroutine mdp(k,ek,tail, v,l, head,last,next, mark)
!***********************************************************************
!  mdp -- purge inactive elements and do mass elimination
!***********************************************************************
integer  ek, tail,  v(*), l(*),  head(*), last(*), next(*),mark(*),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax
!
!----initialize tag
tag = mark(ek)
!
!----for each vertex vi in ek
li = ek
ilpmax = last(ek)
if (ilpmax.le.0)  GOTO 12
do 11 ilp=1,ilpmax
i = li
li = l(i)
vi = v(li)
!
!------remove vi from degree list
if (last(vi).eq.0)  GOTO 3
if (last(vi).gt.0)  GOTO 1
head(-last(vi)) = next(vi)
GOTO 2
1        next(last(vi)) = next(vi)
2      if (next(vi).gt.0)  last(next(vi)) = last(vi)
!
!------remove inactive items from element list of vi
3    ls = vi
4    s = ls
ls = l(s)
if (ls.eq.0)  GOTO 6
es = v(ls)
if (mark(es).lt.tag)  GOTO 5
free = ls
l(s) = l(ls)
ls = s
5      GOTO 4
!
!------if vi is interior vertex, then remove from list and eliminate
6    lvi = l(vi)
if (lvi.ne.0)  GOTO 7
l(i) = l(li)
li = i
!
k = k+1
next(vi) = -k
last(ek) = last(ek) - 1
GOTO 11
!
!------else ...
!--------classify vertex vi
7      if (l(lvi).ne.0)  GOTO 9
evi = v(lvi)
if (next(evi).ge.0)  GOTO 9
if (mark(evi).lt.0)  GOTO 8
!
!----------if vi is prototype vertex, then mark as such, initialize
!----------overlap count for corresponding element, and move vi to end
!----------of boundary list
last(vi) = evi
mark(evi) = -1
l(tail) = li
tail = li
l(i) = l(li)
li = i
GOTO 10
!
!----------else if vi is duplicate vertex, then mark as such and adjust
!----------overlap count for corresponding element
8            last(vi) = 0
mark(evi) = mark(evi) - 1
GOTO 10
!
!----------else mark vi to compute degree
9            last(vi) = -ek
!
!--------insert ek in element list of vi
10      v(free) = ek
l(free) = l(vi)
l(vi) = free
11    continue
!
!----terminate boundary list
12  l(tail) = 0
!
return
end subroutine mdp

subroutine mdu(ek,dmin, v,l, head,last,next, mark)
!***********************************************************************
!  mdu -- update degrees of uneliminated vertices in ek
!***********************************************************************
integer  ek, dmin,  v(*), l(*),  head(*), last(*), next(*),mark(*),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,blp,blpmax
equivalence  (vs, es)
!
!----initialize tag
tag = mark(ek) - last(ek)
!
!----for each vertex vi in ek
i = ek
ilpmax = last(ek)
if (ilpmax.le.0)  GOTO 11
do 10 ilp=1,ilpmax
i = l(i)
vi = v(i)
if (last(vi))  1, 10, 8
!
!------if vi neither prototype nor duplicate vertex, then merge elements
!------to compute degree
1      tag = tag + 1
dvi = last(ek)
!
!--------for each vertex/element vs/es in element list of vi
s = l(vi)
2      s = l(s)
if (s.eq.0)  GOTO 9
vs = v(s)
if (next(vs).lt.0)  GOTO 3
!
!----------if vs is uneliminated vertex, then tag and adjust degree
mark(vs) = tag
dvi = dvi + 1
GOTO 5
!
!----------if es is active element, then expand
!------------check for outmatched vertex
3          if (mark(es).lt.0)  GOTO 6
!
!------------for each vertex vb in es
b = es
blpmax = last(es)
do 4 blp=1,blpmax
b = l(b)
vb = v(b)
!
!--------------if vb is untagged, then tag and adjust degree
if (mark(vb).ge.tag)  GOTO 4
mark(vb) = tag
dvi = dvi + 1
4            continue
!
5        GOTO 2
!
!------else if vi is outmatched vertex, then adjust overlaps but do not
!------compute degree
6      last(vi) = 0
mark(es) = mark(es) - 1
7      s = l(s)
if (s.eq.0)  GOTO 10
es = v(s)
if (mark(es).lt.0)  mark(es) = mark(es) - 1
GOTO 7
!
!------else if vi is prototype vertex, then calculate degree by
!------inclusion/exclusion and reset overlap count
8      evi = last(vi)
dvi = last(ek) + last(evi) + mark(evi)
mark(evi) = 0
!
!------insert vi in appropriate degree list
9    next(vi) = head(dvi)
head(dvi) = vi
last(vi) = -dvi
if (next(vi).gt.0)  last(next(vi)) = vi
if (dvi.lt.dmin)  dmin = dvi
!
10    continue
!
11  return
end subroutine mdu

SUBROUTINE ADJLR (N, ISP, LDIF)
INTEGER N, ISP, LDIF
DIMENSION ISP(*)
!-----------------------------------------------------------------------
! This routine computes an adjustment, LDIF, to the required
! integer storage space in IWK (sparse matrix work space).
! It is called only if the word length ratio is LRAT = 1.
! This is to account for the possibility that the symbolic LU phase
! may require more storage than the numerical LU and solution phases.
!-----------------------------------------------------------------------
INTEGER IP, JLMAX, JUMAX, LNFC, LSFC, NZLU
!
IP = 2*N + 1
! Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
JLMAX = ISP(IP)
JUMAX = ISP(IP+IP)
! NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
NZLU = ISP(N+1) - ISP(1) + ISP(IP+N+1) - ISP(IP+1)
LSFC = 12*N + 3 + 2*MAX(JLMAX,JUMAX)
LNFC = 9*N + 2 + JLMAX + JUMAX + NZLU
LDIF = MAX(0, LSFC - LNFC)
RETURN
!----------------------- End of Subroutine ADJLR -----------------------
END subroutine adjlr

SUBROUTINE CNTNZU (N, IA, JA, NZSUT)
INTEGER N, IA, JA, NZSUT
DIMENSION IA(*), JA(*)
!-----------------------------------------------------------------------
! This routine counts the number of nonzero elements in the strict
! upper triangle of the matrix M + M(transpose), where the sparsity
! structure of M is given by pointer arrays IA and JA.
! This is needed to compute the storage requirements for the
! sparse matrix reordering operation in ODRV.
!-----------------------------------------------------------------------
INTEGER II, JJ, J, JMIN, JMAX, K, KMIN, KMAX, NUM
!
NUM = 0
DO 50 II = 1,N
JMIN = IA(II)
JMAX = IA(II+1) - 1
IF (JMIN .GT. JMAX) GOTO 50
DO 40 J = JMIN,JMAX
IF (JA(J) - II) 10, 40, 30
10       JJ =JA(J)
KMIN = IA(JJ)
KMAX = IA(JJ+1) - 1
IF (KMIN .GT. KMAX) GOTO 30
DO 20 K = KMIN,KMAX
IF (JA(K) .EQ. II) GOTO 40
20         CONTINUE
30       NUM = NUM + 1
40       CONTINUE
50     CONTINUE
NZSUT = NUM
RETURN
!----------------------- End of Subroutine CNTNZU ----------------------
END subroutine CNTNZU

SUBROUTINE DUMSUM(A,B,C)
!     Routine to force normal storing of A + B, for DUMACH.
DOUBLE PRECISION A, B, C
C = A + B
RETURN
END subroutine dumsum

SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
!***BEGIN PROLOGUE  DINTDY
!***SUBSIDIARY
!***PURPOSE  Interpolate solution derivatives.
!***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DINTDY computes interpolated values of the K-th derivative of the
!  dependent variable vector y, and stores it in DKY.  This routine
!  is called within the package with K = 0 and T = TOUT, but may
!  also be called by the user for any K up to the current order.
!  (See detailed instructions in the usage documentation.)
!
!  The computed values in DKY are gotten by interpolation using the
!  Nordsieck history array YH.  This array corresponds uniquely to a
!  vector-valued polynomial of degree NQCUR or less, and DKY is set
!  to the K-th derivative of this polynomial at T.
!  The formula for DKY is:
!               q
!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
!              j=K
!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
!  communicated by COMMON.  The above sum is done in reverse order.
!  IFLAG is returned negative if either K or T is out of bounds.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  XERRWD
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   050427  Corrected roundoff decrement in TP. (ACH)
!***END PROLOGUE  DINTDY
!**End
INTEGER K, NYH, IFLAG
DOUBLE PRECISION T, YH, DKY
DIMENSION YH(NYH,*), DKY(*)
INTEGER IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF
INTEGER LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
DOUBLE PRECISION ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(6), IOWNS(6),&
ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,&
MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
DOUBLE PRECISION C, R, S, TP
CHARACTER*80 MSG
!
!***FIRST EXECUTABLE STATEMENT  DINTDY
IFLAG = 0
IF (K .LT. 0 .OR. K .GT. NQ) GOTO 80
TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
IF ((T-TP)*(T-TN) .GT. 0.0D0) GOTO 90
!
S = (T - TN)/H
IC = 1
IF (K .EQ. 0) GOTO 15
JJ1 = L - K
DO 10 JJ = JJ1,NQ
10     IC = IC*JJ
15   C = IC
DO 20 I = 1,N
20     DKY(I) = C*YH(I,L)
IF (K .EQ. NQ) GOTO 55
JB2 = NQ - K
DO 50 JB = 1,JB2
J = NQ - JB
JP1 = J + 1
IC = 1
IF (K .EQ. 0) GOTO 35
JJ1 = JP1 - K
DO 30 JJ = JJ1,J
30       IC = IC*JJ
35     C = IC
DO 40 I = 1,N
40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
50     CONTINUE
IF (K .EQ. 0) RETURN
55   R = H**(-K)
DO 60 I = 1,N
60     DKY(I) = R*DKY(I)
RETURN
!
80   MSG = 'DINTDY-  K (=I1) illegal      '
CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
IFLAG = -1
RETURN
90   MSG = 'DINTDY-  T (=R1) illegal      '
CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
IFLAG = -2
RETURN
!----------------------- END OF SUBROUTINE DINTDY ----------------------
END subroutine dintdy

SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.
!
DOUBLE PRECISION R1, R2
INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
CHARACTER*(*) MSG
!
!  Declare local variables.
!
INTEGER LUNIT, MESFLG
!
!  Get logical unit number and message print flag.
!
!***FIRST EXECUTABLE STATEMENT  XERRWD
LUNIT = IXSAV (1, 0, .FALSE.)
MESFLG = IXSAV (2, 0, .FALSE.)
IF (MESFLG .EQ. 0) GOTO 100
!
!  Write the message.
!
WRITE (LUNIT,10)  MSG
10   FORMAT(1X,A)
IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
20   FORMAT(6X,'In above message,  I1 =',I10)
IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
40   FORMAT(6X,'In above message,  R1 =',D21.13)
IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
!
!  Abort the run if LEVEL = 2.
!
100  IF (LEVEL .NE. 2) RETURN
STOP
!----------------------- End of Subroutine XERRWD ----------------------
END subroutine XERRWD

INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!    LUNIT, the logical unit number to which messages are printed, and
!    MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!   LUNIT  = Logical unit number for messages.  The default is obtained
!            by a call to IUMACH (may be machine-dependent).
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  On input..
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!    ISET   = Logical flag to indicate whether to read or write.
!             If ISET = .TRUE., the parameter will be given
!             the value IVALUE.  If ISET = .FALSE., the parameter
!             will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!    IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Modified prologue to SLATEC format. (FNF)
!   930915  Added IUMACH call to get default output unit.  (ACH)
!   930922  Minor cosmetic changes. (FNF)
!   010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
LOGICAL ISET
INTEGER IPAR, IVALUE
!-----------------------------------------------------------------------
INTEGER LUNIT, MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
SAVE LUNIT, MESFLG
DATA LUNIT/-1/, MESFLG/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
IF (IPAR .EQ. 1) THEN
IF (LUNIT .EQ. -1) LUNIT = IUMACH()
IXSAV = LUNIT
IF (ISET) LUNIT = IVALUE
ENDIF
!
IF (IPAR .EQ. 2) THEN
IXSAV = MESFLG
IF (ISET) MESFLG = IVALUE
ENDIF
!
RETURN
!----------------------- End of Function IXSAV -------------------------
END function ixsav

DOUBLE PRECISION FUNCTION DUMACH ()
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
DOUBLE PRECISION U, COMP
!***FIRST EXECUTABLE STATEMENT  DUMACH
U = 1.0D0
10   U = U*0.5D0
CALL DUMSUM(1.0D0, U, COMP)
IF (COMP .NE. 1.0D0) GOTO 10
DUMACH = U*2.0D0
RETURN
!----------------------- End of Function DUMACH ------------------------
END function dumach

INTEGER FUNCTION IUMACH()
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
IUMACH = 6
!
RETURN
!----------------------- End of Function IUMACH ------------------------
END function iumach

end module ode_solve
