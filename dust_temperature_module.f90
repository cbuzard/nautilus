!******************************************************************************
! MODULE: dust_temperature_module
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines needed to compute consistently
!! the dust temperature. Only dust_temp is public.
!!\n\n The grain temperature is only a function of UV flux and visual extinction
!!. If those parameters are constant, you only need to calculate the grain 
!! temperature once, during initialisation.
!
!******************************************************************************

module dust_temperature_module

use iso_fortran_env, only : error_unit
use global_variables

implicit none

integer, parameter :: nb_wavelengths=500
real(double_precision), dimension(nb_wavelengths) :: wavelength = [& !< Wavelength [Angstrom]
9.118d+02, 9.433d+02, 9.759d+02, 1.010d+03, 1.045d+03, 1.081d+03, 1.118d+03, 1.157d+03, 1.197d+03, 1.238d+03, 1.281d+03, &
1.325d+03, 1.371d+03, 1.418d+03, 1.468d+03, 1.518d+03, 1.571d+03, 1.625d+03, 1.681d+03, 1.739d+03, 1.800d+03, 1.862d+03, &
1.926d+03, 1.993d+03, 2.062d+03, 2.133d+03, 2.207d+03, 2.283d+03, 2.362d+03, 2.444d+03, 2.528d+03, 2.616d+03, 2.706d+03, &
2.800d+03, 2.897d+03, 2.997d+03, 3.100d+03, 3.208d+03, 3.319d+03, 3.433d+03, 3.552d+03, 3.675d+03, 3.802d+03, 3.933d+03, &
4.070d+03, 4.210d+03, 4.356d+03, 4.506d+03, 4.662d+03, 4.824d+03, 4.990d+03, 5.163d+03, 5.342d+03, 5.526d+03, 5.717d+03, &
5.915d+03, 6.120d+03, 6.331d+03, 6.550d+03, 6.777d+03, 7.011d+03, 7.254d+03, 7.504d+03, 7.764d+03, 8.032d+03, 8.310d+03, &
8.598d+03, 8.895d+03, 9.202d+03, 9.521d+03, 9.850d+03, 1.019d+04, 1.054d+04, 1.091d+04, 1.128d+04, 1.168d+04, 1.208d+04, &
1.250d+04, 1.293d+04, 1.338d+04, 1.384d+04, 1.432d+04, 1.481d+04, 1.532d+04, 1.585d+04, 1.640d+04, 1.697d+04, 1.756d+04, &
1.816d+04, 1.879d+04, 1.944d+04, 2.011d+04, 2.081d+04, 2.153d+04, 2.227d+04, 2.304d+04, 2.384d+04, 2.467d+04, 2.552d+04, &
2.640d+04, 2.731d+04, 2.826d+04, 2.924d+04, 3.025d+04, 3.129d+04, 3.238d+04, 3.349d+04, 3.465d+04, 3.585d+04, 3.709d+04, &
3.837d+04, 3.970d+04, 4.107d+04, 4.249d+04, 4.396d+04, 4.548d+04, 4.706d+04, 4.868d+04, 5.037d+04, 5.211d+04, 5.391d+04, &
5.578d+04, 5.771d+04, 5.970d+04, 6.177d+04, 6.390d+04, 6.611d+04, 6.840d+04, 7.076d+04, 7.321d+04, 7.574d+04, 7.836d+04, &
8.107d+04, 8.388d+04, 8.678d+04, 8.978d+04, 9.288d+04, 9.609d+04, 9.942d+04, 1.029d+05, 1.064d+05, 1.101d+05, 1.139d+05, &
1.178d+05, 1.219d+05, 1.261d+05, 1.305d+05, 1.350d+05, 1.397d+05, 1.445d+05, 1.495d+05, 1.547d+05, 1.600d+05, 1.656d+05, &
1.713d+05, 1.772d+05, 1.833d+05, 1.897d+05, 1.962d+05, 2.030d+05, 2.100d+05, 2.173d+05, 2.248d+05, 2.326d+05, 2.406d+05, &
2.490d+05, 2.576d+05, 2.665d+05, 2.757d+05, 2.852d+05, 2.951d+05, 3.053d+05, 3.158d+05, 3.268d+05, 3.381d+05, 3.498d+05, &
3.619d+05, 3.744d+05, 3.873d+05, 4.007d+05, 4.146d+05, 4.289d+05, 4.437d+05, 4.591d+05, 4.750d+05, 4.914d+05, 5.084d+05, &
5.260d+05, 5.441d+05, 5.630d+05, 5.824d+05, 6.026d+05, 6.234d+05, 6.450d+05, 6.673d+05, 6.904d+05, 7.142d+05, 7.389d+05, &
7.645d+05, 7.909d+05, 8.183d+05, 8.466d+05, 8.758d+05, 9.061d+05, 9.375d+05, 9.699d+05, 1.003d+06, 1.038d+06, 1.074d+06, &
1.111d+06, 1.150d+06, 1.189d+06, 1.230d+06, 1.273d+06, 1.317d+06, 1.363d+06, 1.410d+06, 1.458d+06, 1.509d+06, 1.561d+06, &
1.615d+06, 1.671d+06, 1.729d+06, 1.789d+06, 1.850d+06, 1.914d+06, 1.981d+06, 2.049d+06, 2.120d+06, 2.193d+06, 2.269d+06, &
2.348d+06, 2.429d+06, 2.513d+06, 2.600d+06, 2.690d+06, 2.783d+06, 2.879d+06, 2.978d+06, 3.081d+06, 3.188d+06, 3.298d+06, &
3.412d+06, 3.530d+06, 3.652d+06, 3.779d+06, 3.909d+06, 4.044d+06, 4.184d+06, 4.329d+06, 4.479d+06, 4.634d+06, 4.794d+06, &
4.960d+06, 5.131d+06, 5.309d+06, 5.492d+06, 5.682d+06, 5.879d+06, 6.082d+06, 6.292d+06, 6.510d+06, 6.735d+06, 6.968d+06, &
7.209d+06, 7.458d+06, 7.716d+06, 7.983d+06, 8.259d+06, 8.544d+06, 8.840d+06, 9.146d+06, 9.462d+06, 9.789d+06, 1.013d+07, &
1.048d+07, 1.084d+07, 1.122d+07, 1.160d+07, 1.200d+07, 1.242d+07, 1.285d+07, 1.329d+07, 1.375d+07, 1.423d+07, 1.472d+07, &
1.523d+07, 1.576d+07, 1.630d+07, 1.687d+07, 1.745d+07, 1.805d+07, 1.868d+07, 1.932d+07, 1.999d+07, 2.068d+07, 2.140d+07, &
2.214d+07, 2.290d+07, 2.369d+07, 2.451d+07, 2.536d+07, 2.624d+07, 2.715d+07, 2.808d+07, 2.906d+07, 3.006d+07, 3.110d+07, &
3.218d+07, 3.329d+07, 3.444d+07, 3.563d+07, 3.686d+07, 3.814d+07, 3.946d+07, 4.082d+07, 4.223d+07, 4.369d+07, 4.520d+07, &
4.677d+07, 4.838d+07, 5.006d+07, 5.179d+07, 5.358d+07, 5.543d+07, 5.735d+07, 5.933d+07, 6.138d+07, 6.351d+07, 6.570d+07, &
6.798d+07, 7.033d+07, 7.276d+07, 7.528d+07, 7.788d+07, 8.057d+07, 8.336d+07, 8.624d+07, 8.922d+07, 9.231d+07, 9.550d+07, &
9.880d+07, 1.022d+08, 1.058d+08, 1.094d+08, 1.132d+08, 1.171d+08, 1.212d+08, 1.254d+08, 1.297d+08, 1.342d+08, 1.388d+08, &
1.436d+08, 1.486d+08, 1.537d+08, 1.590d+08, 1.645d+08, 1.702d+08, 1.761d+08, 1.822d+08, 1.885d+08, 1.950d+08, 2.018d+08, &
2.087d+08, 2.160d+08, 2.234d+08, 2.312d+08, 2.391d+08, 2.474d+08, 2.560d+08, 2.648d+08, 2.740d+08, 2.835d+08, 2.933d+08, &
3.034d+08, 3.139d+08, 3.248d+08, 3.360d+08, 3.476d+08, 3.596d+08, 3.721d+08, 3.849d+08, 3.982d+08, 4.120d+08, 4.263d+08, &
4.410d+08, 4.562d+08, 4.720d+08, 4.883d+08, 5.052d+08, 5.227d+08, 5.408d+08, 5.595d+08, 5.788d+08, 5.989d+08, 6.196d+08, &
6.410d+08, 6.632d+08, 6.861d+08, 7.098d+08, 7.344d+08, 7.598d+08, 7.860d+08, 8.132d+08, 8.413d+08, 8.704d+08, 9.005d+08, &
9.317d+08, 9.639d+08, 9.972d+08, 1.032d+09, 1.067d+09, 1.104d+09, 1.143d+09, 1.182d+09, 1.223d+09, 1.265d+09, 1.309d+09, &
1.354d+09, 1.401d+09, 1.449d+09, 1.500d+09, 1.551d+09, 1.605d+09, 1.661d+09, 1.718d+09, 1.777d+09, 1.839d+09, 1.903d+09, &
1.968d+09, 2.036d+09, 2.107d+09, 2.180d+09, 2.255d+09, 2.333d+09, 2.414d+09, 2.497d+09, 2.584d+09, 2.673d+09, 2.765d+09, &
2.861d+09, 2.960d+09, 3.062d+09, 3.168d+09, 3.278d+09, 3.391d+09, 3.508d+09, 3.630d+09, 3.755d+09, 3.885d+09, 4.019d+09, &
4.158d+09, 4.302d+09, 4.451d+09, 4.605d+09, 4.764d+09, 4.929d+09, 5.099d+09, 5.276d+09, 5.458d+09, 5.647d+09, 5.842d+09, &
6.044d+09, 6.253d+09, 6.470d+09, 6.693d+09, 6.925d+09, 7.164d+09, 7.412d+09, 7.668d+09, 7.934d+09, 8.208d+09, 8.492d+09, &
8.785d+09, 9.089d+09, 9.404d+09, 9.729d+09, 1.007d+10, 1.041d+10, 1.077d+10, 1.115d+10, 1.153d+10, 1.193d+10, 1.234d+10, &
1.277d+10, 1.321d+10, 1.367d+10, 1.414d+10, 1.463d+10, 1.514d+10, 1.566d+10, 1.620d+10, 1.676d+10, 1.734d+10, 1.794d+10, &
1.856d+10, 1.920d+10, 1.987d+10, 2.055d+10, 2.126d+10]

real(double_precision), dimension(nb_wavelengths) :: incident_dust_flux = [& !< incident dust flux [erg.cm-2.s-1.A-1.sr-1]
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 5.888d-10, 1.222d-09, 6.861d-10, 4.963d-10, 5.079d-10, 6.569d-10, &
9.016d-10, 8.933d-10, 7.119d-10, 5.954d-10, 5.372d-10, 5.085d-10, 4.942d-10, 4.871d-10, 4.836d-10, 4.801d-10, 4.733d-10, &
4.543d-10, 4.145d-10, 4.165d-10, 4.551d-10, 4.763d-10, 4.863d-10, 4.920d-10, 4.950d-10, 4.973d-10, 4.996d-10, 5.014d-10, &
5.041d-10, 5.087d-10, 5.183d-10, 5.394d-10, 5.960d-10, 8.052d-10, 3.220d-09, 4.763d-09, 8.774d-10, 6.064d-10, 5.304d-10, &
4.990d-10, 4.842d-10, 4.776d-10, 4.766d-10, 4.812d-10, 4.927d-10, 5.141d-10, 5.552d-10, 6.561d-10, 1.143d-09, 1.052d-09, &
1.397d-09, 2.041d-09, 2.426d-09, 1.120d-08, 4.654d-09, 3.403d-09, 3.383d-09, 3.913d-09, 6.347d-09, 1.307d-08, 1.244d-08, &
6.131d-09, 4.507d-09, 4.921d-09, 2.088d-09, 1.457d-09, 1.246d-09, 1.168d-09, 1.193d-09, 1.422d-09, 2.722d-09, 6.300d-09, &
2.909d-09, 2.923d-09, 3.828d-09, 2.029d-09, 1.929d-09, 1.410d-09, 1.226d-09, 1.120d-09, 1.100d-09, 1.157d-09, 1.524d-09, &
1.441d-09, 1.125d-09, 9.238d-10, 8.874d-10, 7.588d-10, 6.986d-10, 6.503d-10, 6.084d-10, 5.721d-10, 5.410d-10, 5.146d-10, &
4.928d-10, 4.747d-10, 4.608d-10, 4.504d-10, 4.427d-10, 4.381d-10, 4.356d-10, 4.344d-10, 4.341d-10, 4.338d-10, 4.338d-10, &
4.335d-10, 4.331d-10, 4.331d-10, 4.334d-10, 4.353d-10, 4.393d-10, 4.456d-10, 4.548d-10, 4.675d-10, 4.838d-10, 5.040d-10, &
5.285d-10, 5.567d-10, 5.886d-10, 6.239d-10, 6.625d-10, 7.034d-10, 7.464d-10, 7.900d-10, 8.346d-10, 8.777d-10, 9.207d-10, &
9.611d-10, 9.986d-10, 1.032d-09, 1.061d-09, 1.086d-09, 1.105d-09, 1.119d-09, 1.127d-09, 1.130d-09, 1.126d-09, 1.118d-09, &
1.104d-09, 1.088d-09, 1.070d-09, 1.050d-09, 1.030d-09, 1.011d-09, 9.866d-10, 9.447d-10, 8.836d-10, 8.174d-10, 7.537d-10, &
6.946d-10, 6.396d-10, 5.874d-10, 5.376d-10, 4.913d-10, 4.467d-10, 4.052d-10, 3.661d-10, 3.299d-10, 2.962d-10, 2.649d-10, &
2.362d-10, 2.099d-10, 1.854d-10, 1.626d-10, 1.421d-10, 1.240d-10, 1.080d-10, 9.379d-11, 8.130d-11, 7.037d-11, 6.074d-11, &
5.237d-11, 4.507d-11, 3.872d-11, 3.324d-11, 2.847d-11, 2.437d-11, 2.084d-11, 1.779d-11, 1.517d-11, 1.294d-11, 1.102d-11, &
9.370d-12, 7.962d-12, 6.764d-12, 5.737d-12, 4.869d-12, 4.126d-12, 3.495d-12, 2.961d-12, 2.505d-12, 2.121d-12, 1.793d-12, &
1.516d-12, 1.280d-12, 1.082d-12, 9.128d-13, 7.698d-13, 6.491d-13, 5.468d-13, 4.600d-13, 3.868d-13, 3.252d-13, 2.734d-13, &
2.296d-13, 1.925d-13, 1.619d-13, 1.359d-13, 1.139d-13, 9.553d-14, 8.005d-14, 6.709d-14, 5.628d-14, 4.718d-14, 3.950d-14, &
3.317d-14, 2.774d-14, 2.323d-14, 1.946d-14, 1.633d-14, 1.367d-14, 1.147d-14, 9.601d-15, 8.044d-15, 6.741d-15, 5.653d-15, &
4.738d-15, 3.974d-15, 3.333d-15, 2.794d-15, 2.344d-15, 1.967d-15, 1.652d-15, 1.386d-15, 1.166d-15, 9.791d-16, 8.229d-16, &
6.913d-16, 5.810d-16, 4.888d-16, 4.112d-16, 3.457d-16, 2.911d-16, 2.451d-16, 2.064d-16, 1.738d-16, 1.465d-16, 1.235d-16, &
1.041d-16, 8.780d-17, 7.406d-17, 6.248d-17, 5.274d-17, 4.450d-17, 3.758d-17, 3.175d-17, 2.683d-17, 2.267d-17, 1.916d-17, &
1.620d-17, 1.371d-17, 1.159d-17, 9.814d-18, 8.305d-18, 7.032d-18, 5.953d-18, 5.045d-18, 4.276d-18, 3.620d-18, 3.066d-18, &
2.596d-18, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, &
0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00]

real(double_precision), dimension(nb_wavelengths) :: Qabs = [&
1.176d+00, 1.176d+00, 1.200d+00, 1.226d+00, 1.226d+00, 1.261d+00, 1.261d+00, 1.282d+00, 1.314d+00, 1.314d+00, 1.336d+00, &
1.336d+00, 1.345d+00, 1.343d+00, 1.343d+00, 1.405d+00, 1.405d+00, 1.420d+00, 1.332d+00, 1.332d+00, 1.151d+00, 1.151d+00, &
7.968d-01, 7.966d-01, 7.173d-01, 4.604d-01, 4.603d-01, 3.486d-01, 3.485d-01, 3.906d-01, 4.159d-01, 4.159d-01, 3.543d-01, &
3.542d-01, 2.850d-01, 2.556d-01, 2.556d-01, 2.548d-01, 2.548d-01, 2.626d-01, 2.584d-01, 2.584d-01, 2.355d-01, 2.355d-01, &
2.036d-01, 2.035d-01, 1.736d-01, 1.499d-01, 1.498d-01, 1.314d-01, 1.314d-01, 1.170d-01, 1.053d-01, 1.053d-01, 9.535d-02, &
9.532d-02, 8.668d-02, 7.903d-02, 7.901d-02, 7.223d-02, 7.221d-02, 6.614d-02, 6.066d-02, 6.064d-02, 5.572d-02, 5.569d-02, &
5.126d-02, 5.123d-02, 4.723d-02, 4.363d-02, 4.360d-02, 4.037d-02, 4.034d-02, 3.743d-02, 3.481d-02, 3.479d-02, 3.241d-02, &
3.239d-02, 3.033d-02, 2.846d-02, 2.843d-02, 2.672d-02, 2.669d-02, 2.511d-02, 2.366d-02, 2.362d-02, 2.233d-02, 2.229d-02, &
2.114d-02, 2.110d-02, 2.006d-02, 1.911d-02, 1.906d-02, 1.819d-02, 1.814d-02, 1.736d-02, 1.666d-02, 1.658d-02, 1.598d-02, &
1.586d-02, 1.540d-02, 2.652d-02, 1.920d-02, 1.378d-02, 1.400d-02, 1.351d-02, 1.361d-02, 1.305d-02, 1.253d-02, 1.259d-02, &
1.211d-02, 1.216d-02, 1.168d-02, 1.122d-02, 1.125d-02, 1.087d-02, 1.090d-02, 1.062d-02, 1.047d-02, 1.048d-02, 1.048d-02, &
1.048d-02, 1.065d-02, 1.090d-02, 1.089d-02, 1.094d-02, 1.094d-02, 1.104d-02, 1.104d-02, 1.412d-02, 3.264d-02, 3.187d-02, &
6.965d-02, 6.806d-02, 1.064d-01, 1.436d-01, 1.421d-01, 1.261d-01, 1.267d-01, 1.037d-01, 7.760d-02, 7.861d-02, 5.998d-02, &
6.071d-02, 4.526d-02, 3.397d-02, 3.440d-02, 2.872d-02, 2.893d-02, 3.232d-02, 3.220d-02, 3.812d-02, 4.472d-02, 4.448d-02, &
4.895d-02, 4.879d-02, 4.856d-02, 4.341d-02, 4.359d-02, 3.793d-02, 3.813d-02, 3.339d-02, 2.931d-02, 2.946d-02, 2.597d-02, &
2.610d-02, 2.294d-02, 2.023d-02, 2.033d-02, 1.798d-02, 1.806d-02, 1.593d-02, 1.601d-02, 1.417d-02, 1.262d-02, 1.268d-02, &
1.132d-02, 1.136d-02, 1.008d-02, 8.930d-03, 8.969d-03, 7.926d-03, 7.963d-03, 7.032d-03, 6.208d-03, 6.236d-03, 5.503d-03, &
5.529d-03, 4.873d-03, 4.294d-03, 4.313d-03, 3.802d-03, 3.820d-03, 3.364d-03, 3.380d-03, 2.979d-03, 2.624d-03, 2.637d-03, &
2.328d-03, 2.339d-03, 2.062d-03, 1.820d-03, 1.829d-03, 1.615d-03, 1.622d-03, 1.434d-03, 1.268d-03, 1.274d-03, 1.125d-03, &
1.131d-03, 1.002d-03, 8.864d-04, 8.903d-04, 7.885d-04, 7.920d-04, 7.016d-04, 7.047d-04, 6.243d-04, 5.534d-04, 5.558d-04, &
4.927d-04, 4.949d-04, 4.389d-04, 3.894d-04, 3.911d-04, 3.471d-04, 3.486d-04, 3.092d-04, 2.744d-04, 2.756d-04, 2.445d-04, &
2.456d-04, 2.180d-04, 1.935d-04, 1.943d-04, 1.726d-04, 1.733d-04, 1.539d-04, 1.546d-04, 1.372d-04, 1.218d-04, 1.224d-04, &
1.087d-04, 1.091d-04, 9.697d-05, 8.614d-05, 8.651d-05, 7.685d-05, 7.718d-05, 6.855d-05, 6.089d-05, 6.115d-05, 5.433d-05, &
5.457d-05, 4.848d-05, 4.869d-05, 4.325d-05, 3.843d-05, 3.859d-05, 3.429d-05, 3.444d-05, 3.060d-05, 2.719d-05, 2.730d-05, &
2.426d-05, 2.436d-05, 2.165d-05, 1.923d-05, 1.931d-05, 1.717d-05, 1.724d-05, 1.532d-05, 1.361d-05, 1.367d-05, 1.242d-05, &
1.160d-05, 1.083d-05, 1.011d-05, 9.443d-06, 8.817d-06, 8.232d-06, 7.686d-06, 7.177d-06, 6.701d-06, 6.257d-06, 5.842d-06, &
5.455d-06, 5.093d-06, 4.756d-06, 4.440d-06, 4.146d-06, 3.871d-06, 3.615d-06, 3.375d-06, 3.151d-06, 2.942d-06, 2.747d-06, &
2.565d-06, 2.395d-06, 2.236d-06, 2.088d-06, 1.950d-06, 1.820d-06, 1.700d-06, 1.587d-06, 1.482d-06, 1.384d-06, 1.292d-06, &
1.206d-06, 1.126d-06, 1.052d-06, 9.819d-07, 9.168d-07, 8.560d-07, 7.993d-07, 7.463d-07, 6.968d-07, 6.506d-07, 6.075d-07, &
5.672d-07, 5.296d-07, 4.945d-07, 4.617d-07, 4.311d-07, 4.025d-07, 3.759d-07, 3.509d-07, 3.277d-07, 3.060d-07, 2.857d-07, &
2.667d-07, 2.491d-07, 2.325d-07, 2.171d-07, 2.027d-07, 1.893d-07, 1.767d-07, 1.650d-07, 1.541d-07, 1.439d-07, 1.343d-07, &
1.254d-07, 1.171d-07, 1.094d-07, 1.021d-07, 9.534d-08, 8.902d-08, 8.312d-08, 7.761d-08, 7.246d-08, 6.766d-08, 6.317d-08, &
5.898d-08, 5.507d-08, 5.142d-08, 4.801d-08, 4.483d-08, 4.186d-08, 3.908d-08, 3.649d-08, 3.407d-08, 3.182d-08, 2.971d-08, &
2.774d-08, 2.590d-08, 2.418d-08, 2.258d-08, 2.108d-08, 1.968d-08, 1.838d-08, 1.716d-08, 1.602d-08, 1.496d-08, 1.397d-08, &
1.304d-08, 1.218d-08, 1.137d-08, 1.062d-08, 9.914d-09, 9.256d-09, 8.643d-09, 8.070d-09, 7.535d-09, 7.035d-09, 6.569d-09, &
6.134d-09, 5.727d-09, 5.347d-09, 4.993d-09, 4.662d-09, 4.353d-09, 4.064d-09, 3.795d-09, 3.543d-09, 3.308d-09, 3.089d-09, &
2.884d-09, 2.693d-09, 2.515d-09, 2.348d-09, 2.192d-09, 2.047d-09, 1.911d-09, 1.785d-09, 1.666d-09, 1.556d-09, 1.453d-09, &
1.356d-09, 1.266d-09, 1.182d-09, 1.104d-09, 1.031d-09, 9.626d-10, 8.987d-10, 8.392d-10, 7.835d-10, 7.316d-10, 6.831d-10, &
6.378d-10, 5.955d-10, 5.561d-10, 5.192d-10, 4.848d-10, 4.526d-10, 4.226d-10, 3.946d-10, 3.685d-10, 3.440d-10, 3.212d-10, &
2.999d-10, 2.800d-10, 2.615d-10, 2.441d-10, 2.280d-10, 2.129d-10, 1.987d-10, 1.856d-10, 1.733d-10, 1.618d-10, 1.511d-10, &
1.410d-10, 1.317d-10, 1.230d-10, 1.148d-10, 1.072d-10, 1.001d-10, 9.346d-11, 8.726d-11, 8.148d-11, 7.608d-11, 7.103d-11, &
6.632d-11, 6.193d-11, 5.782d-11, 5.399d-11, 5.041d-11, 4.707d-11, 4.395d-11, 4.103d-11, 3.831d-11, 3.577d-11, 3.340d-11, &
3.119d-11, 2.912d-11, 2.719d-11, 2.539d-11, 2.371d-11, 2.213d-11, 2.067d-11, 1.930d-11, 1.802d-11, 1.682d-11, 1.571d-11, &
1.467d-11, 1.369d-11, 1.279d-11, 1.194d-11, 1.115d-11, 1.041d-11, 9.718d-12, 9.074d-12, 8.473d-12, 7.911d-12, 7.387d-12, &
6.897d-12, 6.440d-12, 6.013d-12, 5.614d-12, 5.242d-12, 4.894d-12, 4.570d-12, 4.267d-12, 3.984d-12, 3.720d-12, 3.473d-12, &
3.243d-12, 3.028d-12, 2.827d-12, 2.640d-12, 2.465d-12]

real(double_precision), dimension(nb_wavelengths) :: Iinc !< sum of the 4 composant in the incident flux, 3 parts of the spectrum, and dust
real(double_precision) :: rate_abs

private

public :: get_grain_temperature_computed !< Only get_grain_temperature_computed routine will be available outside the module.

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud & Christophe Cossou
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief At a given time, compute the dust temperature using Uv_flux and 
!! Visual_extinction.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_computed(time, gas_temperature, av, grain_temperature)

implicit none

! Inputs
real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
real(double_precision), intent(in) :: av !<[in] visual extinction [mag]

! Outputs
real(double_precision), intent(out) :: grain_temperature !<[out] Steady state dust temperature [K]

! Locals
real(double_precision) :: Tleft, Tright ! Temperature boundary
real(double_precision), parameter :: eps = 1.0D-08 ! Tolerance
real(double_precision) :: fss ! Rq: fss should be equal to zero
integer :: etat ! flag

!----
! Compute the incident radiation filed in term of specific
! intensity (erg cm-2 s-1 Angstrom-1 sr-1) after screening 
! by dust
call compute_Iinc(wavelength, av, Iinc)

!----
! Compute the absorption rate by dust
call compute_abs(wavelength, Iinc, Qabs, rate_abs)

!----
! Find the steady state dust temperature using a 
! dichotomy algorithm.

! Temperature boundary
Tleft = 2.73d00
Tright = 1.0d03

call dicho(Tleft, Tright, eps, grain_temperature, fss, etat)

end subroutine get_grain_temperature_computed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief Compute the incident radiation filed in term of specific
!! intensity (erg cm-2 s-1 Angstrom-1 sr-1) after screening 
!! by dust
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine compute_Iinc(x, av, Iout)

implicit none

! Inputs
real(double_precision), dimension(nb_wavelengths), intent(in) :: x !<[in] Wavelength [Angstrom]
real(double_precision), intent(in) :: av !< Visual extinction [mag]

! Outputs
real(double_precision), dimension(nb_wavelengths), intent(out) :: Iout !<[out] Incident flux after dust screening [erg cm-2 s-1 Angstrom-1 sr-1]

! Locals
integer :: i
real(double_precision) :: Iin

Iout = 0.d0
do i=1, nb_wavelengths
   Iin = max(get_local_UV_flux(x(i)), get_starlight_flux(x(i)), get_CMB_flux(x(i)), incident_dust_flux(i))
   Iout(i) = Iin * 10.d0**(-ext(x(i)) * av / 2.5d0)
enddo


end subroutine compute_Iinc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine compute_abs(x, Iin, Q, dustabs)

implicit none

! Inputs
real(double_precision), dimension(nb_wavelengths), intent(in) :: x
real(double_precision), dimension(nb_wavelengths), intent(in) :: Iin
real(double_precision), dimension(nb_wavelengths), intent(in) :: Q

! Outputs
real(double_precision), intent(out) :: dustabs

! Locals
real(double_precision), dimension(nb_wavelengths) :: Uin
real(double_precision) :: sigmabs
real(double_precision) :: test
real(double_precision), dimension(nb_wavelengths) :: Utest
integer :: i

! --- Convert the specific intensity into specific energy density
Uin = 4.0 * pi * Iin / SPEED_OF_LIGHT

! --- Test: Compute the integrated energy density
Utest=0.d0
do i = 1, nb_wavelengths
   if(x(i).ge.911.0 .AND. x(i).le. 2460.0) Utest(i) = Uin(i)
enddo

call trap_int(nb_wavelengths, x, Utest, test)

call trap_int(nb_wavelengths, x, Uin, test)
!---

! --- Compute the heating rate by absorption of radiation
call trap_int(nb_wavelengths, x, Uin*Q, sigmabs)
dustabs = sigmabs * pi * grain_radius**2 * SPEED_OF_LIGHT

end subroutine compute_abs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief TODO
!
!> @return TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
real(double_precision)  function get_ss_grain(Tss)

implicit none

! Inputs
real(double_precision), intent(in) :: Tss

! Locals
real(double_precision) :: rate_emi
real(double_precision), dimension(nb_wavelengths) :: dustem
real(double_precision) :: sigmaem
integer :: i

! --- Compute the cooling rate by infrared emission
do i = 1, nb_wavelengths
   dustem(i) = get_black_body_flux(wavelength(i), Tss)
enddo

call trap_int(nb_wavelengths, wavelength, dustem*Qabs, sigmaem)
rate_emi = 4.d0 * pi**2 * grain_radius**2 * sigmaem

! --- cooling balance heating
get_ss_grain = rate_emi - rate_abs

end function get_ss_grain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief 1st component
!!\n   Local UV background between 912 and 2460 Angstroms (From Mathis et al. 1983)
!!\n   The fit parameters come from "Physics of the ISM and IGM", Draine, 2011
!!\n      - x in Angstroms
!!\n      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!
!> @return get the local UV flux at the given wavelength
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(double_precision) function  get_local_UV_flux(x)

implicit none

! Inputs
real(double_precision), intent(in) :: x !<[in] Wavelength (Angstrom)

if (x .GE. 1340.0 .AND. x .LT. 2460.0) then
   get_local_UV_flux = 2.373D-14 * (x/1e4)**(-0.6678)
elseif (x .GE. 1100.0 .AND. x .LT. 1340.0) then
   get_local_UV_flux = 6.825D-13*(x/1e4)
elseif (x .GE. 912.0 .AND. x .LT. 1100.0) then
   get_local_UV_flux = 1.287D-9 * (x/1e4)**(4.4172)
else
   get_local_UV_flux = 0.0d00
endif

get_local_UV_flux = get_local_UV_flux * SPEED_OF_LIGHT / (4.d0*pi) / x * UV_FLUX

end function get_local_UV_flux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief 2nd component
!!\n   Starlight emission from near-UV, visible and near-IR (From Mathis et al. 1983)
!!\n   3 diluted Black Body (From Mathis et al. 1983)
!!\n      - x in Angstroms
!!\n      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!
!> @return Starlight emission at a given wavelength
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(double_precision) function  get_starlight_flux(x)

implicit none

! Inputs
real(double_precision), intent(in) :: x !<[in] Wavelength (Angstrom)

get_starlight_flux =  1.0d-14*get_black_body_flux(x, 7500.d0) &
       + 1.0d-13*get_black_body_flux(x, 4000.d0) &
       + 4.0d-13*get_black_body_flux(x, 3000.d0)

end function get_starlight_flux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief 3rd component
!!\n   CMB component at T=2.725 K
!!\n      - x in Angstroms
!!\n      - starlight(x) in erg cm-2 s-1 Angstrom-1 sr-1
!
!> @return Return the Black body spectrum of the CMB at the given wavelength
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(double_precision) function  get_CMB_flux(x)

implicit none

! Inputs
real(double_precision), intent(in) :: x !<[in] Wavelength (Angstrom)

get_CMB_flux = get_black_body_flux(x, 2.725d0)

end function get_CMB_flux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief Black body. x in Angstrom, T in K
!
!
!> @return get_black_body_flux in erg cm-2 s-1 Angstrom-1 sr-1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(double_precision) function  get_black_body_flux(x, T)

implicit none

! Inputs
real(double_precision), intent(in) :: x !<[in] Wavelength  [Angstrom]
real(double_precision), intent(in) :: T !<[in] Temperature [K]

get_black_body_flux = 2.d0 * PLANCK_CONSTANT * SPEED_OF_LIGHT**2 * 1.0d32 / &
       (x**5 * (exp((PLANCK_CONSTANT*SPEED_OF_LIGHT/K_B)*1.0d8 / (x*T)) - 1.d0))

end function get_black_body_flux
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief Fit parameter from Cardelli et al. 1989 (1989ApJ...345..245C)
!! -- The x range is optimized to have a smooth transition between each
!!    components of the extinction curve
!
!> @return TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(double_precision) function  ext(wl)

implicit none

! Inputs
real(double_precision), intent(in) :: wl !< Wavelength [Angstrom]

! Locals
real(double_precision) :: x
real(double_precision) :: a, b, fa, fb, y

! convert x which is in Angstrom to micrometer-1
x = 1.0e4/wl

a = 0.d0
b = 0.d0
! --- Infrared: 0.0 to 1.4 micrometer-1
if(x .lt. 1.4) then
   a =  0.574 * x**1.61
   b = -0.527 * x**1.61
! Optical/NIR: 1.4 to 2.7 micrometer-1
elseif(x .ge. 1.4 .and. x .le. 2.7) then
  y = x - 1.82
  a = 1.00 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + &
      0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
  b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - &
      0.62251*y**5 + 5.30260*y**5 - 2.09002*y**7
! UV: 2.7 to 8.0 micrometer-1
elseif(x .ge. 2.7 .and. x .le. 8.0) then
  if(x .ge. 5.9 .and. x .le. 8.0) then
     fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
     fb =  0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
  elseif(x .le. 5.9) then
     fa = 0.d0
     fb = 0.d0
  endif
  a = 1.752 - 0.316*x - 0.104/((x-4.47)**2 + 0.341) + fa
  b = -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263) + fb
! Far-UV: 8.0 to ...  micrometer-1
elseif(x .ge. 8.0) then
  a = -1.073 - 0.628*(x-8.0) + 0.137*(x-8.0)**2 - 0.070*(x-8.0)**3
  b = 13.670 + 4.257*(x-8.0) - 0.420*(x-8.0)**2 + 0.374*(x-8.0)**3
endif

ext = a + b/3.1

end function ext

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine trap_int(n, x, f, res)

implicit none

! Inputs
integer, intent(in) :: n
real(double_precision), dimension(n), intent(in) :: x, f

! Outputs
real(double_precision), intent(out) :: res

! Locals
integer :: i
real(double_precision) :: h, int

res = 0.d0
int = 0.d0

do i=1, n-1
  h = abs((x(i+1) - x(i))/2.0)
  int = h * (f(i) + f(i+1))
  res = res + int
enddo

end subroutine trap_int

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Maxime Ruaud
!
!> @date may 2014
!
! DESCRIPTION: 
!> @brief TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dicho(xl_in, xr_in, eps_in, xtry, ftry, etat)

implicit none

! Inputs
real(double_precision), intent(in) :: xl_in, xr_in !<[in] Borne gauche et droite initiale
real(double_precision), intent(in) :: eps_in !<[in] Precision

! Outputs
real(double_precision), intent(out) :: xtry, ftry !<[out] Point milieu
integer, intent(out):: etat !<[out] Flag

! Locals
real(double_precision) :: xl, xr ! Borne gauche et droite
real(double_precision) :: fl, fr ! f(xl) et f(xr)
real(double_precision) :: dxrel ! Variation
integer :: i ! Compteur
integer, parameter :: imax = 10000

!--- Initialisation

xl = xl_in
xr = xr_in

fl = get_ss_grain(xl)
fr = get_ss_grain(xr)

!--- consistency test
if (xl >= xr) then
   write(Error_unit,*) "Error: In dust_temperature:dicho"
   write(Error_unit,*) "We need xl < xr but xl=", xl, " and xr=", xr
   stop
endif

if (fl*fr < 0.0d00) then ! The zero is between xl and xr
   etat = 0 ! No solution found yet
   i    = 0
   do while (etat == 0)
      i = i + 1
      xtry = 0.5d00 * (xl + xr)
      ftry = get_ss_grain(xtry)
      if (fl*ftry > 0.0d00) then
         xl = xtry
         fl = ftry
      else
         xr = xtry
         fr = ftry
      endif


      if (xtry /= 0.0d00) then
         dxrel = ABS((xr-xl)/xtry)
      else
         dxrel = ABS(xr-xl)
      endif
      if (dxrel < eps_in .AND. ABS(ftry) < eps_in) then
         etat = 1
      else if (i > imax) then
         etat = 2
      endif
   enddo

else
   etat = 3
endif

if (etat == 3) then
   write(Error_unit, *) " "
   write(Error_unit, *) "--------------------------------------------------------  "
   write(Error_unit, *) " Error in subroutine DICHO of dust_temperature.f90:"
   write(Error_unit, *) " the sign of the function doesn't change in the interval:"
   write(Error_unit, *) " x(left)  :", xl, " f(left)   :", fl
   write(Error_unit, *) " x(right) :", xr, " f(right)  :", fr
   write(Error_unit, *) "--------------------------------------------------------  "
   write(Error_unit, *) " "
   call exit(11)
endif

end subroutine dicho

end module dust_temperature_module