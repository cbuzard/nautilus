PROGRAM extract_ab

IMPLICIT NONE
INTEGER :: i
INTEGER,PARAMETER :: ntime = 9
INTEGER,PARAMETER :: nsmax = 684
INTEGER,DIMENSION(nsmax+2) :: comp
CHARACTER (len=30) :: filename
CHARACTER (len=30) :: file_output
CHARACTER (len=11), DIMENSION(ntime,nsmax) :: spec
REAL(KIND=8), DIMENSION(ntime,nsmax) :: ab
REAL(KIND=8), DIMENSION(ntime) :: time, temp,dens, tau, zeta
REAL(KIND=8) :: zspace

open(10, file='filename',status='old')
DO i = 1, ntime

  READ(10,*) filename

  OPEN(20, file = filename, form='unformatted')
  READ(20) time(i), zspace, spec(i,:)
  READ(20) temp(i), dens(i), tau(i), zeta(i)
  READ(20) ab(i,:)
  !PRINT*,spec(i,1), ab(i,1), time(i)
  CLOSE(20)

ENDDO
CLOSE(10)

DO i = 1,nsmax+2
   comp(i) = i
ENDDO

PRINT*, 'Entrez nom de fichier de sorti:'
READ*, file_output
PRINT*, 'Sorti vers:', TRIM(file_output)//'.dat'

OPEN(30, file=TRIM(file_output)//'.dat',status='unknown')

WRITE(30,FMT='(I2,9X,I4,685(11x,I5))') comp(:)
WRITE(30,FMT='(A2,10x,A4,685(7x,A9))') 'I', 'TIME', spec(1,:)
WRITE(30,*) ' '
DO i = 1, ntime
WRITE(30,FMT='(I2,5x,ES14.8,685(2x,ES14.8))') i, time(i), ab(i,:)
ENDDO
CLOSE(30)


END PROGRAM extract_ab
