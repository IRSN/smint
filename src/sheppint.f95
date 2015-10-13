!******************************************************************************
! M-dimensional interpolation with no gradient.
!
! Written by Yves Deville, not part of SHEPPACK.
!
! This code is directly inspired by the SHEPPACK driver "qshepmd_test.f95"
! See the src/905A direzctory in the package for the original f95 source 
!
!******************************************************************************

subroutine SHEPQM(M, N, X, F, NNEW, XNEW, FNEW, NQ, NW, NR, RMAX, IER, &
     NTT, A, DX, XMIN, RSQ, LCELL, LNEXT, WS, IW)
  
  implicit NONE
 
  integer M, N, NNEW, NQ, NW, NR, IER, IER1, INEW, L
  integer NTT
  
  double precision X(M, N), F(N), XNEW(M, NNEW), XNEW1(M, NNEW), FNEW(NNEW)
  double precision RMAX
  double precision QSMVAL

  ! new args
  double precision A(N, NTT), DX(M), XMIN(M), RSQ(N), WS(NTT*NTT), IW(M, 5)
  integer LCELL(NR**M), LNEXT(N) 

  ! first call QSHEMP. This gives the value of RMAX 
  ! write ( *, * ) 'RMAX before', RMAX
  call QSHEPM(M, N, X, F, NQ, NW, NR, RMAX, IER, &
              NTT, A, DX, XMIN, RSQ, LCELL, LNEXT, WS, IW)
  ! write ( *, * ) 'RMAX after', RMAX

  if (IER /= 0) then
!!$     write ( *, * ) ' '
!!$     write ( *, * ) 'QSHEPM'
!!$     write ( *, * ) 'Error return from QSHEPM, IER = ', IER
     return
  end if
  
  ! copy vector and eval
  do INEW = 1, NNEW
     do L = 1, M
        XNEW1(L, 1) = XNEW(L, INEW)
!        write ( *, * ) L, XNEW1(L, 1)
     end do
     FNEW(INEW) = QSMVAL(M, N, XNEW1, X, F, NR, RMAX, &
          NTT, A, DX, XMIN, RSQ, LCELL, LNEXT, WS, IW)
!     write ( *, * ) FNEW(INEW)
  end do
  
end subroutine SHEPQM

!******************************************************************************
! 2D interpolation with no gradient.
!
! Written by Yves Deville, not part of SHEPPACK.
!
! This code is directly inspired by the SHEPPACK driver "qshep2d_test.f95"
! See the src/905A directory in the package for the original f95 source 
!
!******************************************************************************

subroutine SHEPQ2(N, X, Y, F, NNEW, XNEW, YNEW, FNEW, &
     NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
     DX, DY, RMAX, RSQ, A, IER)
  
  implicit NONE
  
  integer N, NNEW, NQ, NW, NR, IER, INEW
  
  double precision X(N), Y(N), F(N), RSQ(N)
  double precision XNEW(NNEW), YNEW(NNEW), FNEW(NNEW)
  double precision RMAX, XNEW1, YNEW1
  double precision QS2VAL

  ! new args
  double precision A(5,N), DX, DY, XMIN, YMIN 
  integer LCELL(NR, NR), LNEXT(N) 

  ! first call QSHEMP. This gives the value of RMAX 
  ! write ( *, * ) 'RMAX before', RMAX
  call QSHEP2(N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
              DX, DY, RMAX, RSQ, A, IER)
  ! write ( *, * ) 'RMAX after', RMAX
  
  do INEW = 1, NNEW
     XNEW1 = XNEW(INEW)
     YNEW1 = YNEW(INEW)
     FNEW(INEW) = QS2VAL(XNEW1, YNEW1, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
          YMIN, DX, DY, RMAX, RSQ, A) 
  end do

end subroutine SHEPQ2

!******************************************************************************
! 2D interpolation with gradient computation
!
! written by Yves Deville, not part of SHEPPACK.
!
! This code is directly inspired by the SHEPPACK driver "qshep2d_test.f95"
! See the src/905A direzctory in the package for the original f95 source 
!
!******************************************************************************

subroutine SHEPQ2G(N, X, Y, F, NNEW, XNEW, YNEW, FNEW, &
     NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
     DX, DY, RMAX, RSQ, A, IER)
  
  implicit NONE
  
  integer N, NNEW, NQ, NW, NR, IER, DERIV, INEW
  
  double precision X(N), Y(N), F(N), RSQ(N)
  double precision XNEW(NNEW), YNEW(NNEW), FNEW(NNEW, 3)
  double precision RMAX, XNEW1, YNEW1, Q, QX, QY
  double precision QS2VAL

  ! new args
  double precision A(5,N), DX, DY, XMIN, YMIN 
  integer LCELL(NR, NR), LNEXT(N) 

  ! first call QSHEMP. This gives the value of RMAX 
  ! write ( *, * ) 'RMAX before', RMAX
  call QSHEP2(N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
              DX, DY, RMAX, RSQ, A, IER)
  ! write ( *, * ) 'RMAX after', RMAX

  do INEW = 1, NNEW
     XNEW1 = XNEW(INEW)
     YNEW1 = YNEW(INEW)
     call QS2GRD(XNEW1, YNEW1, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
          YMIN, DX, DY, RMAX, RSQ, A, Q, QX, QY, IER) 
     FNEW(INEW, 1) = Q
     FNEW(INEW, 2) = QX
     FNEW(INEW, 3) = QY
  end do
  
end subroutine SHEPQ2G

!******************************************************************************
! 3D interpolation with no gradient
!
! Written by Yves Deville, not part of SHEPPACK.
!
! This code is directly inspired by the SHEPPACK driver "qshep3d_test.f95"
! See the src/905A directory in the package for the original f95 source 
!
!******************************************************************************

subroutine SHEPQ3(N, X, Y, Z, F, NNEW, XNEW, YNEW, ZNEW, FNEW, &
     NQ, NW, NR, LCELL, LNEXT, XYZMIN, XYZDEL, &
     RMAX, RSQ, A, IER)
  
  implicit NONE
  
  integer N, NNEW, NQ, NW, NR, IER, INEW
  
  double precision X(N), Y(N), Z(N), F(N), RSQ(N)
  double precision XNEW(NNEW), YNEW(NNEW), ZNEW(NNEW), FNEW(NNEW)
  double precision RMAX, XNEW1, YNEW1, ZNEW1
  double precision QS3VAL

  ! new args
  double precision A(9,N), XYZMIN(3), XYZDEL(3) 
  integer LCELL(NR, NR, NR), LNEXT(N) 

  ! first call QSHEMP. This gives the value of RMAX 
  ! write ( *, * ) 'RMAX before', RMAX
  call QSHEP3(N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, XYZMIN, &
              XYZDEL, RMAX, RSQ, A, IER)
!!$  write ( *, * ) 'RMAX after', RMAX
  
  do INEW = 1, NNEW
     XNEW1 = XNEW(INEW)
     YNEW1 = YNEW(INEW)
     ZNEW1 = ZNEW(INEW)
     FNEW(INEW) = QS3VAL(XNEW1, YNEW1, ZNEW1, N, X, Y, Z, F, NR, LCELL, LNEXT, &
          XYZMIN, XYZDEL, RMAX, RSQ, A) 
  end do

end subroutine SHEPQ3

!******************************************************************************
! 3D interpolation with gradient computation
!
! Written by Yves Deville, not part of SHEPPACK.
!
! This code is directly inspired by the SHEPPACK driver "qshep3d_test.f95"
! See the src/905A directory in the package for the original f95 source 
!
!******************************************************************************

subroutine SHEPQ3G(N, X, Y, Z, F, NNEW, XNEW, YNEW, ZNEW, FNEW, &
     NQ, NW, NR, LCELL, LNEXT, XYZMIN, XYZDEL, &
     RMAX, RSQ, A, IER)
  
  implicit NONE
  
  integer N, NNEW, NQ, NW, NR, IER, INEW
  
  double precision X(N), Y(N), Z(N), F(N), RSQ(N)
  double precision XNEW(NNEW), YNEW(NNEW), ZNEW(NNEW), FNEW(NNEW, 4)
  double precision RMAX, XNEW1, YNEW1, ZNEW1, Q, QX, QY, QZ
  double precision QS2VAL

  double precision A(9,N), XYZMIN(3), XYZDEL(3) 
  integer LCELL(NR, NR, NR), LNEXT(N) 

  ! first call QSHEMP. This gives the value of RMAX 
  ! write ( *, * ) 'RMAX before', RMAX
  call QSHEP3(N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, XYZMIN, &
       XYZDEL, RMAX, RSQ, A, IER)
  
  ! write ( *, * ) 'RMAX after', RMAX
  
  do INEW = 1, NNEW
     XNEW1 = XNEW(INEW)
     YNEW1 = YNEW(INEW)
     ZNEW1 = ZNEW(INEW)
     call QS3GRD(XNEW1, YNEW1, ZNEW1, N, X, Y, Z, F, NR, LCELL, LNEXT, & 
          XYZMIN, XYZDEL, RMAX, RSQ, A, Q, QX, QY, QZ, IER) 
     FNEW(INEW, 1) = Q
     FNEW(INEW, 2) = QX
     FNEW(INEW, 3) = QY
     FNEW(INEW, 4) = QZ
  end do
  
end subroutine SHEPQ3G
