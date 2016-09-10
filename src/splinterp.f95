!     DECK SPLINTERP
SUBROUTINE SPLINTERP(LCUBIC, N, X, VAL, LDERIV, TERP)
  !     
  !-----------------------------------------------------------------------
  !     
  !  Purpose:
  !
  !     determination of the TERP interpolation/derivation components using
  !     the "natural" cubic (order 4) spline. The second order derivative at
  !     end-points X(1) and X(N)$ is zero.
  !     
  !  Copyright:
  !
  !     Copyright (C) 2006 Ecole Polytechnique de Montreal
  !     This library is free software; you can redistribute it and/or
  !     modify it under the terms of the GNU Lesser General Public
  !     License as published by the Free Software Foundation; either
  !     version 2.1 of the License, or (at your option) any later version
  !     
  !  Author(s): A. Hebert and Y. Deville (adaptation for cubic spline)
  !     
  !  Parameters: input
  !
  !     LCUBIC 
  !
  !        .TRUE.: natural cubic spline interpolation; 
  !
  !        .FALSE: linear Lagrange interpolation.  
  !
  !     N number of points. MUST BE >= 4.  
  !
  !     X abscissas VAL abscissa of the interpolated point.
  !
  !     LDERIV 
  !   
  !         if 0 the interpolation factors are computed.  
  !
  !         If 1 or 2, the derivative factors of order 1 or 2 are
  !         computed.
  !     
  !  Parameters: output
  !     
  !     TERP    interpolation/derivation components.
  !     
  !-----------------------------------------------------------------------
  !
  !---- 
  !     SUBROUTINE ARGUMENTS
  !----

  implicit NONE
  integer N, LCUBIC, LDERIV
  double precision X(N), VAL, TERP(N)
  
  !
  !     LOCAL VARIABLES
  !
 
  integer I0, IL, IR, I, J
  double precision WK(2, N-2), B(N, N) 
  double precision DX, HP, HM, HL, HR, GAR, PMX, XL, XR
  
  I0 = 0
  if (LCUBIC .LE. 0) then
     if (N .LE. 1) return
  else
     IF (N .LE. 2) return
  end if
  
  do I = 1, N
     TERP(I) = 0.0D00
  end do
  
  !
  !     INTERVAL IDENTIFICATION.
  !
  
  do I = 1, N-1
     if ((VAL .GE. X(I)) .AND. (VAL .LE. X(I + 1))) then
        I0 = I
        DX = X(I0 + 1) - X(I0)
        exit
     end if
  end do
  
  !
  !     LINEAR LAGRANGE POLYNOMIAL.
  ! 
  
  if (LCUBIC .LE. 0) then

     if (LDERIV .GT. 0) then
        TERP(I0) = -1.0D00 / DX
        TERP(I0 + 1) = 1.0D00 / DX
     else
        TERP(I0) = (X(I0 + 1) - VAL) / DX
        TERP(I0 + 1) = 1.0D00 - TERP(I0)
     end if

     return

  end if
  
  !
  !     Initial values for 'B'
  !
  
  do J = 1, N   
     do I = 1, N
        B(I, J) = 0.0 
     end do
  end do
  
  !     Compute the non-zero elements in 'A' and 'B'
  !     main diag of 'A' is in the 1-st row of WK
  !
  
  HP =  X(2) - X(1)
  do I = 2, N-1
     HM = HP
     HP = X(I + 1) - X(I)
     WK(1, I - 1) = (HM + HP) / 3.0D00
     WK(2, I - 1) = HP / 6.0D00
     B(I, I - 1) = 1.0D00 / HM
     B(I, I) = - 1.0D00 / HP - 1.0D00 / HM
     B(I, I + 1) = 1.0D00 / HP
  end do
  
  !     Solve A C = B in 'C' and store 'C' in place of 'B'
  !     FORWARD ELIMINATION.
  !
  PMX = WK(1, 1)
  
  do J = 1, N
     B(2, J) = B(2, J) / PMX
  end do
  
  do I = 2, N-2
     GAR = WK(2, I - 1)
     WK(2, I - 1) = WK(2, I - 1) / PMX
     PMX = WK(1, I) - GAR * WK(2, I - 1)

     do J = 1, N
        B(I + 1, J) = (B(I + 1, J) - GAR * B(I, J)) / PMX
     end do

  end do
  
  !
  !     BACK SUBSTITUTION.
  !
   
   do I = N - 3, 1, -1

      do J = 1, N
         B(I + 1, J) = B(I + 1, J) - WK(2, I) * B(I + 2, J)
      end do

   end do

   !
   ! Compute the interpolation factors.
   !

   IL = I0
   IR = I0 + 1
   XL = X(IL)
   XR = X(IR)
   HL = (VAL - XL) / DX
   HR = (XR - VAL) / DX
   
   IF (LDERIV .EQ. 0) THEN
      
      do J = 1, N
         TERP(J) = (B(IL, J) * HR * (HR * HR - 1.0D00) + &
              B(IR, J) * HL * (HL * HL - 1.0D00)) * &
              DX * DX / 6.0D00
      end do
      
      TERP(IL) = TERP(IL) + HR
      TERP(IR) = TERP(IR) + HL
   
   else if (LDERIV .EQ. 1) then
      
      do J = 1,N
         TERP(J) = (B(IL, J) * (-3.0D00 * HR * HR + 1.0D00) + &
              B(IR, J) * (3.0D00 * HL * HL - 1.0D00)) * &
              DX / 6.0D00
      end do
      
      TERP(IL) = TERP(IL) - 1.0D00   / DX
      TERP(IR) = TERP(IR) + 1.0D00 / DX
   
   else if (LDERIV .EQ. 2) then

      do J = 1,N
         TERP(J) = (B(IL, J) * HR + B(IR, J) * HL)
      end do

   end if
   
   return
   
 end SUBROUTINE SPLINTERP
    
