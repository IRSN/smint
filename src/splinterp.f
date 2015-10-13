*     DECK SPLINTERP
      SUBROUTINE SPLINTERP(LCUBIC, N, X, VAL, LDERIV, TERP)
*     
*-----------------------------------------------------------------------
*     
*     Purpose:
*     determination of the TERP interpolation/derivation components using
*     the "natural" cubic (order 4) spline. The second order derivative at
*     end-points X(1) and X(N)$ is zero.
*     
*     Copyright:
*     Copyright (C) 2006 Ecole Polytechnique de Montreal
*     This library is free software; you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation; either
*     version 2.1 of the License, or (at your option) any later version
*     
*     Author(s): A. Hebert and Y. Deville (adaptation for cubic spline)
*     
*     Parameters: input
*     LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*     Lagrange interpolation.
*     N       number of points. MUST BE >= 4.
*     X       abscissas
*     VAL     abscissa of the interpolated point.
*     LDERIV  if 0 the interpolation factors are computed. If 1 or 2,
*     the derivative factors of order 1 or 2 are computed.
*     
*     Parameters: output
*     TERP    interpolation/derivation components.
*     
*-----------------------------------------------------------------------
*
*---- 
*     SUBROUTINE ARGUMENTS
*----
      IMPLICIT NONE
      INTEGER N, LCUBIC, LDERIV, I0, IL, IR, I, J
      double precision X(N), VAL, TERP(N), DX, HP, HM, HL, HR, GAR
      double precision PMX, XL, XR
*---  
*     LOCAL VARIABLES
*---  
      CHARACTER HSMG*131
      double precision WK(2, N-2), B(N, N) 
  
*     
      I0 = 0
      IF(LCUBIC .LE. 0) THEN
         IF (N .LE. 1) RETURN
      ELSE 
         IF (N .LE. 2) RETURN
      ENDIF

      DO 10 I = 1, N
         TERP(I) = 0.0D00
 10   CONTINUE

*---- 
*     INTERVAL IDENTIFICATION.
*---- 
      DO 20 I = 1, N-1
         IF ((VAL .GE. X(I)) .AND. (VAL .LE. X(I + 1))) THEN
            I0 = I
            GO TO 30
         ENDIF
 20   CONTINUE
c$$$      WRITE(HSMG,'(35HALTERP: UNABLE TO INTERPOLATE (VAL=,1P,E11.4,
c$$$     18H LIMITS=,E11.4,2H, ,E11.4,2H).)') VAL, X(1), X(N)
      return
 30   DX = X(I0 + 1) - X(I0)
*---- 
*     LINEAR LAGRANGE POLYNOMIAL.
*---- 
      IF(LCUBIC .LE. 0) THEN
         IF (LDERIV .GT. 0) THEN
            TERP(I0) = -1.0D00 / DX
            TERP(I0 + 1) = 1.0D00 / DX
         ELSE
            TERP(I0) = (X(I0 + 1) - VAL) / DX
            TERP(I0 + 1) = 1.0D00 - TERP(I0)
         ENDIF
         RETURN
      ENDIF
*-----
*     Initial values for 'B'
*---- 
      DO 41 J = 1, N   
         DO 42 I = 1, N
            B(I, J) = 0.0 
 42      CONTINUE
 41   CONTINUE      
*-----
*     Compute the non-zero elements in 'A' and 'B'
*     main diag of 'A' is in the 1-st row of WK
*---- 
      HP =  X(2) - X(1)
      DO 43 I = 2, N-1
         HM = HP
         HP = X(I + 1) - X(I)
         WK(1, I - 1) = (HM + HP) / 3.0D00
         WK(2, I - 1) = HP / 6.0D00
         B(I, I - 1) = 1.0D00 / HM
         B(I, I) = - 1.0D00 / HP - 1.0D00 / HM
         B(I, I + 1) = 1.0D00 / HP
 43   CONTINUE
*----
*     Solve A C = B in 'C' and store 'C' in place of 'B'
*     FORWARD ELIMINATION.
*----
      PMX = WK(1, 1)
      DO 44 J = 1, N
         B(2, J) = B(2, J) / PMX
 44   CONTINUE
      DO 45 I = 2, N-2
         GAR = WK(2, I - 1)
         WK(2, I - 1) = WK(2, I - 1) / PMX
         PMX = WK(1, I) - GAR * WK(2, I - 1)
         DO 46 J = 1, N
            B(I + 1, J) = (B(I + 1, J) - GAR * B(I, J)) / PMX
 46      CONTINUE
 45   CONTINUE
*---- 
*     BACK SUBSTITUTION.
*----
      DO 47 I = N-3, 1, -1
         DO 48 J = 1, N
            B(I + 1, J) = B(I + 1, J) - WK(2, I) * B(I + 2, J)
 48      CONTINUE
 47   CONTINUE
*---- 
*  COMPUTE THE INTERPOLATION FACTORS.
*----
      IL = I0
      IR = I0 + 1
      XL = X(IL)
      XR = X(IR)
      HL = (VAL - XL) / DX
      HR = (XR - VAL) / DX
      IF (LDERIV .EQ. 0) THEN
         DO 49 J = 1, N
            TERP(J) = B(IL, J) * HR * (HR * HR - 1.0D00)  
            TERP(J) = TERP(J) + B(IR, J) * HL * (HL * HL - 1.0D00) 
            TERP(J) = TERP(J) * DX * DX / 6.0D00
 49      CONTINUE
         TERP(IL) = TERP(IL) + HR
         TERP(IR) = TERP(IR) + HL
      ELSE IF (LDERIV .EQ. 1) THEN
         DO 50 J = 1,N
            TERP(J) = B(IL, J) * (-3.0D00 * HR * HR + 1.0D00) 
            TERP(J) = TERP(J) + B(IR, J) * (3.0D00 * HL * HL - 1.0D00) 
            TERP(J) = TERP(J) * DX / 6.0D00
 50      CONTINUE
         TERP(IL) = TERP(IL) - 1.0D00   / DX
         TERP(IR) = TERP(IR) + 1.0D00 / DX
      ELSE IF (LDERIV .EQ. 2) THEN
         DO 51 J = 1,N
            TERP(J) = ( B(IL, J) * HR + B(IR, J) * HL )
 51      CONTINUE
      ENDIF

      RETURN
*     
      END
      
