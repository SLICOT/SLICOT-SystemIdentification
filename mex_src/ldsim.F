#include "fintrf.h"
C LDSIM.F - Gateway function for computing the output response of a
C           linear discrete-time system using SLICOT routine TF01MY.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Y(,x)] = ldsim(A,B,C,D,U(,x,ldwork))
C
C Purpose:
C   To compute the output vector sequence y(1), y(2),..., y(t)
C
C        x(k+1) = A x(k) + B u(k)
C        y(k)   = C x(k) + D u(k),
C
C   given an initial state vector x(1), and the input vector sequence
C   u(1), u(2),..., u(t), where y(k) and u(k) are vectors of length p
C   and m, respectively. The input trajectories are given as
C
C            ( u(1)' )
C            ( u(2)' )
C        U = (   :   ),  where ' denotes the transposition,
C            (   :   )
C            ( u(t)' )
C
C   and the output trajectories result in similarly.
C
C Input parameters:
C   A      - the n-by-n state matrix A.
C   B      - the n-by-m input matrix B.
C   C      - the p-by-n output matrix C.
C   D      - the p-by-m input-output matrix D.
C   U      - the t-by-m matrix U.
C   x      - (optional) the initial state x(1).
C            Default: x = 0.
C   ldwork - (optional) the length of working array.
C            Default: ldwork = n.
C            Larger values could increase the efficiency.
C
C Output parameters:
C   Y      - the t-by-p output matrix Y.
C   x      - (optional) the final state x(t+1).
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C   June 2012.
C
C Revisions:
C   Mar. 2002, Apr. 2009.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      mwPointer         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCalloc, mxCreateDoubleMatrix,
     $                  mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDU, LDWORK, LDY, M,
     $                  N, NY, P
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      mwPointer         A, B, C, D, DWORK, U, X, Y
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ITMP, LDWMIN
      DOUBLE PRECISION  DUM(1), TEMP
C
C .. External subroutines ..
      EXTERNAL          DCOPY, TF01MY
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LDSIM requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LDSIM requires at least 1 output argument' )
      END IF
C
C   A(nxn), B(nxm), C(pxn), D(pxm), and U(txm).
C
      N  = mxGetM( PRHS(1) )
      M  = mxGetN( PRHS(2) )
      P  = mxGetM( PRHS(3) )
      NY = mxGetM( PRHS(5) )
C
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(2) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same number of rows as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $           ( 'C must have the same number of columns as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(4) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same number of rows as C' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $           ( 'D must have the same number of columns as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $           ( 'U must have the same number of columns as B' )
      END IF
C
      IF ( NRHS.GE.6 ) THEN
C
C   x(1)
C
         IF ( mxGetM( PRHS(6) )*mxGetN( PRHS(6) ).LT.N ) THEN
            WRITE( TEXT, '(''X must have '',I7,'' entries'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
      END IF
C
C   ldwork
C
      LDWMIN = MAX( 1, N )
      IF ( NRHS.GE.7 ) THEN
         IF ( mxGetM( PRHS(7) ).NE.1 .OR.
     $        mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'LDWORK must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(7) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'LDWORK must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TEMP, 1 )
         ITMP = INT( TEMP )
         IF ( ITMP.LT.LDWMIN ) THEN
            LDWORK = LDWMIN
         ELSE
            LDWORK = ITMP
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
      LDU = MAX( 1, NY )
      LDY = LDU
      IF ( NRHS.LT.7 ) THEN
         CALL TF01MY( N, M, P, NY, DUM, LDA, DUM, LDB, DUM, LDC, DUM,
     $             LDD, DUM, LDU, DUM, DUM, LDY, DUM, -1, INFO )
         LDWORK = DUM(1)
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      A     = mxCalloc(  LDA*N, 8 )
      B     = mxCalloc(  LDB*M, 8 )
      C     = mxCalloc(  LDC*N, 8 )
      D     = mxCalloc(  LDD*M, 8 )
      DWORK = mxCalloc( LDWORK, 8 )
      U     = mxCalloc(  LDU*M, 8 )
      X     = mxCalloc(      N, 8 )
      Y     = mxCalloc(  LDY*P, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL(A), N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL(B), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL(C), P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL(D), P*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL(U), NY*M )
      IF ( NRHS.GE.6 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), %VAL(X), N )
      ELSEIF ( N.GT.0 ) THEN
         DUM(1) = ZERO
         CALL DCOPY( N, DUM, 0, %VAL(X), 1 )
      END IF
C
C Do the actual computations.
C
      CALL TF01MY( N, M, P, NY, %VAL(A), LDA, %VAL(B), LDB, %VAL(C),
     $             LDC, %VAL(D), LDD, %VAL(U), LDU, %VAL(X), %VAL(Y),
     $             LDY, %VAL(DWORK), LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( NY, P, 0 )
      CALL mxCopyReal8ToPtr( %VAL(Y), mxGetPr( PLHS(1) ), NY*P )
      IF ( NLHS.GE.2 ) THEN
         PLHS(2) = mxCreateDoubleMatrix( N, MIN( N, 1 ), 0 )
         CALL mxCopyReal8ToPtr( %VAL(X), mxGetPr( PLHS(2) ), N )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( D     )
      CALL mxFree( DWORK )
      CALL mxFree( U     )
      CALL mxFree( X     )
      CALL mxFree( Y     )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TF01MY'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of LDSIM ***
      END
