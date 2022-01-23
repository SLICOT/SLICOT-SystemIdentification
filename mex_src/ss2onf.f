C SS2ONF.F - Gateway function for the conversion of a linear
C            discrete-time system into the output normal form
C            using SLICOT routine TB01VD.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   theta = ss2onf(A,B,C,D,x0,apply)
C
C Purpose:
C  To convert the linear discrete-time system given as (A, B, C, D),
C  with initial state x0, into the output normal form, with parameter
C  vector theta. The matrix A is assumed to be stable.
C  Optionally, a bijective transformation is used, guaranteeing that
C  the inverse conversion (using ONF2SS) can always be performed.
C
C Input parameters:
C   A      - the n-by-n state matrix A.
C   B      - the n-by-m input matrix B.
C   C      - the l-by-n output matrix C.
C   D      - the l-by-m input-output matrix D.
C   x0     - the initial state.
C   apply  - (optional) integer specifying whether or not the parameter
C            vector should be transformed using a bijective mapping:
C            = 1: apply the bijective mapping to the n vectors in
C                 theta corresponding to the matrices A and C;
C            = 0: do not apply the bijective mapping.
C            Default:  apply = 0.
C
C Output parameters:
C   theta  - the n*(l+m+1)+l*m parameter vector.
C
C Contributors:
C   A. Riedel, R. Schneider, Chemnitz Univ. of Technology, Mar. 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C   Feb. 2002, Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         APPLY
      INTEGER           INFO, L, LDA, LDB, LDC, LDD, LDWORK, LTHETA, M,
     $                  N
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                 DWORK(:), THETA(:), X0(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ITMP
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          TB01VD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SS2ONF requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SS2ONF requires 1 output argument' )
      END IF
C
C   A(nxn), B(nxm), C(lxn), D(lxm), x0(n), apply.
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(2) )
      L = mxGetM( PRHS(3) )
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
      IF ( mxGetM( PRHS(4) ).NE.L ) THEN
         CALL mexErrMsgTxt( 'D must have the same number of rows as C' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $           ( 'D must have the same number of columns as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'X0 must be a real matrix' )
      END IF
      IF ( mxGetM( PRHS(5) )*mxGetN( PRHS(5) ).LT.N ) THEN
         WRITE( TEXT, '(''X0 must have '',I7,'' entries'')' ) N
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      ITMP = 0
      IF ( NRHS.GE.6 ) THEN
         IF ( mxGetM( PRHS(6) ).NE.1 .OR.
     $        mxGetN( PRHS(6) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'APPLY must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'APPLY must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
            CALL mexErrMsgTxt
     $              ( 'APPLY has 0 and 1 the only admissible values' )
         END IF
      END IF
      IF ( ITMP.EQ.0 ) THEN
         APPLY = 'N'
      ELSE
         APPLY = 'A'
      END IF
C
C Determine the lenghts of working arrays.
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, L )
      LDD = LDC
C
      LTHETA = N*( L + M + 1 ) + L*M
      LDWORK = MAX( 1, N*N*L + N*L + N,
     $              N*N + MAX( N*( N + MAX( N, L ) + 6 ) + MIN( N, L ),
     $                         N*M ) )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           DWORK( LDWORK ), THETA( LTHETA ), X0( N ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), C, L*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), D, L*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), X0, N )
C
C Do the actual computations.
C
      CALL TB01VD( APPLY, N, M, L, A, LDA, B, LDB, C, LDC, D, LDD, X0,
     $             THETA, LTHETA, DWORK, LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( LTHETA, 1, 0 )
      CALL mxCopyReal8ToPtr( THETA, mxGetPr( PLHS(1) ), LTHETA )
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, DWORK, THETA, X0 )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TB01VD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SS2ONF ***
      END
