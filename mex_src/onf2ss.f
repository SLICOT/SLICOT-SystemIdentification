C ONF2SS.F - Gateway function for the conversion of a linear
C            discrete-time system given in the output normal form to
C            the state-space representation using SLICOT routine TB01VY.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [A,B,C,D,x0] = onf2ss(n,m,l,theta,apply)  % if l > 0;
C     [B,C,D,x0] = onf2ss(n,m,l,theta,apply)  % if l = 0.
C
C Purpose:
C  To convert the linear discrete-time system given as its output
C  normal form, with parameter vector THETA, into the state-space
C  representation (A, B, C, D), with the initial state x0.
C
C Input parameters:
C   n      - the order of the system.
C   m      - the number of the system inputs.
C   l      - the number of the system outputs.
C   theta  - the n*(l+m+1)+l*m parameter vector.
C            If l = 0, matrix A cannot be recovered.
C   apply  - (optional) integer specifying whether or not the parameter
C            vector should be transformed using a bijective mapping:
C            = 1: apply the bijective mapping to the n vectors in
C                 theta corresponding to the matrices A and C;
C            = 0: do not apply the bijective mapping.
C            Default:  apply = 0.
C            The value of apply should be the same as that used in the
C            corresponding call of the paired SS2ONF function.
C
C Output parameters:
C   A      - the n-by-n state matrix A.
C   B      - the n-by-m input matrix B.
C   C      - the l-by-n output matrix C.
C   D      - the l-by-m input-output matrix D.
C   x0     - the initial state.
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
      INTEGER           IP, ITMP
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          TB01VY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ONF2SS requires at least 4 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ONF2SS requires at least 1 output arguments' )
      END IF
C
C   n, m, l, theta(n*(l+m+1)+l*m), apply.
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'N must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'N must be a real scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      N = TEMP
      IF ( N.LT.0 ) THEN
         CALL mexErrMsgTxt( 'N must be a non-negative integer' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR.
     $     mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'M must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'M must be a real scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      M = TEMP
      IF ( M.LT.0 ) THEN
         CALL mexErrMsgTxt( 'M must be a non-negative integer' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $     mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'L must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'L must be a real scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      L = TEMP
      IF ( L.LT.0 ) THEN
         CALL mexErrMsgTxt( 'L must be a non-negative integer' )
      END IF
C
      LTHETA = N*( M + L + 1 ) + L*M
C
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'THETA must be a real vector' )
      END IF
      IF ( mxGetM( PRHS(4) )*mxGetN( PRHS(4) ).NE.LTHETA ) THEN
         WRITE( TEXT, '('' THETA must have = '',I7,'' elements'')' )
     $          LTHETA
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      ITMP = 0
      IF ( NRHS.GE.5 ) THEN
         IF ( mxGetM( PRHS(5) ).NE.1 .OR.
     $        mxGetN( PRHS(5) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'APPLY must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'APPLY must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, 1 )
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
      LDWORK = N*( N + L + 1 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           DWORK( LDWORK ), THETA( LTHETA ), X0( N ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), THETA, LTHETA )
C
C Do the actual computations.
C
      CALL TB01VY( APPLY, N, M, L, THETA, LTHETA, A, LDA, B, LDB, C,
     $             LDC, D, LDD, X0, DWORK, LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      IF ( L.GT.0 ) THEN
         PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
         CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), N*N )
         IP = 2
      ELSE
         IP = 1
      END IF
      PLHS(IP) = mxCreateDoubleMatrix( N, M, 0 )
      CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(IP) ), N*M )
      PLHS(IP+1) = mxCreateDoubleMatrix( L, N, 0 )
      CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(IP+1) ), L*N )
      PLHS(IP+2) = mxCreateDoubleMatrix( L, M, 0 )
      CALL mxCopyReal8ToPtr( D, mxGetPr( PLHS(IP+2) ), L*M )
      PLHS(IP+3) = mxCreateDoubleMatrix( N, 1, 0 )
      CALL mxCopyReal8ToPtr( X0, mxGetPr( PLHS(IP+3) ), N )
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, DWORK, THETA, X0 )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM TB01VY'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of ONF2SS ***
      END
