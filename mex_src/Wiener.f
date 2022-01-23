C WIENER.F - Gateway function for computing the output of a Wiener
C            system using SLICOT routine NF01AD.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   y = Wiener(n,l,nn,theta,u(,ldwork))
C
C Purpose:
C  To calculate the output y of the Wiener system
C
C        x(k+1) = A*x(k) + B*u(k)
C        z(k)   = C*x(k) + D*u(k),
C
C        y(k)   = f(z(k),wb(1:l)),
C
C  where the linear discrete-time system is given as its output normal
C  form, with parameter vector THETA(l1+1:l1+l2), and the parameters of
C  the nonlinear part are contained in THETA(1:l1), with
C  l1 = (nn*(l+2)+1)*l, l2 = n*(l+m+1)+l*m.
C
C Input parameters:
C   n      - the order of the linear system.
C   l      - the number of the system outputs.
C   nn     - the number of neurons of the nonlinear part.
C   theta  - the (nn*(l+2)+1)*l+n*(l+m+1)+l*m parameter vector.
C   u      - the t-by-m input trajectory.
C   ldwork - (optional) the length of working array.
C            Default: ldwork = t*l+ max( 2*nn,(n+l)*(n+m)+2*n+
C                                             max(n*(n+l),nml)),
C                              where nml = n+m+l, if m > 0, and
C                                    nml = l,     if m = 0.
C            Larger values could increase the efficiency.
C
C Output parameters:
C   y      - the t-by-l output trajectory.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
C
C Revisions:
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
      INTEGER           INFO, L, LDU, LDWORK, LDY, LIPAR, LTHETA, M,
     $                  NSMP
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER, ALLOCATABLE          :: IPAR(:)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:), THETA(:), U(:,:),
     $                                 Y(:,:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      INTEGER           ITMP, LDWMIN, N, NN
      DOUBLE PRECISION  TEMP
C
C .. External subroutines ..
      EXTERNAL          NF01AD
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'WIENER requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'WIENER requires 1 output arguments' )
      END IF
C
C   n, l, nn, theta((nn*(l+2)+1)*l+n*(l+m+1)+l*m), u(t*m).
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'N must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'N must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      N = TEMP
      IF ( N.LT.0 ) THEN
         CALL mexErrMsgTxt( 'N must be a non-negative integer' )
      END IF
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR.
     $     mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'L must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'L must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      L = TEMP
      IF ( L.LT.0 ) THEN
         CALL mexErrMsgTxt( 'L must be a non-negative integer' )
      END IF
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR.
     $     mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NN must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NN must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      NN = TEMP
      IF ( NN.LT.0 ) THEN
         CALL mexErrMsgTxt( 'NN must be a non-negative integer' )
      END IF
C
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
      NSMP = mxGetM( PRHS(5) )
      M    = mxGetN( PRHS(5) )
      IF ( NSMP.LT.0 .OR. M.LT.0 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
C
      LTHETA = ( NN*( L + 2 ) + 1 )*L + N*( M + L + 1 ) + L*M
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
C   ldwork
C
      IF ( M.GT.0 ) THEN
         LDWMIN = MAX( N*( N + L ), N + M + L )
      ELSE
         LDWMIN = MAX( N*( N + L ), L )
      END IF
      LDWMIN = MAX( 1, NSMP*L + MAX( 2*NN, ( N + L )*( N + M ) + 2*N +
     $                              LDWMIN ) )
      IF ( NRHS.GE.6 ) THEN
         IF ( mxGetM( PRHS(6) ).NE.1 .OR.
     $        mxGetN( PRHS(6) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'LDWORK must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(6) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'LDWORK must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TEMP, 1 )
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
      LDU = MAX( 1, NSMP )
      LDY = LDU
      IF ( NRHS.LT.6 )
     $   LDWORK = LDWMIN
      LIPAR = 2
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( DWORK( LDWORK ), IPAR( LIPAR ), THETA( LTHETA ),
     $           U( LDU, M ), Y( LDY, L ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), THETA, LTHETA )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), U, NSMP*M )
C
C Do the actual computations.
C
      IPAR(1) = N
      IPAR(2) = NN
      CALL NF01AD( NSMP, M, L, IPAR, LIPAR, THETA, LTHETA, U, LDU,
     $             Y, LDY, DWORK, LDWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( NSMP, L, 0 )
      CALL mxCopyReal8ToPtr( Y, mxGetPr( PLHS(1) ), NSMP*L )
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( DWORK, IPAR, THETA, U, Y )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM NF01AD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of WIENER ***
      END
