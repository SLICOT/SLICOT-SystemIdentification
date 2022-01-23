C SIDENT.F - Gateway function for computation of a discrete-time
C            state-space realization (A,B,C,D) and Kalman gain
C            using SLICOT routine IB01BD.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [(A,C)(,B(,D))(,K,Q,Ry,S,rcnd)] = sident(meth,job,s,n,l,R(,tol,t,A,
C                                            C,printw))
C
C Purpose:
C   To compute a state-space realization (A,B,C,D) and the Kalman
C   predictor gain K of a discrete-time system, given the system
C   order and the relevant part of the R factor of the concatenated
C   block-Hankel matrices, using subspace identification techniques
C   (MOESP and N4SID).
C
C Input parameters:
C   meth  - integer option to determine the method to use:
C           = 1 : MOESP method with past inputs and outputs;
C           = 2 : N4SID method;
C           = 3 : combined method: A and C via MOESP, B and D via N4SID.
C   job   - integer option to determine the calculation to be performed:
C           = 1 : compute all system matrices, A, B, C, D;
C           = 2 : compute the matrices A and C only;
C           = 3 : compute the matrix B only;
C           = 4 : compute the matrices B and D only.
C   s     - the number of block rows in the processed input and output
C           block Hankel matrices.  s > 0.
C   n     - the order of the system.
C   l     - the number of system outputs.
C   R     - the 2*(m+l)*s-by-2*(m+l)*s part of  R  contains the
C           processed upper triangular factor  R  from the
C           QR factorization of the concatenated block-Hankel matrices,
C           and further details needed for computing system matrices.
C           (Above, m denotes the number of system inputs, determined
C           by s, l, and the size of R.)
C   tol   - (optional) tolerance used for estimating the rank of
C           matrices. If  tol > 0,  then the given value of  tol  is
C           used as a lower bound for the reciprocal condition number;
C           an m-by-n matrix whose estimated condition number is less
C           than  1/tol  is considered to be of full rank.
C           Default:    m*n*epsilon_machine where epsilon_machine is
C           the relative machine precision.
C   t     - (optional) the total number of samples used for calculating
C           the covariance matrices.  Either  t = 0, or  t >= 2*(m+l)*s.
C           This parameter is not needed if the covariance matrices
C           and/or the Kalman predictor gain matrix are not desired.
C           If t = 0, then K, Q, Ry, and S are not computed.
C           Default:    t = 0.
C   A     - (optional) the n-by-n system state matrix A.
C           This parameter is needed if meth >= 2 and job >= 3.
C   C     - (optional) the l-by-n system output matrix C.
C           This parameter is needed if meth >= 2 and job >= 3.
C   printw- (optional) switch for printing the warning messages.
C           = 1:  print warning messages;
C           = 0:  do not print warning messages.
C           Default:    printw = 0.
C
C Output parameters:
C   A     - if job <= 2, the n-by-n system state matrix A.
C   C     - if job <= 2, the l-by-n system output matrix C.
C   B     - if job <> 2, the n-by-m system input matrix B.
C   D     - if job = 1 or 4, the l-by-m system matrix D.
C   K     - (optional) the n-by-l Kalman predictor gain matrix K.
C   Q     - (optional) the n-by-n positive semidefinite state covariance
C           matrix used as state weighting matrix when computing the
C           Kalman gain.
C   Ry    - (optional) the l-by-l positive (semi)definite output
C           covariance matrix used as output weighting matrix when
C           computing the Kalman gain.
C   S     - (optional) the n-by-l state-output cross-covariance matrix
C           used as cross-weighting matrix when computing the Kalman
C           gain.
C   rcnd  - (optional) vector of length lr, containing estimates of the
C           reciprocal condition numbers of the matrices involved in
C           rank decisions, least squares or Riccati equation solutions,
C           where lr = 4,  if Kalman gain matrix K is not required, and
C                 lr = 12, if Kalman gain matrix K is required.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 1999.
C
C Revisions:
C   V. Sima, May 2000, July 2000, March 2005, Apr. 2009, Dec. 2012.
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
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      DOUBLE PRECISION  TOL
      INTEGER           INFO, IWARN, L, LDA, LDB, LDC, LDD, LDK, LDO,
     $                  LDQ, LDR, LDRY, LDS, LDWORK, M, N, NOBR, NSMPL
      CHARACTER         JOB, JOBCK, METH
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                 DWORK(:), K(:,:), Q(:,:), R(:,:),
     $                                 RY(:,:), S(:,:)
C
C .. Local variables and constant dimension arrays ..
      DOUBLE PRECISION  TEMP
      CHARACTER*120     TEXT
      INTEGER           ID, IJOB, IP, ITMP, LBWORK, LDUNN, LIWORK, LL,
     $                  LNOBR, MA, MNOBR, MNOBRN, N2, NA, NCOL, NL, NN,
     $                  NPL, NR, NRC, TASK
      LOGICAL           PRINTW
C
C .. External subroutines ..
      EXTERNAL          IB01BD
C
C ..Intrinsic functions..
      INTRINSIC         MAX
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SIDENT requires at least 6 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SIDENT requires at least 1 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   meth
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      TASK = TEMP
      IF ( TASK.LT.1 .OR. TASK.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1, 2, or 3 the only admissible values' )
      END IF
C
      IF ( TASK.EQ.1 ) THEN
         METH = 'M'
      ELSE IF ( TASK.EQ.2 ) THEN
         METH = 'N'
      ELSE
         METH = 'C'
      END IF
C
C   job
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      IJOB = TEMP
      IF ( IJOB.LT.1 .OR. IJOB.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'JOB has 1, 2, 3 or 4 the only admissible values' )
      END IF
C
      IF ( IJOB.EQ.1 ) THEN
         JOB = 'A'
      ELSE IF ( IJOB.EQ.2 ) THEN
         JOB = 'C'
      ELSE IF ( IJOB.EQ.3 ) THEN
         JOB = 'B'
      ELSE
         JOB = 'D'
      END IF
C
C   s
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 's must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 's must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      NOBR = TEMP
      IF ( NOBR.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 's must be a positive integer' )
      END IF
C
C   n
C
      IF ( mxGetM( PRHS(4) ).NE.1 .OR. mxGetN( PRHS(4) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'n must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'n must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), TEMP, 1 )
      N = TEMP
      IF ( N.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'n must be a positive integer' )
      END IF
      IF ( N.GE.NOBR ) THEN
         WRITE( TEXT, '('' The order should be at most '', I6)' )
     $          NOBR - 1
         CALL mexErrMsgTxt( TEXT )
      END IF
C
C   l
C
      IF ( mxGetM( PRHS(5) ).NE.1 .OR. mxGetN( PRHS(5) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'l must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'l must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, 1 )
      L = TEMP
      IF ( L.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'The system has no outputs' )
      END IF
C
C   R(nr,nr)
C
      NR   = mxGetM( PRHS(6) )
      NCOL = mxGetN( PRHS(6) )
      IF ( NR.LT.2*L ) THEN
         WRITE( TEXT, '(''R must have at least '',I5,'' rows'')' ) 2*L
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( NCOL.LT.NR ) THEN
         WRITE( TEXT, '(''R must have at least '',I5,'' columns'')' )
     $          NCOL
         CALL mexErrMsgTxt( TEXT )
      END IF
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'R must be a real matrix' )
      END IF
C
C   m
C
      M = NR / ( 2*NOBR ) - L
C
C   tol
C
      TOL = ZERO
      IF ( NRHS.GT.6 ) THEN
         IF ( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(7) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TOL, 1 )
      END IF
C
C   t
C
      NSMPL = 0
      JOBCK = 'N'
      IF ( NRHS.GT.7 ) THEN
         JOBCK = 'K'
         IF ( mxGetM( PRHS(8) ).NE.1 .OR. mxGetN( PRHS(8) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 't must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(8) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 't must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), TEMP, 1 )
         NSMPL = TEMP
C
         IF ( NSMPL.NE.0 .AND. NSMPL.LT.NR ) THEN
            WRITE( TEXT,
     $             '('' The number of samples should be at least '',
     $                  I10)' ) NR
            CALL mexErrMsgTxt( TEXT )
         ELSE IF ( NSMPL.EQ.0 ) THEN
            JOBCK = 'N'
         END IF
      END IF
C
C   A(n,n)
C
      IF ( TASK.GE.2 .AND. IJOB.GE.3 ) THEN
         MA = mxGetM( PRHS(9) )
         NA = mxGetN( PRHS(9) )
         IF ( MA.NE.N .OR. NA.NE.N ) THEN
            WRITE( TEXT, '(''A must have '',I5,'' rows and columns'')' )
     $             N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'A must be a real matrix' )
         END IF
C
C   C(l,n)
C
         MA = mxGetM( PRHS(10) )
         NA = mxGetN( PRHS(10) )
         IF ( MA.NE.L ) THEN
            WRITE( TEXT, '(''C must have '',I5,'' rows'')' ) L
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( NA.NE.N ) THEN
            WRITE( TEXT, '(''C must have '',I5,'' columns'')' ) N
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'C must be a real matrix' )
         END IF
      END IF
C
C     printw
C
      PRINTW = .FALSE.
      IP = 11
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'PRINTW must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'PRINTW must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.LT.0 .OR. ITMP.GT.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'PRINTW has 0 or 1 the only admissible values' )
         END IF
         PRINTW = ITMP.EQ.1
      END IF
C
C Determine the lenghts of working arrays.
C The value for LDWORK is the minimum value needed by IB01BD for each
C method and algorithm implemented.  Using a larger value could
C increase the efficiency.
C
      MNOBR  = M*NOBR
      LNOBR  = L*NOBR
      MNOBRN = MNOBR + N
      LDUNN  = ( LNOBR - L )*N
      NPL    = N + L
      N2     = N + N
      NN     = N*N
      NL     = N*L
      LL     = L*L
C
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, L )
      LDD = LDC
      LDO = LNOBR
      LDR = NR
      IF ( NSMPL.NE.0 ) THEN
         LDK    = LDA
         LDQ    = LDA
         LDS    = LDA
         LDRY   = LDC
         LBWORK = N2
      ELSE
         LDK    = 1
         LDQ    = 1
         LDS    = 1
         LDRY   = 1
         LBWORK = 1
      END IF
C
      LIWORK = MNOBR + N
      IF ( TASK.EQ.1 ) THEN
         LIWORK = MAX( LIWORK, LNOBR )
      ELSE IF ( TASK.EQ.2 ) THEN
         LIWORK = MAX( LIWORK, M*NPL )
      ELSE
         LIWORK = MAX( LIWORK, LNOBR, M*NPL )
      END IF
      IF ( NSMPL.GT.0 )
     $   LIWORK = MAX( LIWORK, NN )
C
      IAW    = 0
      LDWORK = LDUNN + 4*N
      IF ( TASK.EQ.1 ) THEN
         ID = 0
      ELSE
         ID = N
      END IF
C
      IF ( TASK.NE.2 ) THEN
         IF ( IJOB.LE.2 ) THEN
            LDWORK = MAX( LDWORK, 2*LDUNN + N2, LDUNN + NN + 7*N )
         END IF
      END IF
C
      IF ( ( M.GT.0 .AND. IJOB.NE.2 ) .OR. TASK.GE.2 ) THEN
         LDWORK = MAX( LDWORK, 2*LDUNN + NN + ID + 7*N )
         IF ( TASK.EQ.1 )
     $      LDWORK = MAX( LDWORK, LDUNN + N + 6*MNOBR, LDUNN + N +
     $                    MAX( L + MNOBR, LNOBR +
     $                                    MAX( 3*LNOBR + 1, M ) ) )
      ELSE
         IF ( TASK.NE.2 ) THEN
     $      IAW = N + NN
      END IF
C
      IF ( TASK.NE.1 .OR. NSMPL.GT.0 ) THEN
         LDWORK = MAX( LDWORK, LDUNN + IAW + N2 +
     $                 MAX( 5*N, LNOBR + 2*MNOBR + L ),
     $                 ID + 4*MNOBRN + 1, ID + MNOBRN + NPL )
         IF ( TASK.NE.1 .AND. M.GT.0 .AND. IJOB.NE.2 )
     $      LDWORK = MAX( LDWORK, MNOBR*NPL*( M*NPL + 1 ) +
     $                    MAX( NPL**2, 4*M*NPL + 1 ) )
         LDWORK = LNOBR*N + LDWORK
      END IF
C
      IF ( NSMPL.GT.0 )
     $   LDWORK = MAX( LDWORK, 4*NN + 2*NL + LL + MAX( 3*L, NL ),
     $                 14*NN + 12*N + 5 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( A( LDA, N ), B( LDB, M ), C( LDC, N ), D( LDD, M ),
     $           DWORK( LDWORK ), IWORK( LIWORK ), Q( LDQ, N ),
     $           R( LDR, NCOL ), RY( LDRY, L ), S( LDS, L ) )
      ALLOCATE ( BWORK( LBWORK ), K( LDK, L ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), R, LDR*NCOL )
      IF ( TASK.GE.2 .AND. IJOB.GE.3 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ),  A, LDA*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ), C, LDC*N )
      END IF
C
C Do the actual computations.
C
      CALL IB01BD( METH, JOB, JOBCK, NOBR, N, M, L, NSMPL, R, LDR, A,
     $             LDA, C, LDC, B, LDB, D, LDD, Q, LDQ, RY, LDRY, S,
     $             LDS, K, LDK, TOL, IWORK, DWORK, LDWORK, BWORK, IWARN,
     $             INFO )
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM IB01BD'')'
     $        ) IWARN
      END IF
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '(''INFO = '',I4,'' ON EXIT FROM IB01BD'')'
     $        ) INFO
      ELSE
C
C Copy output to MATLAB workspace.
C
         IF ( IJOB.LE.2 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), NN )
            IP = 1
            IF ( NLHS.GT.1 ) THEN
               PLHS(2) = mxCreateDoubleMatrix( L, N, 0 )
               CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(2) ), NL )
               IP = 2
            END IF
         ELSE
            IP = 0
         END IF
C
         IF ( NLHS.GT.IP ) THEN
            IF ( IJOB.EQ.1 .OR. IJOB.GE.3 ) THEN
               IP = IP + 1
               PLHS(IP) = mxCreateDoubleMatrix( N, M, 0 )
               CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(IP) ), N*M )
            END IF
            IF ( NLHS.GT.IP ) THEN
               IF ( IJOB.EQ.1 .OR. IJOB.EQ.4 ) THEN
                  IP = IP + 1
                  PLHS(IP) = mxCreateDoubleMatrix( L, M, 0 )
                  CALL mxCopyReal8ToPtr( D, mxGetPr( PLHS(IP) ), L*M )
               END IF
            END IF
         END IF
C
         IF ( NSMPL.GT.0 .AND. NLHS.GT.IP ) THEN
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( N, L, 0 )
            CALL mxCopyReal8ToPtr( K, mxGetPr( PLHS(IP) ), NL )
         END IF
C
         IF ( NSMPL.GT.0 .AND. NLHS.GT.IP ) THEN
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( N, N, 0 )
            CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(IP) ), NN )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( L, L, 0 )
            CALL mxCopyReal8ToPtr( RY, mxGetPr( PLHS(IP) ), LL )
            IP = IP + 1
            PLHS(IP) = mxCreateDoubleMatrix( N, L, 0 )
            CALL mxCopyReal8ToPtr( S, mxGetPr( PLHS(IP) ), NL )
         END IF
C
         IF ( NLHS.GT.IP ) THEN
            IP = IP + 1
            IF ( NSMPL.EQ.0 ) THEN
               NRC = 4
            ELSE
               NRC = 12
            END IF
            PLHS(IP) = mxCreateDoubleMatrix( NRC, 1, 0 )
            CALL mxCopyReal8ToPtr( DWORK(2), mxGetPr( PLHS(IP) ), NRC )
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( A, B, C, D, R, Q, RY, S, IWORK, DWORK )
      DEALLOCATE( K, BWORK )
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SIDENT ***
      END
