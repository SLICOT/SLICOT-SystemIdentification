C WIDENTC.F - Gateway function for computation of a discrete-time
C             model of a Wiener system using SLICOT routine IB03AD.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [xopt(,perf,nf,rcnd)] = widentc(job,u,y,nn(,s,n,x,alg,stor,iter,
C                                   nprint,tol,seed,printw,ldwork))
C
C   [xopt(,perf,nf,rcnd)] = widentc(1,u,y,nn,s,n,x(,alg,stor,iter,
C                                   nprint,tol,printw,ldwork))
C   [xopt(,perf,nf)]      = widentc(2,u,y,nn,n,x(,alg,stor,iter,nprint,
C                                   tol,seed,printw,ldwork))
C   [xopt(,perf,nf,rcnd)] = widentc(3,u,y,nn,s(,n,[],alg,stor,iter,
C                                   nprint,tol,seed,printw,ldwork))
C   [xopt(,perf,nf)]      = widentc(4,u,y,nn,n,x(,alg,stor,iter,nprint,
C                                   tol,printw,ldwork))
C
C Purpose:
C   To compute a set of parameters for approximating a Wiener system
C   in a least-squares sense, using a neural network approach and a
C   Levenberg-Marquardt algorithm, and a Cholesky or conjugate gradients
C   algorithm for solving linear systems of equations.
C
C Input parameters:
C   job    - integer specifying which parameters must be initialized:
C            = 1 : initialize the linear part only;
C            = 2 : initialize the static nonlinearity only;
C            = 3 : initialize both linear and nonlinear parts;
C            = 4 : do not initialize anything, x already contains
C                  an initial approximation.
C   u      - the t-by-m input trajectory, where t is the number of
C            samples.
C   y      - the t-by-l output trajectory.
C   nn     - the number of neurons used to approximate the nonlinear
C            part.
C   s      - if job = 1 or 3, the number of block rows in the input and
C            output block Hankel matrices to be processed.  s > 0.
C            The argument s must not be specified if job = 2 or 4.
C   n      - (optional) the order of the linear part of the system or
C            an option for how to compute it.
C            If job = 1 or 3, and n < 0, the order will be found by the
C            program. Otherwise, n should be non-negative (or positive,
C            but less than s, if job = 1 or 3).
C            The value of n should be given if job = 2 or 4.
C            Default for job = 1 or 3: the order found automatically.
C   x      - the lx0 vector of initial parameters, where
C            lx0 = l1,  if job = 1,
C            lx0 = l2,  if job = 2,
C            lx0 = lx,  if job = 4,
C            with l1 = (nn*(l+2)+1)*l, l2 = n*(l+m+1)+l*m, lx = l2+l1.
C            Here, l1 is the number of parameters for the nonlinear part
C            and l2 is the number of parameters for the linear part.
C            If n < 0 on entry, then s-1 is used instead of n for
C            setting l2 and lx, and for allocating the workspace.
C            If job = 3, the vector x is not needed.
C   alg    - (optional) integer specifying the algorithm used for
C            solving the linear systems involving a Jacobian matrix J:
C            = 1 :  a direct algorithm, which computes the Cholesky
C                   factor of the matrix J'*J + par*I is used, where
C                   par is the Levenberg factor;
C            = 2 :  an iterative Conjugate Gradients algorithm, which
C                   only needs the matrix J, is used.
C            Default: alg = 1.
C   stor   - (optional) integer specifying the storage scheme for the
C            symmetric matrix J'*J, if alg = 1:
C            = 1 :  full storage is used;
C            = 2 :  packed storage is used.
C            Default: stor = 1.
C            The option stor = 1 usually ensures a faster execution.
C            This parameter is not relevant if alg = 2.
C   iter   - (optional) vector of length 2 containing the maximal
C            numbers of iterations for the Levenberg-Marquardt algorithm
C            for the initialization of the static nonlinearity (ignored
C            if job = 1 or 4), and for the optimization process; these
C            numbers are stored in iter(1) and iter(2), respectively.
C            Default: iter = [10*lx 10*lx].
C   nprint - (optional) integer specifying the frequency of printing
C            details of the optimization process. If nprint > 0, the
C            intermediate results are printed in the file IB03AD.prn,
C            which is overwritten at each execution of this mexfile.
C            Default: nprint = 0 (no printing).
C   tol    - (optional) vector of length 2 containing the tolerances
C            for the initialization of the static nonlinearity
C            (ignored if job = 1 or 4), and for the optimization
C            process, respectively.
C            Default: tol = [sqrt(eps) sqrt(eps)], with
C            eps = epsilon_machine.
C   seed   - (optional) if job = 2 or 3, vector of length 4 containing
C            the random number generator seed used to initialize the
C            parameters of the static nonlinearity.
C   printw - (optional) switch for printing the warning messages.
C            = 1:  print warning messages;
C            = 0:  do not print warning messages.
C            Default:    printw = 0.
C   ldwork - (optional) the length of working array.
C            Default: the minimum workspace needed, computed by the
C                     program based on the specified input arguments.
C            Larger values could increase the efficiency.
C
C Output parameters:
C   xopt   - the lx vector of optimal parameters.
C            The corresponding linear system and its initial state can
C            be found as follows:
C            m = size(u,2);  l = size(y,2);  l1 = (nn*(l+2)+1)*l;
C            n = ( length(xopt) - l1 - l*m ) / ( l+m+1 );
C            [A,B,C,D,x0] = onf2ss(n,m,l,xopt(l1+1:end),1);
C   perf   - (optional) vector of length 5 or 10 containing performance
C            results:
C            perf(1)  contains the optimal value of the workspace
C                     length;
C            perf(2)  contains the maximum residual error norm;
C            perf(3)  contains the total number of iterations performed;
C            perf(4)  contains the total number of CG iterations
C                     performed, if alg = 2;
C            perf(5)  contains the final Levenberg factor.
C            If job = 2 or 3, then similar results for the nonlinear
C            part initialization step are returned:
C            perf(6)  contains the optimal value of the workspace
C                     length;
C            perf(7)  contains the maximum residual error norm;
C            perf(8)  contains the total number of iterations performed;
C            perf(9)  contains the total number of CG iterations
C                     performed, if alg = 2;
C            perf(10) contains the final Levenberg factor (maximum over
C                     all outputs).
C   nf     - (optional) 2-vector containing additional performance
C            results:
C            nf(1) contains the (total) number of function evaluations;
C            nf(2) contains the (total) number of Jacobian evaluations.
C   rcnd   - (optional) if job = 1 or 3, vector of suitable length
C            containing the reciprocal condition number estimates for
C            determining the linear part by susbspace techniques.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
C
C Revisions:
C   V. Sima, Mar. 2002, Apr. 2009, Dec. 2012.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
C     .. NOUT is the unit number for printing intermediate results ..
C     .. if NPRINT is positive ..
      INTEGER           NOUT
      PARAMETER         ( NOUT = 6 )
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
      CHARACTER         ALG, INIT, STOR
      INTEGER           INFO, ITMAX1, ITMAX2, IWARN, L, LDU, LDWORK,
     $                  LDY, LX, M, N, NN, NOBR, NPRINT, NSMP
      DOUBLE PRECISION  TOL1, TOL2
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:), U(:,:), X(:), Y(:,:)
      INTEGER,          ALLOCATABLE :: IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           INIT1, INIT2, PRINTW
      INTEGER           IALG, ISTOR, BSN, INI, IP, ISAD, ISIZE, ITMP,
     $                  IW1, IW2, IX, JWORK, L1, L2, LDAC, LDWMIN,
     $                  LIWORK, LNOL, LPAR, LPAR0, M1, M2, ML, MNO, N2,
     $                  NS, NSML, NX, TASK
      DOUBLE PRECISION  RSEED(4), TEMP, TMP(4)
C
C .. External functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C
C .. External subroutines ..
      EXTERNAL          DCOPY, IB03AD
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'WIDENTC requires at least 4 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'WIDENTC requires at least 1 output arguments' )
      END IF
C
C   job, u(t*m), y(t*l), nn, s, n, x(lx0), alg, stor, iter(2), nprint,
C   tol(2), seed, printw, ldwork
C
      IF ( mxGetM( PRHS(1) ).NE.1 .OR.
     $     mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), TEMP, 1 )
      TASK = TEMP
      IF ( TASK.LT.0 .OR. TASK.GT.4 ) THEN
         CALL mexErrMsgTxt( 'JOB has 1 : 4 the only admissible values' )
      END IF
C
      IF ( TASK.EQ.1 ) THEN
         INIT = 'L'
      ELSE IF ( TASK.EQ.2 ) THEN
         INIT = 'S'
      ELSE IF ( TASK.EQ.3 ) THEN
         INIT = 'B'
      ELSE
         INIT = 'N'
      END IF
      INIT1 = MOD( TASK, 2 ).EQ.1
      INIT2 = TASK.EQ.2 .OR. TASK.EQ.3
C
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
      NSMP = mxGetM( PRHS(2) )
      M    = mxGetN( PRHS(2) )
      IF ( NSMP.LT.0 .OR. M.LT.0 ) THEN
         CALL mexErrMsgTxt( 'U must be a real matrix' )
      END IF
C
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Y must be a real matrix' )
      END IF
      ITMP = mxGetM( PRHS(3) )
      L    = mxGetN( PRHS(3) )
      IF ( ITMP.LT.0 .OR. L.LT.0 ) THEN
         CALL mexErrMsgTxt( 'Y must be a real matrix' )
      END IF
      IF ( ITMP.NE.NSMP ) THEN
         CALL mexErrMsgTxt
     $       ( 'Y must have the same row dimension as U' )
      END IF
C
      IF ( mxGetM( PRHS(4) ).NE.1 .OR.
     $     mxGetN( PRHS(4) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'NN must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'NN must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), TEMP, 1 )
      NN = TEMP
      IF ( NN.LT.0 ) THEN
         CALL mexErrMsgTxt( 'NN must be a non-negative integer' )
      END IF
C
      IP = 5
      IF ( NRHS.GE.IP .AND. INIT1 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'NOBR (s) must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'NOBR must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         NOBR = TEMP
         IF ( NOBR.LE.0 ) THEN
            CALL mexErrMsgTxt( 'NOBR must be a positive integer' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'N must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'N must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         N = TEMP
         IF ( N.LT.0 .AND. .NOT.INIT1 ) THEN
            CALL mexErrMsgTxt( 'N must be a non-negative integer' )
         END IF
         IP = IP + 1
      ELSE IF ( .NOT. INIT1 ) THEN
         WRITE( TEXT, '('' WIDENTC requires at least '',I4,
     $                  '' input arguments'')' ) IP
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      NS = N
      IF ( NS.LT.0 )
     $   N  = NOBR - 1
      ML    = M + L
      BSN   = NN*( L + 2 ) + 1
      L1    = BSN*L
      L2    = N*( ML + 1 ) + L*M
      LPAR0 = 0
      LX    = L1 + L2
      INI   = 1
      IF ( TASK.EQ.1 ) THEN
         LPAR = L1
      ELSE IF ( TASK.EQ.2 ) THEN
         INI  = L1 + 1
         LPAR = L2
      ELSE IF ( TASK.EQ.4 ) THEN
         LPAR = LX
      END IF
C
      IX = IP
      IF ( NRHS.GE.IP ) THEN
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'X must be a real vector' )
         END IF
         M1 = mxGetM( PRHS(IP) )
         M2 = mxGetN( PRHS(IP) )
         IP = IP + 1
         LPAR0 = M1*M2
      END IF
C
      IALG = 1
      IF ( NRHS.GE.IP ) THEN
C
C   alg
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.NE.1 ) THEN
            WRITE( TEXT,
     $         '('' ALG must be either void or with 1 element'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALG must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, ISIZE )
         IF ( ISIZE.GT.0 ) IALG = INT( TEMP )
         IF ( IALG.LT.1 .OR. IALG.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'IALG has 1 and 2 the only admissible values' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( IALG.EQ.1 ) THEN
         ALG = 'D'
      ELSE
         ALG = 'I'
      END IF
C
      ISTOR = 1
      IF ( NRHS.GE.IP .AND. IALG.EQ.1 ) THEN
C
C   stor
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.NE.1 ) THEN
            WRITE( TEXT,
     $         '('' STOR must be either void or with 1 element'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'STOR must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, ISIZE )
         IF ( ISIZE.GT.0 ) ISTOR = INT( TEMP )
         IF ( ISTOR.LT.1 .OR. ISTOR.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'ISTOR has 1 and 2 the only admissible values' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( ISTOR.EQ.1 ) THEN
         STOR = 'F'
      ELSE
         STOR = 'P'
      END IF
C
      ITMAX1 = 10*LX
      ITMAX2 = 10*LX
      IF ( NRHS.GE.IP ) THEN
C
C   iter
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.NE.2 ) THEN
            WRITE( TEXT, '('' ITER must have either 0 or 2 elements'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ITER must be an integer vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TMP, ISIZE )
         IF ( ISIZE.GT.0 ) ITMAX1 = INT( TMP(1) )
         IF ( ISIZE.GT.1 ) ITMAX2 = INT( TMP(2) )
         IP = IP + 1
      END IF
C
      NPRINT = 0
      IF ( NRHS.GE.IP ) THEN
C
C   nprint
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.NE.1 ) THEN
            WRITE( TEXT,
     $         '('' NPRINT must be either void or with 1 element'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'NPRINT must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, ISIZE )
         IF ( ISIZE.GT.0 ) NPRINT = INT( TEMP )
         IF ( NPRINT.GT.0 )
     $      OPEN( NOUT, file = 'IB03AD.prn', status = 'replace' )
         IP = IP + 1
      END IF
C
      TOL1 = SQRT( DLAMCH( 'Epsilon' ) )
      TOL2 = TOL1
      IF ( NRHS.GE.IP ) THEN
C
C   tol
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.NE.2 ) THEN
            WRITE( TEXT, '('' TOL must have either 0 or 2 elements'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TMP, ISIZE )
         IF ( ISIZE.GT.0 ) TOL1 = TMP(1)
         IF ( ISIZE.GT.1 ) TOL2 = TMP(2)
         IP = IP + 1
      END IF
C
      RSEED(1) = 1998
      RSEED(2) = 1999
      RSEED(3) = 2000
      RSEED(4) = 2001
      IF ( NRHS.GE.IP .AND. INIT2 ) THEN
C
C   seed
C
         ISIZE = mxGetM( PRHS(IP) )*mxGetN( PRHS(IP) )
         IF ( ISIZE.NE.0 .AND. ISIZE.LT.4 ) THEN
            WRITE( TEXT, '('' SEED must have either 0 or 4 elements'')')
            CALL mexErrMsgTxt( TEXT )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'SEED must be an integer or real vector')
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TMP, ISIZE )
         IF ( ISIZE.GT.0 )
     $      CALL DCOPY( ISIZE, TMP, 1, RSEED, 1 )
         IP = IP + 1
      END IF
C
C   printw
C
      PRINTW = .FALSE.
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
         IP = IP + 1
      END IF
C
      LDAC = N + L
      LNOL = L*NOBR - L
      MNO  = M*NOBR
      NSML = NSMP*L
      ISAD = LDAC*( N + M )
      N2   = N*N
C
      JWORK = 0
      IF ( INIT1 ) THEN
C        Workspace for IB01AD.
         JWORK = 2*ML*NOBR*( 2*ML*( NOBR + 1 ) + 3 ) + L*NOBR
C        Workspace for IB01BD.
         IW1 = MAX( 2*LNOL*N + 2*N, LNOL*N + N2 + 7*N, L*NOBR*N +
     $              MAX( LNOL*N + 2*N + ( M + ML )*NOBR + L,
     $                   2*LNOL*N + N2 + 8*N, N + 4*( MNO + N ),
     $                   MNO + 3*N + L ) )
         IF ( M.GT.0 ) THEN
            IW2 = L*NOBR*N + MNO*LDAC*( M*LDAC + 1 ) +
     $            MAX( LDAC**2, 4*M*LDAC + 1 )
         ELSE
            IW2 = 0
         END IF
         JWORK = MAX( JWORK,
     $                ( 2*ML*NOBR )**2 + ISAD + MAX( IW1, IW2 ) )
C        Workspace for IB01CD.
         IW1   = NSML*( N + 1 ) + 2*N + MAX( 2*N2, 4*N )
         IW2   = N*( N + 1 ) + 2*N +
     $           MAX( N*L*( N + 1 ) + 2*N2 + L*N, 4*N )
         JWORK = MAX( JWORK, ISAD + 2 + N*( N + 1 + LDAC + M ) +
     $                MAX( 5*N, 2, MIN( IW1, IW2 ) ) )
C        Workspace for TF01MX.
         JWORK = MAX( JWORK, NSML + ISAD + LDAC + 2*N + M )
C        Workspace for TB01VD.
         JWORK = MAX( JWORK, NSML + ISAD + N +
     $                MAX( 1, N2*L + N*L + N,
     $                     N2 + MAX( N2 + N*MAX( N, L ) +
     $                               6*N + MIN( N, L ), N*M ) ) )
      END IF
C
      IF ( INIT2 ) THEN
C        Workspace for MD03AD (initialization of the nonlinear part).
         IF ( IALG.EQ.1 ) THEN
            IF ( ISTOR.EQ.1 ) THEN
               IW1 = BSN**2
            ELSE
               IW1 = ( BSN*( BSN + 1 ) )/2
            END IF
         ELSE
            IW1 = 3*BSN + NSMP
         END IF
         JWORK = MAX( JWORK, NSML +
     $                MAX( 5, NSMP + 2*BSN + NSMP*BSN +
     $                        MAX( 2*NN + BSN, IW1 ) ) )
         IF ( .NOT.INIT1 ) THEN
C           Workspace for TB01VY.
            JWORK = MAX( JWORK, NSML + LDAC*( 2*N + M ) + 2*N )
C           Workspace for TF01MX.
            IF ( M.GT.0 ) THEN
               IW1 = N + M
            ELSE
               IW1 = 0
            END IF
            JWORK = MAX( JWORK, NSML + ISAD + IW1 + 2*N + L )
         END IF
      END IF
C     Workspace for MD03AD (whole optimization).
      IF ( M.GT.0 ) THEN
         IW1 = LDAC + M
      ELSE
         IW1 = L
      END IF
      IW1 = NSML + MAX( 2*NN, ISAD + 2*N + MAX( N*LDAC, IW1 ) )
      IF ( IALG.EQ.1 ) THEN
         IF ( ISTOR.EQ.1 ) THEN
            IW2 = LX**2
         ELSE
            IW2 = ( LX*( LX + 1 ) )/2
         END IF
      ELSE
         IW2 = 3*LX + NSML
      END IF
      JWORK = MAX( JWORK,
     $             5, NSML + 2*LX + NSML*( BSN + L2 ) +
     $                MAX( IW1 + LX, NSML + IW1, IW2 ) )
      LDWMIN = JWORK
C
      IF ( NRHS.GE.IP ) THEN
C
C     ldwork
C
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'LDWORK must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'LDWORK must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
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
      IF ( NRHS.LT.IP )
     $   LDWORK = LDWMIN
      IF ( INIT1 ) THEN
         IW1 = ML
         IW2 = MAX( MNO + N, M*LDAC )
      ELSE
         IW1 = 0
         IW2 = 0
      END IF
      LIWORK = MAX( 3, IW1, IW2 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( DWORK( LDWORK ), IWORK( LIWORK ), U( LDU, M ),
     $           X( LX ), Y( NSMP, L ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), U, NSMP*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), Y, NSMP*L )
      IF ( LPAR0.GT.0 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(IX) ), X(INI), LPAR0 )
C
C Copy the seed for random numbers generation.
C
      IF ( INIT2 )
     $   CALL DCOPY( 4, RSEED, 1, DWORK, 1 )
      IF ( NS.LT.0 )
     $   N = NS
C
C Do the actual computations.
C
      CALL IB03AD( INIT, ALG, STOR, NOBR, M, L, NSMP, N, NN, ITMAX1,
     $             ITMAX2, NPRINT, U, LDU, Y, LDY, X, LX, TOL1, TOL2,
     $             IWORK, DWORK, LDWORK, IWARN, INFO )
C
C Copy output to MATLAB workspace.
C
      NX = L1 + N*( ML + 1 ) + L*M
      PLHS(1) = mxCreateDoubleMatrix( NX, 1, 0 )
      CALL mxCopyReal8ToPtr( X, mxGetPr( PLHS(1) ), NX )
      IF ( NLHS.GT.1 ) THEN
         IF ( INIT2 ) THEN
            ISIZE = 10
         ELSE
            ISIZE = 5
         END IF
         PLHS(2) = mxCreateDoubleMatrix( ISIZE, 1, 0 )
         CALL mxCopyReal8ToPtr( DWORK, mxGetPr( PLHS(2) ), ISIZE )
      END IF
      IF ( NLHS.GT.2 ) THEN
         PLHS(3)  = mxCreateDoubleMatrix( 2, 1, 0 )
         DWORK(1) = IWORK(1)
         DWORK(2) = IWORK(2)
         CALL mxCopyReal8ToPtr( DWORK, mxGetPr( PLHS(3) ), 2 )
      END IF
      IF ( NLHS.GT.3 .AND. INIT1 ) THEN
         ISIZE   = IWORK(3)
         PLHS(4) = mxCreateDoubleMatrix( ISIZE, 1, 0 )
         CALL mxCopyReal8ToPtr( DWORK(11), mxGetPr( PLHS(4) ), ISIZE )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( DWORK, IWORK, U, X, Y )
C
C Close the printing file.
C
      IF ( NPRINT.GT.0 )
     $   CLOSE( NOUT )
C
C Error and warning handling.
C
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM IB03AD'')'
     $        ) IWARN
         CALL mexPrintf( TEXT )
      END IF
C
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM IB03AD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of WIDENTC ***
      END
