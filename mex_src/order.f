C ORDER.F  - Gateway function for computation of the order of a
C            discrete-time system using SLICOT routine IB01AD.
C
C RELEASE 2.0 of SLICOT System Identification Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [R(,n,sval,rcnd)] = order(meth,alg,jobd,batch,conct,s,Y(,U,tol,
C                             printw,ldwork,R))
C
C Purpose:
C   To preprocess the input-output data for estimating the matrices
C   of a linear time-invariant dynamical system, using Cholesky or
C   QR factorization and subspace identification techniques (MOESP
C   and N4SID), and to find an estimate of the system order. The
C   input-output data can, optionally, be processed sequentially.
C
C Input parameters:
C   meth  - integer option to determine the method to use:
C           = 1 : MOESP method with past inputs and outputs;
C           = 2 : N4SID method.
C   alg   - integer option to determine the algorithm for computing the
C           triangular factor of the concatenated block-Hankel matrices
C           built from the input-output data:
C           = 1 : Cholesky algorithm on the correlation matrix;
C           = 2 : fast QR algorithm;
C           = 3 : standard QR algorithm.
C   jobd  - integer option to specify if the matrices B and D should
C           later be computed using the MOESP approach:
C           = 1 : the matrices B and D should later be computed using
C                 the MOESP approach;
C           = 2 : the matrices B and D should not be computed using
C                 the MOESP approach.
C           This parameter is not relevant for meth = 2.
C   batch - integer option to specify whether or not sequential data
C           processing is to be used, and, for sequential processing,
C           whether or not the current data block is the first block,
C           an intermediate block, or the last block, as follows:
C           = 1:  the first block in sequential data processing;
C           = 2:  an intermediate block in sequential data processing;
C           = 3:  the last block in sequential data processing;
C           = 4:  one block only (non-sequential data processing).
C   conct - integer option to specify whether or not the successive data
C           blocks in sequential data processing belong to a single
C           experiment, as follows:
C           = 1:  the current data block is a continuation of the
C                 previous data block and/or it will be continued by the
C                 next data block;
C           = 2:  there is no connection between the current data block
C                 and the previous and/or the next ones.
C           This parameter is not used if batch = 4.
C   s     - the number of block rows in the input and output block
C           Hankel matrices to be processed.  s > 0.
C   Y     - the t-by-l output-data sequence matrix.  Column j of  Y
C           contains the  t  values of the j-th output component for
C           consecutive time increments.
C   U     - (optional) the t-by-m input-data sequence matrix.  Column j
C           of  U  contains the  t  values of the j-th input component
C           for consecutive time increments.
C           Default:    U = [].
C   tol   - (optional) vector of length 2 containing tolerances:
C           tol(1) - tolerance used for estimating the rank of matrices.
C           If  tol(1) > 0,  then the given value of  tol(1)  is used
C           as a lower bound for the reciprocal condition number;
C           an m-by-n matrix whose estimated condition number is less
C           than  1/tol(1)  is considered to be of full rank.
C           If  tol(1) <= 0,  then a default value m*n*epsilon_machine
C           is used, where epsilon_machine is the relative machine
C           precision.
C           tol(2) - tolerance used for determining an estimate of the
C           system order. If  tol(2) >= 0,  the estimate is indicated
C           by the index of the last singular value greater than or
C           equal to  tol(2).  (Singular values less than  tol(2) are
C           considered as zero.) When  tol(2) = 0,  an internally
C           computed default value,  tol(2) = s*epsilon_machine*sval(1),
C           is used, where  sval(1)  is the maximal singular value, and
C           epsilon_machine the relative machine precision.
C           When  tol(2) < 0,  the estimate is indicated by the index of
C           the singular value that has the largest logarithmic gap to
C           its successor.
C           Default:    tol(1:2) = [0,-1].
C   printw- (optional) switch for printing the warning messages.
C           = 1:  print warning messages;
C           = 0:  do not print warning messages.
C           Default:    printw = 0.
C   ldwork- (optional) the workspace size.
C           Default : computed by the formulas
C           nr = 2*( m + l )*s
C           LDWORK = ( t - 2*s + 3 + 64 )*nr
C           if ( CSIZE > MAX( nr*nr + t*( m + l ) + 16, 2*nr ) ) then
C              LDWORK = MIN( LDWORK, CSIZE - nr*nr - t*( m + l ) - 16 )
C           else
C              LDWORK = MIN( LDWORK, MAX( 2*nr, CSIZE/2 ) )
C           end if
C           LDWORK = MAX( minimum workspace size needed, LDWORK )
C           where CSIZE is the cache size in double precision words.
C           If LDWORK is specified less than the minimum workspace size
C           needed, that minimum value is used instead.
C   R     - (optional) if batch = 2 or 3, the 2*(m+l)*s-by-2*(m+l)*s
C           (upper triangular, if alg <> 2) part of  R  must contain the
C           (upper triangular) matrix  R  computed at the previous call
C           of this mexfile in sequential data processing. If conct = 1,
C           R  has an additional column, also set at the previous call.
C           If alg = 2,  R  has m+l+1 additional columns, set at the
C           previous call.
C           This parameter is not used for batch = 1 or batch = 4.
C
C Output parameter:
C   R     - if batch = 3 or 4, the 2*(m+l)*s-by-2*(m+l)*s part of  R
C           contains the processed upper triangular factor  R  from the
C           QR factorization of the concatenated block-Hankel matrices,
C           and further details needed for computing system matrices.
C           If batch = 1 or 2, then  R  contains intermediate results
C           needed at the next call of this mexfile. If batch = 1 or 2
C           and conct = 1,  R  has an additional column, also set before
C           return. If batch = 1 or 2 and alg = 2,  R  has m+l+1
C           additional columns, set before return.
C   n     - the order of the system.
C   sval  - (optional) the singular values used for estimating the order
C           of the system.
C   rcnd  - (optional) if meth = 2, vector of length 2 containing the
C           reciprocal condition numbers of the matrices involved in
C           rank decisions or least squares solutions.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Oct. 1999.
C
C Revisions:
C   V. Sima, May 2000, July 2000, February 2001, Apr. 2009, Dec. 2012,
C   May 2020.
C
C **********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
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
      DOUBLE PRECISION  TOL1, TOL2
      INTEGER           INFO, IWARN, L, LDR, LDU, LDWORK, LDY, M, N,
     $                  NOBR, NSMP
      CHARACTER         ALG, BATCH, CONCT, CTRL, JOBD, METH
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DWORK(:), R(:,:), SV(:), U(:,:),
     $                                 Y(:,:)
C
C .. Local variables and constant dimension arrays ..
      DOUBLE PRECISION  TEMP, TOL(2)
      INTEGER           CSIZE, IALG, IBCH, ICNC, IJOB, IMTH, IP, ISIZE,
     $                  ITMP, IWRK(3), LDWMIN, LIWORK, NCOL, NR, NRSAVE,
     $                  NSMPMN
      CHARACTER*120     TEXT
      LOGICAL           PRINTW
C
C .. Cache size in double precision words ..
      PARAMETER         ( CSIZE = 64000 )
C
C .. External subroutines ..
      EXTERNAL          DCOPY, DLACPY, IB01AD
C
C ..Intrinsic functions..
      INTRINSIC         INT, MAX, MIN
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ORDER requires at least 7 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ORDER requires at least 1 output arguments' )
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
      IMTH = TEMP
      IF ( IMTH.LT.1 .OR. IMTH.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 or 2 the only admissible values' )
      END IF
C
      IF ( IMTH.EQ.1 ) THEN
         METH = 'M'
      ELSE
         METH = 'N'
      END IF
C
C   alg
C
      IF ( mxGetM( PRHS(2) ).NE.1 .OR. mxGetN( PRHS(2) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'ALG must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'ALG must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), TEMP, 1 )
      IALG = TEMP
      IF ( IALG.LT.1 .OR. IALG.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'ALG has 1, 2 or 3 the only admissible values' )
      END IF
C
      IF ( IALG.EQ.1 ) THEN
         ALG = 'C'
      ELSE IF ( IALG.EQ.2 ) THEN
         ALG = 'F'
      ELSE
         ALG = 'Q'
      END IF
C
C   jobd
C
      IF ( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'JOBD must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'JOBD must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TEMP, 1 )
      IJOB = TEMP
      IF ( IJOB.LT.1 .OR. IJOB.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'JOBD has 1 or 2 the only admissible values' )
      END IF
C
      IF ( IJOB.EQ.1 ) THEN
         JOBD = 'M'
      ELSE
         JOBD = 'N'
      END IF
C
C   batch
C
      IF ( mxGetM( PRHS(4) ).NE.1 .OR. mxGetN( PRHS(4) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'BATCH must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'BATCH must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), TEMP, 1 )
      IBCH = TEMP
      IF ( IBCH.LT.1 .OR. IBCH.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'BATCH has 1, 2, 3 or 4 the only admissible values' )
      END IF
C
      IF ( IBCH.EQ.1 ) THEN
         BATCH = 'F'
      ELSE IF ( IBCH.EQ.2 ) THEN
         BATCH = 'I'
      ELSE IF ( IBCH.EQ.3 ) THEN
         BATCH = 'L'
      ELSE
         BATCH = 'O'
      END IF
C
C   conct
C
      IF ( mxGetM( PRHS(5) ).NE.1 .OR. mxGetN( PRHS(5) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'CONCT must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'CONCT must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TEMP, 1 )
      ICNC = TEMP
      IF ( ICNC.LT.1 .OR. ICNC.GT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONCT has 1 or 2 the only admissible values' )
      END IF
C
      IF ( ICNC.EQ.1 ) THEN
         CONCT = 'C'
      ELSE
         CONCT = 'N'
      END IF
C
C   ctrl is set to 'N'
C
      CTRL = 'N'
C
C   s
C
      IF ( mxGetM( PRHS(6) ).NE.1 .OR. mxGetN( PRHS(6) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 's must be a scalar' )
      END IF
      IF ( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 's must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TEMP, 1 )
      NOBR = TEMP
      IF ( NOBR.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 's must be a positive integer' )
      END IF
C
C   Y(txp), (U(txm))
C
      NSMP = mxGetM( PRHS(7) )
      L    = mxGetN( PRHS(7) )
C
      IF ( L.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no outputs' )
      END IF
C
      IF ( NRHS.GE.8 ) THEN
         M = mxGetN( PRHS(8) )
      ELSE
         M = 0
      END IF
      NR = 2*( M + L )*NOBR
      IF ( IBCH.LE.2 ) THEN
         NSMPMN = 2*NOBR
      ELSE
         NSMPMN = NR + 2*NOBR - 1
      END IF
C
      IF ( NSMP.LT.NSMPMN ) THEN
         WRITE( TEXT, '('' The number of samples should be at least '',
     $                     I10)' ) NSMPMN
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      IF ( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'Y must be a real matrix' )
      END IF
C
      IF ( M.GT.0 ) THEN
         IF ( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'U must be a real matrix' )
         END IF
         IF ( mxGetM( PRHS(8) ).NE.NSMP ) THEN
            CALL mexErrMsgTxt
     $          ( 'U must have the same row dimension as Y' )
         END IF
      END IF
C
C     tol
C
      TOL1 = ZERO
      TOL2 = -ONE
      IF ( NRHS.GE.9 ) THEN
         ISIZE = mxGetM( PRHS(9) ) * mxGetN( PRHS(9) )
         IF ( ISIZE.GT.2 )
     $      CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 2 elements' )
         IF ( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(9) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), TOL, ISIZE )
         IF ( ISIZE.GT.0 ) TOL1 = TOL(1)
         IF ( ISIZE.GT.1 ) TOL2 = TOL(2)
      END IF
C
      PRINTW = .FALSE.
      IP     = 10
      NCOL   = NR
      IF ( IBCH.LT.4 ) THEN
         IF ( ICNC.EQ.1 )
     $      NCOL = NCOL + 1
         IF ( IALG.EQ.2 )
     $      NCOL = NCOL + M + L + 1
      END IF
      NRSAVE  = ( NCOL - NR )*NR
C
C     printw
C
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
C Determine the lenghts of working arrays.
C The quasi-optimal value for LDWORK (assumed block-size 64) is possibly
C modified, so that the sum of the lenghts of DWORK and other needed
C arrays does not exceed the cache size. Specifically, the default value
C for LDWORK is computed using the formulas
C           nr = 2*( m + l )*s
C           LDWORK = ( t - 2*s + 3 + 64 )*nr
C           if ( CSIZE > MAX( nr*nr + t*( m + l ) + 16, 2*nr ) ) then
C              LDWORK = MIN( LDWORK, CSIZE - nr*nr - t*( m + l ) - 16 )
C           else
C              LDWORK = MIN( LDWORK, MAX( 2*nr, CSIZE/2 ) )
C           end if
C           LDWORK = MAX( minimum workspace size needed, LDWORK )
C where CSIZE is the cache size in double precision words.
C If LDWORK is specified, but it is less than the minimum workspace size
C needed, that minimum value is used instead.
C
      LDY = NSMP
      IF ( M.GT.0 ) THEN
         LDU = LDY
      ELSE
         LDU = 1
      END IF
      LDR = NR
      IF ( IMTH.EQ.1 .AND. IJOB.EQ.1 )
     $   LDR = MAX( NR, 3*M*NOBR )
      IF ( IMTH.EQ.2 ) THEN
         LIWORK = ( M + L )*NOBR
      ELSE IF ( IALG.EQ.2) THEN
         LIWORK = M + L
      ELSE
         LIWORK = 1
      END IF
C
C     The next statement is included in order to allow switching
C     from Cholesky or fast QR to QR algorithm.
C
      LDWORK = ( NSMP - 2*NOBR + 3 + 64 )*NR
      IF ( CSIZE .GT. MAX( NR*NR + NSMP*( M + L ) + 16, 2*NR ) ) THEN
         LDWORK = MIN( LDWORK, CSIZE - NR*NR - NSMP*( M + L ) - 16 )
      ELSE
         LDWORK = MIN( LDWORK, MAX( 2*NR, CSIZE/2 ) )
      END IF
      LDWMIN = 2*NR*( NOBR + 1 )
      IF ( IALG.EQ.2 )
     $   LDWMIN = MAX( LDWMIN, 2*NR*( M + L + 1 ) + NR,
     $                 NR*( M + L + 3 ) )
      LDWORK = MAX( LDWORK, LDWMIN )
C
C   ldwork
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'LDWORK must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'LDWORK must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = INT( TEMP )
         IF ( ITMP.LT.LDWMIN ) THEN
            LDWORK = LDWMIN
         ELSE
            LDWORK = ITMP
         END IF
         IP = IP + 1
C
C     R(2*(m+p)*s,2*(m+p)*s(+c)), where
C     c = 1,       if conct = 1 and batch < 4;
C     c = 0,       if conct = 2 or  batch = 4;
C     c = c+m+l+1, if  alg  = 2 and batch < 4.
C
         IF ( NRHS.GE.IP .AND. ( IBCH.EQ.2 .OR. IBCH.EQ.3 ) ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.NR ) THEN
               WRITE( TEXT, '(''R must have '',I5,'' rows'')' ) NR
               CALL mexErrMsgTxt( TEXT )
            END IF
            IF ( mxGetN( PRHS(IP) ).NE.NCOL ) THEN
               WRITE( TEXT, '(''R must have '',I5,'' columns'')' ) NCOL
               CALL mexErrMsgTxt( TEXT )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
               CALL mexErrMsgTxt( 'R must be a real matrix' )
            END IF
         END IF
      END IF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( R( LDR, NCOL ), SV( L*NOBR ), U( LDU, M ), Y( LDY, L )
     $         )
      ALLOCATE ( IWORK( LIWORK ), DWORK( LDWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), Y, NSMP*L )
      IF ( M.GT.0 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), U, NSMP*M )
      END IF
      IF ( NRHS.GE.IP .AND. ( IBCH.EQ.2 .OR. IBCH.EQ.3 ) ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), R, NR*NCOL )
         IF ( IALG.EQ.2 .OR. ICNC.EQ.1 ) THEN
            IF ( LDR.EQ.NR ) THEN
               CALL DCOPY( NRSAVE, R(1,NR+1), 1, DWORK, 1 )
            ELSE
               CALL DCOPY( NRSAVE, R(NR*NR+1,1), 1, DWORK, 1 )
            END IF
         END IF
      END IF
C
C Do the actual computations.
C
      IF ( IBCH.EQ.2 .OR. IBCH.EQ.3 ) THEN
         IWORK(1) = IWRK(1)
         IWORK(2) = IWRK(2)
         IWORK(3) = IWRK(3)
      END IF
      CALL IB01AD( METH, ALG, JOBD, BATCH, CONCT, CTRL, NOBR, M, L,
     $             NSMP, U, LDU, Y, LDY, N, R, LDR, SV, TOL1, TOL2,
     $             IWORK, DWORK, LDWORK, IWARN, INFO )
      IF ( IWARN.NE.0 .AND. PRINTW ) THEN
         WRITE( TEXT, '(''  IWARN = '',I4,'' ON EXIT FROM IB01AD'')' )
     $          IWARN
      END IF
      IF ( INFO.NE.0 ) THEN
         WRITE( TEXT, '(''INFO = '',I4,'' ON EXIT FROM IB01AD'')' )
     $          INFO
         GO TO 20
      END IF
      IF ( IBCH.EQ.1 .OR. IBCH.EQ.2 ) THEN
         IWRK(1) = IWORK(1)
         IWRK(2) = IWORK(2)
         IWRK(3) = IWORK(3)
      END IF
C
C Copy output to MATLAB workspace.
C
      IF ( LDR.GT.NR )
     $   CALL DLACPY( 'Full', NR, NR, R, LDR, R, NR )
      IF ( ( IALG.EQ.2 .OR. ICNC.EQ.1 ) .AND. IBCH.LE.2 ) THEN
         IF ( LDR.EQ.NR ) THEN
            CALL DCOPY( NRSAVE, DWORK, 1, R(1,NR+1), 1 )
         ELSE
            CALL DCOPY( NRSAVE, DWORK, 1, R(NR*NR+1,1), 1 )
         END IF
      END IF
      PLHS(1) = mxCreateDoubleMatrix( NR, NCOL, 0 )
      CALL mxCopyReal8ToPtr( R, mxGetPr( PLHS(1) ), NR*NCOL )
      IF ( IBCH.GT.2 ) THEN
         IF ( NLHS.GT.1 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( 1, 1, 0 )
            TEMP = N
            CALL mxCopyReal8ToPtr( TEMP, mxGetPr( PLHS(2) ), 1 )
            IF ( NLHS.GT.2 ) THEN
               PLHS(3) = mxCreateDoubleMatrix( L*NOBR, 1, 0 )
               CALL mxCopyReal8ToPtr( SV, mxGetPr( PLHS(3) ), L*NOBR )
            END IF
            IF ( IMTH.EQ.2 .AND. NLHS.GT.3 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( 2, 1, 0 )
               CALL mxCopyReal8ToPtr( DWORK(2), mxGetPr( PLHS(4) ), 2 )
            END IF
         END IF
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
   20 CONTINUE
C
      DEALLOCATE( R, SV, U, Y, IWORK, DWORK )
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
C *** Last line of ORDER ***
      END
