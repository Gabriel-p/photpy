C
C**********************************************************************
C
C                              MASTER.FOR
C
C A program for intercomparing star lists obtained from different CCD
C frames to cross-identify stars.  General linear functions of (x,y)
C are used to transform from one set of coordinates to another.  
C Eventual inclusion of color and magnitude terms is possible.  Input 
C is expected in the form of DAOPHOT ASCII data files and output in a 
C variety of forms is available, at the user's discretion.
C
C                   Official DAO version:  1996 June 28
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER MAXDAT, MAXFRM, MAXMTR, MAXCON
      PARAMETER (MAXDAT=60 000 000, MAXFRM=400, 
     .           MAXMTR=10 000 000, MAXCON=20)
C
      CHARACTER*30 FILE(MAXFRM), EXTEND
      CHARACTER*3 CASE
      REAL WORK(MAXDAT)
      REAL CON(MAXCON,MAXFRM), RCOL(MAXFRM), RROW(MAXFRM)
      INTEGER*4 IDMAST(MAXMTR), LASTAR(0:MAXFRM)
      INTEGER IWHICH2(MAXMTR*6),NLINE2(MAXMTR)
      BYTE LINE2(MAXMTR*6)
      INTEGER NOBS2(MAXMTR),INDEX2(MAXMTR)
      INTEGER*4 id1,id2,id3,id4,id5,id6,id7,id8,id9,id10
C
      CHARACTER*255 TEXT
      CHARACTER*30 LFILE
      REAL DUM, CHIMAX
      INTEGER I, J, NN, NNN, IFRM, NFRM, MAXSTR, MAX, ISTAT, A
C
C=======================================================================
C
      CALL FABORT
      DO I=1,MAXFRM
         DO J=1,MAXCON
            CON(J,I) = 0.
         END DO
         CON(3,I) = 1.
         CON(6,I) = 1.
      END DO
C
C Open the file containing the names of the star lists, get the names
C and count the number of input files, and get the estimates of the
C transformation constants.
C
      A=1
      DO A=1,5
C THIS DO RUNS THROUGH ALL THE .mch FILES
      
      IF (A .EQ. 1) THEN
      WRITE (*,*) ' '
      WRITE (*,*) '                 NOW WORKING ON FILTER U'
      ELSE
          IF (A .EQ. 2) THEN
          WRITE (*,*) ' '
          WRITE (*,*) '                 NOW WORKING ON FILTER B'
          ELSE
              IF (A .EQ. 3) THEN
              WRITE (*,*) ' '
              WRITE (*,*) '                 NOW WORKING ON FILTER V'
              ELSE
                  IF (A .EQ. 4) THEN
                  WRITE (*,*) ' '
                  WRITE (*,*) '                 NOW WORKING ON FILTER I'
                  ELSE
                      IF (A .EQ. 5) THEN
                      WRITE (*,*) ' '
             WRITE (*,*) '                 NOW WORKING ON FILE daom.mch'
                      ELSE
                          WRITE (*,*) ' '
                          WRITE (*,*) '       Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF               

      CALL TBLANK
C  890 CALL GETNAM ('File with list of input files:', LFILE)

 890   IF (A .EQ. 1) THEN
      LFILE = 'ufilter.mch'
      ELSE
          IF (A .EQ. 2) THEN
          LFILE = 'bfilter.mch'
          ELSE
              IF (A .EQ. 3) THEN
              LFILE = 'vfilter.mch'
              ELSE
                  IF (A .EQ. 4) THEN
                  LFILE = 'ifilter.mch'
                  ELSE
                      IF (A .EQ. 5) THEN
                      LFILE = 'daom.mch'
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 

      IF ((LFILE .EQ. 'EXIT') .OR. (LFILE .EQ. 'END OF FILE')) 
     .     CALL BYEBYE
      DO I=2,30
         IF (LFILE(I:I) .EQ. ';') THEN
            READ (LFILE(I+1:30),*) CHIMAX
            LFILE = LFILE(1:I-1)
            GO TO 895
         END IF
      END DO
      CHIMAX = 1.E20
  895 CONTINUE
      LFILE = EXTEND(LFILE, CASE('mch'))
      CALL INFILE (1, LFILE, ISTAT)
      IF (ISTAT .NE. 0) THEN 
         CALL STUPID ('Error opening input file '//LFILE)
         LFILE = 'EXIT'
C         GO TO 890
      END IF
C      print*,'(reading ',lfile,')' ! AWS
      IFRM = 0
  900 IFRM = IFRM+1
      CALL RDCHAR (1, TEXT, NN, ISTAT)
      IF (ISTAT .GT. 0) GO TO 950
      IF (IFRM .GT. MAXFRM) THEN
         WRITE (TEXT,6) MAXFRM
    6    FORMAT ('Too many input files.  Maximum allowed is', I4)
         CALL STUPID (TEXT)
         CALL CLFILE (1)
         CALL OOPS
      END IF
      IF (NN .GT. 200) THEN
         READ (TEXT,*,IOSTAT=ISTAT) FILE(IFRM), 
     .        (CON(I,IFRM), I=1,6), DUM, DUM,
     .        (CON(I,IFRM), I=7,20)
      ELSE IF (NN .GT. 120) THEN
         READ (TEXT,*,IOSTAT=ISTAT) FILE(IFRM), 
     .        (CON(I,IFRM), I=1,6), DUM, DUM,
     .        (CON(I,IFRM), I=7,12)
      ELSE
         READ (TEXT,*) FILE(IFRM), 
     .        (CON(I,IFRM), I=1,6)
      END IF
      IF (ISTAT .NE. 0) THEN
         CALL STUPID ('Error reading input from file.')
         LFILE = 'EXIT'
         CALL CLFILE (1)
C         GO TO 890
      END IF
      GO TO 900
  950 NFRM = IFRM-1
C
C Now read all the input data into the array WORK.  WORK must be able 
C to contain six 4-byte numbers per stellar observation: 
C
C     X, Y  ---  coordinates of the star
C     M, S  ---  the magnitude and its standard error (sigma)
C      CHI  ---  the CHI value from the profile fit
C    SHARP  ---  the SHARP value from the profile fit
C
C Thus, the absolute maximum number of stars there is room for is the 
C size of WORK divided by 6.
C
      MAXSTR = MAXDAT/6

      print*,'  Currently dimensioned for a maximum file size of',maxmtr
      print*,'  and a maximum total number of stars of', maxstr

      CALL RDFRMS (FILE, CHIMAX, NFRM, MAXMTR, WORK, MAXSTR, IDMAST, 
     .     LASTAR, RCOL, RROW)
C
C At this point, we have the following input arrays stored at
C the following locations:
C
C       X(ISTAR)  ---  WORK(1)
C       Y(ISTAR)  ---  WORK(MAXDAT/6+1)
C       M(ISTAR)  ---  WORK(2*MAXDAT/6+1)
C       S(ISTAR)  ---  WORK(3*MAXDAT/6+1)
C     CHI(ISTAR)  ---  WORK(4*MAXDAT/6+1)
C   SHARP(ISTAR)  ---  WORK(5*MAXDAT/6+1)
C
C RDFRMS has changed the value of the variable MAXSTR from MAXDAT/6 to 
C the number of stars in the largest single file.  Now we must pack 
C these data as tightly as possible in the upper part of the array 
C WORK, so that we can use the remainder of the array for accumulating 
C and storing the master list.
C
      CALL REPACK (WORK, MAXDAT/6, LASTAR(NFRM) )
C
C If NTOT = LASTAR(NFRM) is the total number of stellar observations in 
C all frames, then the input arrays are now stored at the following 
C locations:
C
C       X(ISTAR)  ---  WORK(1)
C       Y(ISTAR)  ---  WORK(NTOT+1)
C       M(ISTAR)  ---  WORK(2*NTOT+1)
C       S(ISTAR)  ---  WORK(3*NTOT+1)
C     CHI(ISTAR)  ---  WORK(4*NTOT+1)
C   SHARP(ISTAR)  ---  WORK(5*NTOT+1)
C
C This leaves us with a total of MAXDAT-6*NTOT storage locations 
C starting at WORK(6*NTOT+1) to divide up among the following arrays:
C
C   IWHICH of dimension NFRM*MAX, where MAX is the maximum number of 
C          stars for which there will be room for in the provisional 
C          master list
C
C and either
C
C   LINE of dimension MAXSTR*32 (to hold 127 characters per star, in the
C        single frame containing the largest number of stars = MAXSTR),
C        plus
C
C   NLINE of dimension MAXSTR
C   
C
C or
C
C   the seven arrays XM, YM, MM, SM, RMIN, NOBS, and INDEX, each of 
C        dimension MAX
C
C whichever is larger (because the array LINE will not be needed until
C after we are completely finished with XM, YM, MM, ..., they may share
C the same storage locations).  We are also restricted to having MAX no
C greater than MAXMTR by the existence of array IDMAST.  Hence
C
      MAX = MIN0 ( MAXMTR                                  ,
     .            ( MAXDAT-6*LASTAR(NFRM)-33*MAXSTR)/NFRM   ,
     .            ( MAXDAT-6*LASTAR(NFRM))/(NFRM+7)         )
C
C MAX now contains the largest permissible number of stars ever to
C appear in the provisional master list (the same star recurring in
C several frames counts once).  
C
      WRITE (TEXT,7) MAX
    7 FORMAT (' The maximum permissible number of stars in the ',
     .     'master list: ', I9)
      CALL OVRWRT(TEXT(1:73), 3)
C
C Now we can actually run MASTER, now that we can specify where in 
C WORK the storage for each of the arrays is to start.
C
      NN=6*LASTAR(NFRM)+1                    ! End of input data + 1
      NNN=NN+MAX*NFRM                        ! End of IWHICH storage + 1
C
C Don't use Real work array for IWHICH
      IF (NFRM*MAX.GT.MAXMTR*6) THEN
          WRITE(*,*) NFRM,MAX
          CALL STUPID ('IWHICH2 IS NOT BIG ENOUGH')
          CALL OOPS
      END IF
C Don't use Real work array for nline
      IF ((MAXSTR+1).GT.MAXMTR) THEN
          CALL STUPID ('NLINE2 IS NOT BIG ENOUGH')
          CALL OOPS
      END IF
C Don't use Real work array for line
      IF ((MAXSTR*128).GT.MAXMTR*6) THEN
          CALL STUPID ('LINE2 IS NOT BIG ENOUGH')
          CALL OOPS
      END IF
C Don't use Real work array for nobs or index
      IF (MAX.GT.MAXMTR) THEN
          CALL STUPID ('NOBS2 OR INDEX2 IS NOT BIG ENOUGH')
          CALL OOPS
      END IF
      id1=LASTAR(NFRM)+1
      id2=2*LASTAR(NFRM)+1
      id3=3*LASTAR(NFRM)+1
      id4=4*LASTAR(NFRM)+1
      id5=5*LASTAR(NFRM)+1
      id6=NNN
      id7=NNN+MAX
      id8=NNN+2*MAX
      id9=NNN+3*MAX
      id10=NNN+4*MAX
      CALL   MASTER (LFILE, FILE, NFRM, CON, MAX, IDMAST, LASTAR,
     .  RCOL,RROW,WORK,WORK(id1),WORK(id2),WORK(id3),WORK(id4),
     .  WORK(id5),IWHICH2,NLINE2,LINE2,WORK(id6),WORK(id7),
     .  WORK(id8),WORK(id9),WORK(id10),NOBS2,INDEX2)
C
C SUBROUTINE  MASTER (LFILE, FILE, NFRM, CON, MAX, IDMAST, LASTAR,
C
c     .     RCOL, RROW,
c     .     WORK(1), WORK(LASTAR(NFRM)+1), WORK(2*LASTAR(NFRM)+1), 
c     .     WORK(3*LASTAR(NFRM)+1), WORK(4*LASTAR(NFRM)+1),  ! Input data
c     .     WORK(5*LASTAR(NFRM)+1),
C
C          X, Y, M, S, CHI, SHARP                           ! Input data
C
C     .     WORK(NN), 
c     .     IWHICH2, 
C
C           IWHICH,
C
c     .     WORK(NNN), WORK(NNN+MAXSTR),                  ! MAXSTR arrays
c     .     NLINE2, LINE2,                  ! MAXSTR arrays
C
C            NLINE,        LINE,                         ! MAXSTR arrays
C
c     .     WORK(NNN), WORK(NNN+MAX), WORK(NNN+2*MAX),
c     .     WORK(NNN+3*MAX), WORK(NNN+4*MAX), WORK(NNN+5*MAX),
c     .     WORK(NNN+6*MAX)    )                             ! MAX arrays
c     .     WORK(NNN+3*MAX), WORK(NNN+4*MAX), NOBS2, INDEX2) ! MAX arrays
C
C          XM, YM, MM, SM, RMIN, NOBS, INDEX)               ! MAX arrays
C
   
      END DO
      
      CALL BYEBYE
      END!
C
C###############################################################################
C
      SUBROUTINE  MASTER  (LFILE, FILE, NFRM, CON, MAX, IDMAST, LASTAR,
     .     RCOL, RROW,
     .     X, Y, M, S, CHI, SHARP,                       ! Input data
C
     .     IWHICH,                                       ! IWHICH
C
     .     NLINE, LINE,                                  ! MAXSTR arrays
C
     .     XM, YM, MM, SM, RMIN, NOBS, INDEX)            ! MAX arrays
C
c      IMPLICIT NONE
C
C Input data:
C 
C    MAX is the maximum number of stars that may ever appear in the
C        provisional master list (the same star reappearing in several 
C        frames counts as one.)
C
C   NFRM is the number of different data frames (input files).
C
C Parameters:
C
C MAXFRM is the maximum permitted number of frames.
C
C MAXCON is the maximum number of constants, or coefficients, in the 
C        pair of transformation equations.
C
C MINY, MAXY are the minimum and maximum y-coordinates (referred to the
C        standard data frame) for which space will be reserved in some
C        auxiliary arrays.  Larger and smaller y-values are permitted
C        by the program, however.
C
C NOTE that NLINE and LINE share storage with XM, YM, ...
C
      INTEGER MAX, MAXFRM, MINY, MAXY, MAXCON
      PARAMETER (MAXFRM=400, MAXCON=20, MINY=-11 384, MAXY=15 000)
C
C Arrays and functions.
C
      CHARACTER*132 TEXT
      CHARACTER*81 HEAD(2)
      CHARACTER*30 FILE(*), OFILE
      CHARACTER*4 CASE
      DOUBLE PRECISION C(MAXCON,MAXCON), V(MAXCON), DCON(MAXCON)
      DOUBLE PRECISION CLAMP(MAXCON,MAXFRM)!, DABS
      REAL X(*), Y(*), M(*), S(*), CHI(*), SHARP(*), RCOL(*), RROW(*)
      REAL XM(*), YM(*), MM(*), SM(*), RMIN(*)
      INTEGER IRMIN(10 000 000)
      REAL RMIN2(10 000 000)
      REAL OLD(MAXCON,MAXFRM), CON(MAXCON,*), TERM(MAXCON), 
     .     NOC(MAXCON,MAXFRM)
      REAL DATUM(2,MAXFRM), DMAG(MAXFRM), FLIP(MAXFRM)
      CHARACTER(LEN=10) DATUM2(2,MAXFRM)
      REAL MAG90(MAXFRM), MAG95(MAXFRM), SMAG(MAXFRM), WWW(MAXFRM)
      INTEGER IDMAST(MAX), LASTAR(0:MAXFRM), NOK(MAXFRM)
      INTEGER IWHICH(MAX,*), NOBS(*), INDEX(*), NLINE(*)
      INTEGER LAST1(MINY:MAXY), LAST2(MINY:MAXY), NTERM(MAXFRM)
      INTEGER IBNRY, INT, NINT, ICLOSE, LENGTH
      BYTE LINE(128,*)

C
C Variables and function names.
C
      CHARACTER*30 LFILE, SWITCH
      CHARACTER*1 ANSWER
      REAL AMIN1, AMAX1, PCTILE
      REAL XX, YY, RADIUS, RADSQ, SIGMAX, ZZ, A, B, SUM, SUMWT
      REAL WT, DIFF, SUMCHI, SUMSHP, SADX, SADY, SAD, SIGN, RR, FRAC
      REAL RCLAMP, ROLD, XS, YS, SCALE, DIAG, SOFT!, DTRMNT
      INTEGER I, J, K, L, IFIRST, ISTAR, IMASTR, IFRM, NFRM, NL
      INTEGER IBUMP, III, IFLAG, MODE, NEW, NMASTR, LEGIT
      INTEGER LIMIT, ISTAT, MINFRM, ENOUGH, NITER, MIN, MINUSE
      DATA DMAG/MAXFRM*99.999/, NOK/MAXFRM*0/, MIN/10/
      DATA SOFT /1.E-6/
C
C=======================================================================
C
      CALL TBLANK
C      CALL GETDAT ('Minimum number, minimum fraction, enough frames:', 
C     .     RMIN, 3)
     
     
      IF (LFILE .EQ. 'ufilter.mch') THEN
      RMIN(1) = 1
      RMIN(2) = 0
      RMIN(3) = 1
      ELSE
          IF (LFILE .EQ. 'bfilter.mch') THEN
      RMIN(1) = 1
      RMIN(2) = 0
      RMIN(3) = 1
          ELSE
              IF (LFILE .EQ. 'vfilter.mch') THEN
      RMIN(1) = 1
      RMIN(2) = 0
      RMIN(3) = 1
              ELSE
                  IF (LFILE .EQ. 'ifilter.mch') THEN
      RMIN(1) = 1
      RMIN(2) = 0
      RMIN(3) = 1
                  ELSE
                      IF (LFILE .EQ. 'daom.mch') THEN
      RMIN(1) = 2
      RMIN(2) = 0
      RMIN(3) = 2
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 
       
      WRITE (*,*) 'Minimum number, minimum fraction, enough frames: ',
     .       RMIN(1), RMIN(2), RMIN(3)
       
     
      IF (RMIN(1) .LE. 0.5) CALL BYEBYE
      MINFRM = NINT(RMIN(1))
      MINUSE = MAX0(2, MINFRM)
      FRAC = AMIN1(1., AMAX1(0., RMIN(2)))
      ENOUGH = NINT(RMIN(3))
C
C      CALL GETDAT ('Maximum sigma:', SIGMAX, 1)

      IF (LFILE .EQ. 'ufilter.mch') THEN
      SIGMAX = 99.
      ELSE
          IF (LFILE .EQ. 'bfilter.mch') THEN
          SIGMAX = 99.
          ELSE
              IF (LFILE .EQ. 'vfilter.mch') THEN
              SIGMAX = 99.
              ELSE
                  IF (LFILE .EQ. 'ifilter.mch') THEN
                  SIGMAX = 99.
                  ELSE
                      IF (LFILE .EQ. 'daom.mch') THEN
                      SIGMAX = 99.
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 
       
       WRITE (*,*) 'Maximum sigma: ',SIGMAX

      IF (SIGMAX .LE. 0) CALL BYEBYE
      SIGMAX=SIGMAX**2
C
   80 CALL TBLANK
      WRITE (6,81)
   81 FORMAT (15X, 'Desired degrees of freedom ---'//
     .        18X, '2:  Translations only'/
     .        18X, '4:  Translations, rotation, and scale'/
     .        18X, '6:  Full, six-constant linear transformation')
      CALL TBLANK
C      CALL GETDAT ('Your choice:', WT, 1)
      
      IF (LFILE .EQ. 'ufilter.mch') THEN
      WT = 2
      ELSE
          IF (LFILE .EQ. 'bfilter.mch') THEN
          WT = 2
          ELSE
              IF (LFILE .EQ. 'vfilter.mch') THEN
              WT = 2
              ELSE
                  IF (LFILE .EQ. 'ifilter.mch') THEN
                  WT = 2
                  ELSE
                      IF (LFILE .EQ. 'daom.mch') THEN
                      WT = 2
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 
       
       WRITE (*,*) 'Desired degrees of freedom: ',WT
      
      IF ((WT .LE. 0.5) .OR. (WT .GT. 20.5)) CALL BYEBYE
      MODE=NINT(WT)
      IF ((MODE .NE. 2) .AND. (MODE .NE. 4) .AND. (MODE .NE. 6) .AND.
     .     (MODE .NE. 12) .AND. (MODE .NE. 20)) THEN
           CALL STUPID ('Invalid response.')
           GO TO 80
      END IF
C
      DO I=1,NFRM
         XX = (RCOL(I)+1.)/2.
         YY = (RROW(I)+1.)/2.
         CALL TRFM (CON(1,I), RCOL(I), RROW(I), FLIP(I), 
     .        MODE, XX, YY, MAG90(I), MAG95(I))
         ZZ = ABS(MAG90(I)-MAG90(1))+ABS(MAG95(I)-MAG95(1))
         WWW(I) = 1./(2.+ZZ/AMAX1(512., AMIN1(RCOL(I), RCOL(1), 
     .        RROW(I), RROW(1))))
      END DO
C
      WRITE (*,*) ' '
      CALL GETDAT ('Critical match-up radius:', RADIUS, 1)
      IF (RADIUS .LE. 0) GO TO 80
C
   
      DO I=1,NFRM
         L = LASTAR(I-1)
         XX = 0.
         YY = 0.
         DO ISTAR = L+1, LASTAR(I)
            XX = AMAX1(XX, ABS(X(ISTAR)))
            YY = AMAX1(YY, ABS(Y(ISTAR)))
         END DO
         SCALE = SQRT(ABS(CON(3,I)*CON(6,I)) + ABS(CON(4,I)*CON(5,I)))
         DIAG = RADIUS/(2.*SCALE)
c        A = DIAG /1.5
c        B = A /1.5
         IF (MODE .LT. 20) THEN
            DO J=MAX0(7,MODE+1),20
               CON(J,I) = 0.0
            END DO
         END IF
         DMAG(I) = 99.999
         DO J=1,MODE
            OLD(J,I) = 0.0
         END DO
         clamp(1,i) = diag
         clamp(2,i) = diag
         clamp(3,i) = 4.*diag/xx
         clamp(4,i) = 4.*diag/xx
         clamp(5,i) = 4.*diag/yy
         clamp(6,i) = 4.*diag/yy
         do j=7,12
            clamp(j,i) = amin1(1.,diag)
         end do
         do j=13,20
            clamp(j,i) = amin1(1.,diag)
         end do
         IF (MODE .EQ. 2) THEN
            CON(3,I)=1.
            CON(4,I)=0.
            CON(5,I)=0.
            CON(6,I)=1.
            DO J=7,20
               CON(J,I)=0.
            END DO
         ELSE IF (MODE .EQ. 4) THEN
            FLIP(I)=SIGN(1.,
     .           CON(3,I)*CON(6,I)-CON(4,I)*CON(5,I))
            CON(3,I) = 0.5*(CON(3,I) + FLIP(I)*CON(6,I))
            CON(4,I) = 0.5*(CON(4,I) - FLIP(I)*CON(5,I))
            CON(5,I) = -FLIP(I)*CON(4,I)
            CON(6,I) = FLIP(I)*CON(3,I)
         END IF
      END DO
      DMAG(1)=0.
C
C-----------------------------------------------------------------------
C
C Initialize the master list, by setting it equal to the starlist from
C the first frame:
C
      DO I=1,LASTAR(1)
         XM(I)=X(I)
         YM(I)=Y(I)
         MM(I)=M(I)
         NOBS(I)=1
         IWHICH(I,1)=I
      END DO
      NOC(1,1) = 0.
      NOC(2,1) = 0.
      NOC(3,1) = 1.
      NOC(4,1) = 0.
      NOC(5,1) = 0.
      NOC(6,1) = 1.
      NMASTR=LASTAR(1)
      IFIRST=2
C
C Estimate the magnitude limit for each frame.  This is
C provisionally a sloping line from 1.0 at the 90-th percentile
C to 0.5 at the 95-th percentile and extrapolated to zero.
C
      DO IFRM=1,NFRM
         I = LASTAR(IFRM-1)
         DO ISTAR = I+1, LASTAR(IFRM)
            J = ISTAR-I
            RMIN(J) = M(ISTAR)
         END DO
         I = LASTAR(IFRM) - I
         MAG90(IFRM) = PCTILE(RMIN, I, NINT(0.90*I))
         MAG95(IFRM) = PCTILE(RMIN, I, NINT(0.95*I))
         IF (MAG90(IFRM) .GE. MAG95(IFRM)) THEN
            MAG95(IFRM) = 1.E6
         ELSE
            MAG95(IFRM) = 0.5/(MAG95(IFRM) - MAG90(IFRM))
         END IF
      END DO
C
      RCLAMP = 0.1
      NITER = 0
 2000 NITER = NITER+1
      MIN = MAX0(0,MIN-1)
      RADSQ=RADIUS**2
      CALL TBLANK
      CALL OVRWRT (' ', 4)
C
C Sort master list by y coordinate.
C
      CALL QUICK (XM, NMASTR, INDEX)
      CALL RECTFY (YM, NMASTR, INDEX, RMIN)
      CALL RECTFY (MM, NMASTR, INDEX, RMIN)
C      IF (IFIRST .EQ. 2) CALL RECTFY (IWHICH, NMASTR, INDEX, RMIN)
      IF (IFIRST .EQ. 2) CALL IRECTFY (IWHICH, NMASTR, INDEX, IRMIN)
C
      CALL QUICK (YM, NMASTR, INDEX)
      CALL RECTFY (XM, NMASTR, INDEX, RMIN)
      CALL RECTFY (MM, NMASTR, INDEX, RMIN)
      IF (IFIRST .EQ. 2) CALL IRECTFY (IWHICH, NMASTR, INDEX, IRMIN)
C      IF (IFIRST .EQ. 2) CALL RECTFY (IWHICH, NMASTR, INDEX, RMIN)
C
      DO I=MINY,MAXY
         LAST1(I)=0
      END DO
      DO I=1,NMASTR
         J=MAX0(MINY, MIN0(MAXY, INT(YM(I))))
         LAST1(J)=I
      END DO
      DO I=MINY+1,MAXY
         IF (LAST1(I) .LE. 0) LAST1(I)=LAST1(I-1)
      END DO
      NEW=NMASTR
      LIMIT=NEW
C
C Go through each individual star list, trying to match each star in the
C list with a star on the current master list.  If it cannot be matched
C with any star on the master list, try matching it with any stars in
C previous frames which have been appended to the master list during
C this iteration.  If it still cannot be matched, add it to the end
C of the master list.
C
C NMASTR is the number of stars in the master list at the end of the
C        previous iteration.
C  LIMIT is the number of stars in the master list after checking the
C        previous frame, including stars added provisionally from
C        previous frames during the current iteration, but not
C        including stars added from this frame during this iteration
C    NEW will be the number of stars in the current master list,
C        including stars from this frame added during this iteration
C
      DO 2190 IFRM=IFIRST,NFRM
C
C If there have been any stars added to the end of the master list during 
C this iteration, sort the provisional new stars by y coordinate.
C
      WRITE (TEXT,2) LIMIT, 'Starting  ', IFRM, FILE(IFRM)
    2 FORMAT (I8, 2X, A, I3, 3X, A30)
      CALL OVRWRT (TEXT(1:54), 2)
      IF (NEW .GT. NMASTR) THEN
C
C Perform the quicksort twice, once by X coordinate and once by
C Y.  This is much, much faster than trying to sort by Y when the
C list is already nearly sorted by Y.
C
         CALL  QUICK (XM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1))
         CALL RECTFY (YM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1), RMIN)
         CALL RECTFY (MM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1), RMIN)
         DO I=1,IFRM-1
C            CALL RECTFY (IWHICH(NMASTR+1,I), NEW-NMASTR, 
C     .           INDEX(NMASTR+1), RMIN)
            CALL IRECTFY (IWHICH(NMASTR+1,I), NEW-NMASTR, 
     .           INDEX(NMASTR+1), IRMIN)
         END DO
C
         CALL  QUICK (YM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1))
         CALL RECTFY (XM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1), RMIN)
         CALL RECTFY (MM(NMASTR+1), NEW-NMASTR, INDEX(NMASTR+1), RMIN)
         DO I=1,IFRM-1
C            CALL RECTFY (IWHICH(NMASTR+1,I), NEW-NMASTR, 
C     .           INDEX(NMASTR+1), RMIN)
            CALL IRECTFY (IWHICH(NMASTR+1,I), NEW-NMASTR, 
     .           INDEX(NMASTR+1), IRMIN)
         END DO
C
         DO I=MINY,MAXY
            LAST2(I)=NMASTR
         END DO
         DO I=NMASTR+1,NEW
            J=MAX0(MINY, MIN0(MAXY, INT(YM(I))))
            LAST2(J)=I
         END DO
         DO I=MINY+1,MAXY
            IF (LAST2(I) .LE. NMASTR) LAST2(I)=LAST2(I-1)
         END DO
      END IF
C
      LIMIT=NEW
C     WRITE (TEXT,2) LIMIT, 'Checking ', IFRM, FILE(IFRM)
C     CALL OVRWRT (TEXT(1:54), 2)
C
      DO I=1,LIMIT
         RMIN(I)=RADSQ
         IWHICH(I,IFRM)=-1
      END DO
C
C For each star in this star list, identify the star in the master star
C list it lies closest to.
C
      DO 2180 ISTAR=LASTAR(IFRM-1)+1, LASTAR(IFRM)
C
C Transform the star's coordinates from the system of frame IFRM to the 
C system of the master list.
C
         CALL TRFM (CON(1,IFRM), RCOL(IFRM), RROW(IFRM), FLIP(IFRM), 
     .        MODE, X(ISTAR), Y(ISTAR), XX, YY)
         III=ISTAR
C
C See which star in the master list from the previous iteration (stars
C 1 through NMASTR) is closest to this star.  If no unmatched star in 
C the previous master list is within the critical radius, zero will be 
C returned.
C
 2155    IMASTR=ICLOSE(XX, YY, RADIUS, RADSQ, XM, YM, 
     .        NOBS, LAST1, 1, RMIN)
C
         IF (IMASTR .EQ. 0) THEN
C
C This star is not in the master list produced during the last 
C iteration, so see whether it matches any of the stars in other 
C frames which were provisionally added to the master list during this 
C iteration.
C
            IF (LIMIT .GT. NMASTR) THEN
               IMASTR = ICLOSE (XX, YY, RADIUS, RADSQ, XM, YM, 
     .              NOBS, LAST2, NMASTR+1, RMIN)
               IF (IMASTR .NE. 0) GO TO 2160
            END IF
C 
C This star doesn't match anything here, either, so add it to the end 
C of the list (BUT ONLY IF THERE IS ROOM).
C
            IF (NEW .LT. MAX) THEN
               NEW=NEW+1
               XM(NEW)=XX
               YM(NEW)=YY
               NOBS(NEW)=1
               RMIN(NEW)=RADSQ
               IF (DMAG(IFRM) .LT. 90.) THEN
                  MM(NEW)=M(III)+DMAG(IFRM)
               ELSE
                  MM(NEW)=99.999
               END IF
               IF (IFRM .GT. 1) THEN
                  DO I=1,IFRM-1
                     IWHICH(NEW,I)=-1
                  END DO
               END IF
               IWHICH(NEW,IFRM)=III
               GO TO 2180
            END IF
         END IF
C
C This star matches a star either in the previous master list or
C in the stars added to the master list during this iteration...
C 
 2160    IF (IWHICH(IMASTR,IFRM) .GT. 0) THEN
C 
C ... but another star had already been matched with this star in the
C master list.  We know that this new star is a better match, so attempt
C to re-match the previous star.
C
            IBUMP=IWHICH(IMASTR,IFRM)
            IWHICH(IMASTR,IFRM)=III
            CALL TRFM (CON(1,IFRM), RCOL(IFRM), RROW(IFRM), FLIP(IFRM), 
     .           MODE, X(IBUMP), Y(IBUMP), XX, YY)
            III=IBUMP
            GO TO 2155
         ELSE
C
C ... and no other star had previously been matched with this star in 
C the master list.
C
            IWHICH(IMASTR,IFRM)=III
         END IF
 2180 CONTINUE
C
C Update the mean coordinates of any stars added to the master list
C during this iteration.
C
      IF (LIMIT .GT. NMASTR) THEN
         DO I=NMASTR+1,LIMIT
            J = IWHICH(I,IFRM)
            IF (J .GT. 0) THEN
               A = REAL(NOBS(I))
               B = A+1.
               CALL TRFM (CON(1,IFRM), RCOL(IFRM), RROW(IFRM), 
     .              FLIP(IFRM), MODE, X(J), Y(J), XX, YY)
               XM(I) = (A*XM(I) + XX)/B
               YM(I) = (A*YM(I) + YY)/B
               NOBS(I) = NOBS(I) + 1
            END IF
         END DO
      END IF
 2190 CONTINUE
C
C AT THIS POINT...
C
C IWHICH(IMASTR,IFRM) tells which star in the starlist for the IFRM-th
C frame corresponds to the IMASTR-th star in the master list.
C
      DO IMASTR=1,NMASTR
         NOBS(IMASTR)=0
      END DO
C
      DO IFRM=1,NFRM
         DO IMASTR=1,NMASTR
            IF (IWHICH(IMASTR,IFRM) .GT. 0) 
     .           NOBS(IMASTR)=NOBS(IMASTR)+1
         END DO
      END DO
C
      do imastr=nmastr+1,limit
CXYZ  DO IMASTR=1,NMASTR
         IF (NOBS(IMASTR) .LE. 1) THEN
            L = IWHICH(IMASTR,1)
            IF (L .GT. 0) THEN
               MM(IMASTR) = M(L)
            ELSE
               MM(IMASTR) = 99.999
            END IF
         END IF
      END DO
C
C NOW...
C
C NOBS(IMASTR) holds the number of frames in which the star appears.
C
C Looping over the various input frames, use the current values of 
C the standard (x,y) coordinates and the table of cross-identifications
C to improve the transformation coefficients.  Employ only those stars
C which were in the master list at the beginning of this iteration,
C not those just added.
C
C Begin loop over frames...
C
C     RDSQ01 = RADSQ*0.1
      LEGIT = 0
      DO 2850 IFRM=IFIRST,NFRM
      WRITE (TEXT,9) IFRM, FILE(IFRM)
    9 FORMAT ('  Correcting transformations for', I5, 3X, A30)
      CALL OVRWRT (TEXT(1:73), 2)
      DO I=1,MODE
         V(I)=0.0D0
         DO J=1,MODE
            C(I,J)=0.0D0
         END DO
      END DO
C
      SADX=0.0
      SADY=0.0
      SAD=0.0
C
C Begin loop over stars...
C
      III = 0
      DO 2800 IMASTR=1,NMASTR
         IF (NOBS(IMASTR) .LT. MINUSE) GO TO 2800
         IF (MM(IMASTR) .GT. 90.) GO TO 2800
         L=IWHICH(IMASTR,IFRM)
         IF (L .GT. 0) THEN
            WT = SOFT * (NOBS(IMASTR)-1.) / (NFRM*S(L))
C
C Double the weight of a star in the master frame.
C
            IF (IWHICH(IMASTR,1) .GT. 0) WT = 2.*WT
cC
cC Further reduce the weights of stars near the magnitude limit.
cC
c           WT=WT*AMAX1(0., AMIN1(1., 1.-MAG95(IFRM)*
c    .           (M(L)-MAG90(IFRM))))
            CALL TRFM (CON(1,IFRM), RCOL(IFRM), RROW(IFRM), 
     .           FLIP(IFRM), MODE, X(L), Y(L), XX, YY)
C
C Reduce the weight of a star near the edge of the critical circle.
C
            RR = (XX-XM(IMASTR))**2 + (YY-YM(IMASTR))**2
            WT = WT/(1.+(RR/RADSQ)**2)
            III = III+1
C
C X residual.
C
            DIFF=XM(IMASTR)-XX
            SADX=SADX+ABS(DIFF)
            SAD=SAD+1.
            IF (IFRM .EQ. 1) GO TO 2790
            TERM(1)=1.
            TERM(2)=0.
            IF (MODE .EQ. 4) THEN
              TERM(3)=X(L)
              TERM(4)=-FLIP(IFRM)*Y(L)
            ELSE IF (MODE .GE. 6) THEN
              TERM(3)=X(L)
              TERM(4)=0.
              TERM(5)=Y(L)
              TERM(6)=0.
              IF (MODE .GE. 7) THEN
                XS = 2.*(X(L)-1.)/(RCOL(IFRM)-1.)-1.
                YS = 2.*(Y(L)-1.)/(RROW(IFRM)-1.)-1.
                TERM(7)=1.5*XS**2-0.5
                TERM(8)=0.
                TERM(9)=XS*YS
                TERM(10)=0.
                TERM(11)=1.5*YS**2-0.5
                TERM(12)=0.
                IF (MODE .GE. 13) THEN
                  TERM(13)=XS*TERM(7)
                  TERM(14)=0.
                  TERM(15)=YS*TERM(7)
                  TERM(16)=0.
                  TERM(17)=XS*TERM(11)
                  TERM(18)=0.
                  TERM(19)=YS*TERM(11)
                  TERM(20)=0.
                END IF
              END IF
            END IF
            DO I=1,MODE
               V(I)=V(I)+WT*TERM(I)*DIFF
               DO J=1,MODE
                  C(I,J)=C(I,J)+WT*TERM(I)*TERM(J)
               END DO
            END DO
C
C Y residual.
C
 2790       DIFF=YM(IMASTR)-YY
            SADY=SADY+ABS(DIFF)
            IF (IFRM .EQ. 1) GO TO 2800
            TERM(1)=0.
            TERM(2)=1.
            IF (MODE .EQ. 4) THEN
              TERM(3)=FLIP(IFRM)*Y(L)
              TERM(4)=X(L)
            ELSE IF (MODE .GE. 6) THEN
              TERM(3)=0.
              TERM(4)=X(L)
              TERM(5)=0.
              TERM(6)=Y(L)
              IF (MODE .GE. 7) THEN
                TERM(8)=TERM(7)
                TERM(7) = 0.
                TERM(10)=TERM(9)
                TERM(9) = 0.
                TERM(12)=TERM(11)
                TERM(11) = 0.
                IF (MODE .GE. 13) THEN
                  TERM(14) = TERM(13)
                  TERM(13) = 0.
                  TERM(16) = TERM(15)
                  TERM(15) = 0.
                  TERM(18) = TERM(17)
                  TERM(17) = 0.
                  TERM(20) = TERM(19)
                  TERM(19) = 0.
                END IF
              END IF
            END IF
            DO I=1,MODE
               V(I)=V(I)+WT*TERM(I)*DIFF
               DO J=1,MODE
                  C(I,J)=C(I,J)+WT*TERM(I)*TERM(J)
               END DO
            END DO
         END IF
 2800 CONTINUE
C
      IF (IFRM .EQ. 1) THEN
         III = 0
         GO TO 2840
      END IF
C
      IF (III .LT. MIN) THEN
         NTERM(IFRM) = 0
         GO TO 2850
      ELSE IF (III-MIN .LE. 1) THEN
         III = 2
      ELSE IF ((III-MIN .EQ. 2) .AND. (MODE .NE. 4)) THEN
         III = 2
      ELSE IF ((III-MIN .EQ. 2) .AND. (MODE .EQ. 4)) THEN
         III = 4
      ELSE IF ((III-MIN .LE. 4) .AND. (MODE .GE. 6)) THEN
         III = 6
      ELSE IF ((III-MIN .LE. 8) .AND. (MODE .GE. 12)) THEN
         III = 12
      ELSE
         III = MIN0(MODE,20)
      END IF
      IF (NOK(IFRM) .LT. 1) THEN
         III = MIN0(III, 2)
      ELSE IF (NOK(IFRM) .LT. 2) THEN
         III = MIN0(III, 6)
      ELSE IF (NOK(IFRM) .LT. 3) THEN
         III = MIN0(III, 12)
      ELSE 
         III = MIN0(III, 20)
      END IF
      LEGIT = MAX0(LEGIT, III)
C
C     do i=1,iii
C        c(i,i) = c(i,i) / clamp(i,ifrm)
C     end do
c
      CALL DINVRS (C, MAXCON, III, IFLAG)
      IF (IFLAG .NE. 0) THEN
         CALL STUPID ('Matrix error.')
         GO TO 2850
      END IF
      NOK(IFRM) = NOK(IFRM)+1
      CALL DVMUL (C, MAXCON, III, V, DCON)
      DO I=1,III
         A = SNGL(DCON(I))
         B = A*OLD(I,IFRM)
         IF (B .LT. 0.) THEN
            CLAMP(I,IFRM) = 0.5*CLAMP(I,IFRM)
         ELSE IF (B .GT. 0.0) THEN
C            CLAMP(I,IFRM) = AMIN1(1., 1.05*CLAMP(I,IFRM))
            CLAMP(I,IFRM) = DMIN1(1.D0, 1.05*CLAMP(I,IFRM))
         END IF
         CON(I,IFRM) = CON(I,IFRM)+DCON(I)
cxyz .        /(1.+DABS(DCON(I))/CLAMP(I,IFRM))
         OLD(I,IFRM) = A
      END DO
      IF (III .EQ. 4) THEN
         CON(5,I) = -FLIP(I)*CON(4,I)
         CON(6,I) = FLIP(I)*CON(3,I)
      END IF
C
 2840 IF (SAD .LE. 0.5) GO TO 2850
      DATUM(1,IFRM) = 1.2533*SADX/SAD
      DATUM(2,IFRM) = 1.2533*SADY/SAD
      NTERM(IFRM) = III
 2850 CONTINUE
C
C End of loop over input frames.
C
cC If the coordinate system of the master list has drifted away from that
cC of frame 1, twist it back.
cC
c      DTRMNT = (CON(3,1)*CON(6,1) - CON(4,1)*CON(5,1))
c      DO IFRM=NFRM,1,-1
c         DCON(1) = (CON(6,1)*(CON(1,IFRM)-CON(1,1))-
c     .        CON(5,1)*(CON(2,IFRM)-CON(2,1))) / DTRMNT
c         DCON(2) = (CON(3,1)*(CON(2,IFRM)-CON(2,1))-
c     .        CON(4,1)*(CON(1,IFRM)-CON(1,1))) / DTRMNT
c         DCON(3) = (CON(6,1)*CON(3,IFRM)-CON(5,1)*CON(4,IFRM)) / DTRMNT
c         DCON(4) = (CON(3,1)*CON(4,IFRM)-CON(4,1)*CON(3,IFRM)) / DTRMNT
c         DCON(5) = (CON(6,1)*CON(5,IFRM)-CON(5,1)*CON(6,IFRM)) / DTRMNT
c         DCON(6) = (CON(3,1)*CON(6,IFRM)-CON(4,1)*CON(5,IFRM)) / DTRMNT
c         DO I=1,MIN(MODE,6)
c            CON(I,IFRM) = DCON(I)
c         END DO
c      END DO
c      IF (MODE .GE. 7) THEN
c         DO I=7,MODE
c            CON(I,1) = 0.
c         END DO
c      END IF
C
C Now, again looping over the various input frames, compute the additive
C magnitude corrections to the system of the master list.  First, erase
C all memory of the master-list magnitudes of stars in an insufficient
C number of frames, so they won't be used in the computation.
C
      K = MAX0(2,MINFRM)
      DO IMASTR=1,NMASTR
         IF (NOBS(IMASTR) .LT. K) THEN
            MM(IMASTR) = 99.999
         END IF
      END DO
C
C Begin loop over frames...
C
      DO 2950 IFRM=IFIRST,NFRM
      L=0
C
C Begin loop over stars...
C
      DO 2900 IMASTR=1,NMASTR
         IF (MM(IMASTR) .GT. 90.) GO TO 2900
         J=IWHICH(IMASTR,IFRM)
         IF (J .GT. 0) THEN
            L=L+1
            RMIN(L)=MM(IMASTR)-M(J)
            SM(L) = 1. / S(J)
C
C Weight of observations with really huge errors goes as
C 1/sigma**4.
C
            sm(l) = sm(l) * 1.e-2/(1.e-2+s(j))
C
C Double the weight of a star in the master image.
C
            IF (IWHICH(IMASTR,1) .GT. 0) SM(L) = 2.*SM(L)
C
C Further reduce the weights of stars near the magnitude limit.
C
            SM(L)=SM(L)*AMAX1(0., AMIN1(1., 1.-MAG95(IFRM)*
     .           (M(J)-MAG90(IFRM))))
         END IF
 2900 CONTINUE
C
      IF (L .GE. MIN) THEN
         CALL QUICK (RMIN, L, INDEX)
C         CALL RECTFY (SM, L, INDEX, NOBS)
         CALL RECTFY (SM, L, INDEX, RMIN2)
         DO I=2,L
            SM(I) = SM(I)+SM(I-1)
         END DO
C
C Locate the magnitude difference of median cumulative weight, using a 
C binary search.
C
         J = IBNRY(SM, L, 0.5)
         DMAG(IFRM) = RMIN(J)
C        dmag(ifrm) = 0.   !xyz
         J = IBNRY(SM, L, 0.3085)
         A = RMIN(J)
         J = IBNRY(SM, L, 0.6915)
         B = RMIN(J)
         SUM=DMAG(IFRM)-DMAG(1)
         SMAG(IFRM)=B-A
      ELSE
         DMAG(IFRM) = 99.999
         SUM = 99.99
         SMAG(IFRM) = 0.
      END IF
C
      IF (MODE .EQ. 2) THEN
         WRITE (TEXT,3) (DATUM(I,IFRM), I=1,2), 
     .        (NINT(CON(I,IFRM)), I=1,2), SUM, 
     .        SMAG(IFRM), L, NTERM(IFRM), 
     .        FILE(IFRM)(1:LENGTH(FILE(IFRM)))
    3    FORMAT (1X, 2F5.2, 2I7, F7.2, F6.3, I7, I3, 1X, A)
         CALL OVRWRT (TEXT(1:71), 3)
      ELSE IF (MODE .EQ. 4) THEN
         WRITE (TEXT,1) (DATUM(I,IFRM), I=1,2), 
     .       (NINT(CON(I,IFRM)), I=1,2), (CON(I,IFRM), I=3,4), SUM, 
     .       SMAG(IFRM), L, NTERM(IFRM), 
     .        FILE(IFRM)(1:LENGTH(FILE(IFRM)))
    1    FORMAT (1X, 2F5.2, 2I7, 2F7.3, F7.2, F6.3, I7, I3, 1X, A)
         CALL OVRWRT(TEXT(1:87), 3)
      ELSE
         WRITE (TEXT,4) (DATUM(I,IFRM), I=1,2), 
     .        (NINT(CON(I,IFRM)), I=1,2), (CON(I,IFRM), I=3,6), SUM, 
     .        SMAG(IFRM), L, NTERM(IFRM), 
     .        FILE(IFRM)(1:LENGTH(FILE(IFRM)))
    4    FORMAT (1X, 2F5.2, 2I7, 4F7.3, F7.2, F6.3, I7, I3, 1X, A)
         CALL OVRWRT(TEXT(1:LENGTH(TEXT)), 3)
      END IF
      SMAG(IFRM) = SMAG(IFRM)**2
 2950 CONTINUE
C
C Correct magnitude offsets
C
      DO IFRM=NFRM,2,-1
         IF (DMAG(IFRM) .LE. 90.) DMAG(IFRM)=DMAG(IFRM)-DMAG(1)
      END DO
      DMAG(1) = 0.
C
      IFIRST=1
      CALL TBLANK
C
C Now redetermine mean positions and magnitudes of stars on master list.
C
      DO IMASTR=1,NEW
         MM(IMASTR) = 0.
         SM(IMASTR) = 0.
         RMIN(IMASTR) = 0.
         NOBS(IMASTR) = 0
      END DO
C
C Update the mean coordinates in the system of the master
C list.
C
      DO 2970 IFRM=1,NFRM
C        IF ((IFRM .GT. 1) .AND. (NTERM(IFRM) .LT. MODE)) GO TO 2970
C
C Weight goes as the inverse square of the image scale.
C
         WT = 1./(ABS(CON(3,IFRM)*CON(6,IFRM)) + 
     .        ABS(CON(4,IFRM)*CON(5,IFRM)))
c        IF (IFRM .NE. 1) WT = WT*WWW(IFRM)*
c    .        (REAL(NTERM(IFRM))/MODE)**4
         DO IMASTR=1,NEW
            J=IWHICH(IMASTR,IFRM)
            IF (J .GT. 0) THEN
               CALL TRFM (CON(1,IFRM), RCOL(IFRM), RROW(IFRM),
     .              FLIP(IFRM), MODE, X(J), Y(J), XX, YY)
C
C Reduce the weight of a position near the edge of the
C critical circle.  
C
               RR=(XX-XM(IMASTR))**2 + (YY-YM(IMASTR))**2
               RR = WT*(1.-RR/RADSQ) / S(J)
               IF (RR .GT. 0.) THEN
                  NOBS(IMASTR) = NOBS(IMASTR)+1
                  if (nterm(ifrm) .ge. mode) then
                     MM(IMASTR)=MM(IMASTR)+RR*(XX-XM(IMASTR))
                     SM(IMASTR)=SM(IMASTR)+RR*(YY-YM(IMASTR))
                     RMIN(IMASTR)=RMIN(IMASTR)+RR
                  end if
               END IF
            END IF
         END DO
         IF (NTERM(IFRM) .GT. 0) WWW(IFRM) = (3.*WWW(IFRM)+1.)/4.
 2970 CONTINUE
C
      DO IMASTR=1,NEW
         IF (RMIN(IMASTR) .GT. 0.) THEN
            RR = MM(IMASTR)/RMIN(IMASTR)
            XM(IMASTR)=XM(IMASTR)+RR/(1.+(RR/RCLAMP)**2)
            RR = SM(IMASTR)/RMIN(IMASTR)
            YM(IMASTR)=YM(IMASTR)+RR/(1.+(RR/RCLAMP)**2)
            MM(IMASTR)=0.
            SM(IMASTR)=0.
c        ELSE
c           NOBS(IMASTR) = 0
         END IF
      END DO
C
C Now update the magnitudes on the system of the master list.
C
      DO IMASTR=1,NEW
         DO IFRM=1,NFRM
            J=IWHICH(IMASTR,IFRM)
            IF ((J .GT. 0) .AND. (DMAG(IFRM) .LT. 90.)) THEN
               RR = 1. / S(J)
               MM(IMASTR)=MM(IMASTR)+(M(J)+DMAG(IFRM))*RR
               SM(IMASTR)=SM(IMASTR)+RR
            END IF
         END DO
         IF (SM(IMASTR) .GT. 0.) THEN
            MM(IMASTR)=MM(IMASTR)/SM(IMASTR)
            SM(IMASTR)=1./SM(IMASTR)
         ELSE IF (IWHICH(IMASTR,1) .GT. 0) THEN
            MM(IMASTR) = M(IWHICH(IMASTR,1))
            SM(IMASTR) = S(IWHICH(IMASTR,1))
         ELSE
            MM(IMASTR)=99.999
         END IF
      END DO
C
C Invert the transformation equations.
C
      DO IFRM=2,NFRM
         DO I=1,MAXCON
            NOC(I,IFRM) = CON(I,IFRM)
         END DO
         CALL INVERT (NOC(1,IFRM), RCOL(IFRM), RROW(IFRM), 
     .         RCOL(1), RROW(1), FLIP, MODE)
      END DO
C
      IMASTR=0
 2500 CONTINUE
C
C If the LAST star on the master list is insignificant, then eliminate 
C it.
C
      IF ((SM(NEW) .GT. SIGMAX) .OR. 
     .     (NOBS(NEW) .LT. MINFRM)) THEN
         NEW=NEW-1
         IF (NEW .LE. 0) THEN
            CALL STUPID ('ERROR --- No stars left.')
            CALL BYEBYE
         END IF
         GO TO 2500
      ELSE
         IF (NOBS(NEW) .GE. ENOUGH) GO TO 2600
         A = 0.
         B = 0.
C
C Count the frames in which this star COULD reasonably have
C appeared.  That means its predicted coordinates fall
C within the area of the frame [ 1.0 < x < real(NCOL),
C 1.0 < y < real(NROW) ]. Furthermore, the likelihood
C that a star SHOULD be found in a frame is a function of
C its predicted instrumental magnitude in that frame:
C the "weight" of that frame is 1.0 if the predicted 
C magnitude is brighter than the frame's 90-th percentile, and
C drops linearly to 0.0, passing through 0.5 at the
C frame's 95-th percentile.
C
         DO 2550 IFRM=1,NFRM
c           IF (NTERM(IFRM) .LT. LEGIT) GO TO 2550
            CALL TRFM (NOC(1,IFRM), RCOL(1), RROW(1), 
     .           FLIP(IFRM), MODE, XM(NEW), YM(NEW), XX, YY)
            IF ((XX .GE. 1.) .AND. (XX .LE. RCOL(IFRM)) .AND.
     .          (YY .GE. 1.) .AND. (YY .LE. RROW(IFRM)) .AND.
     .          (NTERM(IFRM) .GE. 2)) THEN
              WT = AMAX1(0., AMIN1(1., 1.-MAG95(IFRM)*
     .             (MM(NEW)-DMAG(IFRM)-MAG90(IFRM))))
              A = A+FRAC*WT
              IF (IWHICH(NEW,IFRM) .GT. 0) B = B+WT
            END IF
 2550    CONTINUE
         IF (B .GE. A) GO TO 2600
C
         NEW=NEW-1
         IF (NEW .LE. 0) THEN
            CALL STUPID ('ERROR --- No stars left.')
            CALL BYEBYE
         END IF
         GO TO 2500
      END IF
C
C Now check the next star down from the top of the list.  If this
C star is insignificant, then overwrite it with the last star in
C the list (which we already know, from above, to be significant).
C
 2600 IMASTR=IMASTR+1
      IF (IMASTR .GE. NEW) GO TO 2700
      IF ((SM(IMASTR) .GT. SIGMAX) .OR.
     .     (NOBS(IMASTR) .LT. MINFRM)) THEN
C
C  This star is insignificant.  Overwrite it with the last star
C  in the list, decrement NEW by one, and go back up to
C  statement 2500 to test the significance of the new last star.
C
         XM(IMASTR)=XM(NEW)
         YM(IMASTR)=YM(NEW)
         MM(IMASTR)=MM(NEW)
         SM(IMASTR)=SM(NEW)
         NOBS(IMASTR)=NOBS(NEW)
         DO IFRM=1,NFRM
            IWHICH(IMASTR,IFRM)=IWHICH(NEW,IFRM)
         END DO
         NEW=NEW-1
         GO TO 2500
      ELSE
         IF (NOBS(IMASTR) .GE. ENOUGH) GO TO 2600
         A = 0.
         B = 0.
C
C Count the frames in which this star COULD reasonably have
C appeared.
C
         DO 2650 IFRM=1,NFRM
c           IF (NTERM(IFRM) .LT. LEGIT) GO TO 2650
            CALL TRFM (NOC(1,IFRM), RCOL(1), RROW(1), 
     .           FLIP(IFRM), MODE, XM(IMASTR), YM(IMASTR), XX, YY)
            IF ((XX .GE. 1.) .AND. (XX .LE. RCOL(IFRM)) .AND.
     .          (YY .GE. 1.) .AND. (YY .LE. RROW(IFRM)) .AND.
     .          (NTERM(IFRM) .GE. 2)) THEN
              WT = AMAX1(0., AMIN1(1., 1.-MAG95(IFRM)*
     .             (MM(IMASTR)-DMAG(IFRM)-MAG90(IFRM))))
              A = A+FRAC*WT
              IF (IWHICH(IMASTR,IFRM) .GT. 0) B = B+WT
            END IF
 2650    CONTINUE
         IF (B .GE. A) GO TO 2600
C
C  This star is insignificant.  Overwrite it with the last star
C  in the list, decrement NEW by one, and go back up to
C  statement 2500 to test the significance of the new last star.
C
         XM(IMASTR)=XM(NEW)
         YM(IMASTR)=YM(NEW)
         MM(IMASTR)=MM(NEW)
         SM(IMASTR)=SM(NEW)
         NOBS(IMASTR)=NOBS(NEW)
         DO IFRM=1,NFRM
            IWHICH(IMASTR,IFRM)=IWHICH(NEW,IFRM)
         END DO
         NEW=NEW-1
         GO TO 2500
      END IF
 2700 NMASTR = NEW
C
C Now we have a pure list of verified significant cross-identifications.
C
      WRITE (6,603) NMASTR
  603 FORMAT (2X, I6, 2X, 'stars.', 33(' '))
      CALL TBLANK
      ROLD = RADIUS
      CALL GETDAT ('New match-up radius (0 to exit):', RADIUS, 1)
      IF (RADIUS .LT. -1.E10) CALL BYEBYE
      IF (RADIUS .GE. ROLD) THEN
         RCLAMP = AMAX1(1.E-6,0.7*RCLAMP)
      ELSE
         RCLAMP = AMIN1(1.,1.5*RCLAMP)
      END IF
      IF (RADIUS .GT. 0.) GO TO 2000
C
      WRITE (6,61)
   61 FORMAT (/5X,
     .   ' Transformations are in the sense  STANDARD = fn(OBSERVED).'
     .   /)
C
C Get ready to put out results.
C
c     ofile = 'junk.dat'
c     call outfil (1, ofile, istat)
c     do imastr=1,nmastr
c        if (iwhich(imastr,1) .eq. 0) write (1,*) imastr, xm(imastr),
c    .        ym(imastr), mm(imastr), sm(imastr), nobs(imastr)
c     end do
c     call clfile (1)
c      CALL GETYN ('Assign new star IDs?', ANSWER)

       ANSWER = 'Y'
       WRITE (*,*) 'Assign new star IDs?: ',ANSWER

      IF (ANSWER .EQ. 'Y') THEN
         CALL QUICK (MM, NMASTR, INDEX)
C
C Sort stars by y coordinate.
C
         CALL RECTFY (XM, NMASTR, INDEX, RMIN)
         CALL RECTFY (YM, NMASTR, INDEX, RMIN)
         CALL RECTFY (SM, NMASTR, INDEX, RMIN)
C         CALL RECTFY (NOBS, NMASTR, INDEX, RMIN)
         CALL IRECTFY (NOBS, NMASTR, INDEX, IRMIN)
         DO IFRM=1,NFRM
C            CALL RECTFY (IWHICH(1,IFRM), NMASTR, INDEX, RMIN)
            CALL IRECTFY (IWHICH(1,IFRM), NMASTR, INDEX, IRMIN)
         END DO
         DO IMASTR=1,NMASTR
            INDEX(IMASTR) = IMASTR
         END DO
      ELSE
         III = 200000
         DO IMASTR=1,NMASTR
            J = IWHICH(IMASTR,1)
            IF (J .LE. 0) THEN
               III = III+1
               INDEX(IMASTR) = III
            ELSE
               INDEX(IMASTR) = IDMAST(J)
            END IF
         END DO
      END IF
C
      DO IMASTR=1,NMASTR
         IDMAST(IMASTR) = INDEX(IMASTR)
      END DO
C
      CALL TBLANK
      CALL INQUIR ('Now, do you want...', 45)
      CALL TBLANK
c      CALL GETYN ('A file with mean magnitudes and scatter?', ANSWER)

      IF (LFILE .EQ. 'ufilter.mch') THEN
      ANSWER = 'Y'
      ELSE
          IF (LFILE .EQ. 'bfilter.mch') THEN
          ANSWER = 'Y'
          ELSE
              IF (LFILE .EQ. 'vfilter.mch') THEN
              ANSWER = 'Y'
              ELSE
                  IF (LFILE .EQ. 'ifilter.mch') THEN
                  ANSWER = 'Y'
                  ELSE
                      IF (LFILE .EQ. 'daom.mch') THEN
                      ANSWER = 'N'
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 
       
       WRITE (*,*) 'A file with mean magnitudes and scatter?: ',ANSWER

      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         OFILE=SWITCH(LFILE, CASE('.mag'))
         CALL GETFIL ('Output file name:', OFILE, 1, 'NEW')
         CALL INFILE (2, FILE(1), ISTAT)
         CALL RDCHAR (2, HEAD(1), K, ISTAT)
         CALL RDCHAR (2, HEAD(2), L, ISTAT)
         CALL CLFILE (2)
         IF ((K .GE. 3) .AND. (HEAD(1)(2:3) .EQ. 'NL')) THEN
            WRITE (1,115) HEAD(1)(1:K)
            IF (L .GE. 1) THEN
               WRITE (HEAD(2)(1:3),115) '  1'
               WRITE (1,115) HEAD(2)(1:L)
            ELSE
               WRITE (1,*)
            END IF
            WRITE (1,*)
         END IF
C
         DO IMASTR=1,NMASTR
            IF (NOBS(IMASTR) .GE. 2) THEN
C
C Determine robust weighted mean magnitudes of stars on master list.
C
               MM(IMASTR) = EXP(0.921034*(15.-MM(IMASTR)))
               SUMCHI=0.
               SUMSHP=0.
               SUMWT=0.
               DO IFRM=1,NFRM
                  J=IWHICH(IMASTR,IFRM)
                  IF (J .GT. 0) THEN
                     M(J) = EXP(0.921034*(15.-M(J)-DMAG(IFRM)))
                     S(J) = S(J)*(0.921034*M(J))**2
                     WT = 1. / S(J)
                     SUMWT=SUMWT+WT
                     SUMCHI=SUMCHI+CHI(J)*WT
                     SUMSHP=SUMSHP+SHARP(J)*WT
                  END IF
               END DO
               SUMCHI=SUMCHI/SUMWT
               SUMSHP=SUMSHP/SUMWT
 3010          CONTINUE
               SUM=0.
               SUMWT=0.
               DO IFRM=1,NFRM
                  J=IWHICH(IMASTR,IFRM)
                  IF (J .GT. 0) THEN
                     B=M(J)-MM(IMASTR)
C
C Give half weight to 2-sigma residuals.
C
                     WT=(1. / S(J)) * (4. / (4.+(B**2 / S(J))))
                     SUM=SUM+B*WT
                     SUMWT=SUMWT+WT
                  END IF
               END DO
               B=SUM/SUMWT
               MM(IMASTR)=MM(IMASTR)+B
               IF (ABS(B) .GE. 0.00005*MM(IMASTR)) GO TO 3010
C
               A=0.
               RR=0.
               SUM=0.
               SUMWT=0.
               DO IFRM=1,NFRM
                  J=IWHICH(IMASTR,IFRM)
                  IF (J .GT. 0) THEN
                     B=M(J)-MM(IMASTR)
                     WT = 1. / S(J)
                     A=A+WT*ABS(B)
                     IF (B .GT. 0.) RR = RR+1.
                     SUM=SUM+WT*(ABS(B)/SQRT(S(J)))
                     SUMWT=SUMWT+WT
                     S(J)=S(J)/(0.921034*M(J))**2
                     M(J)=15.-1.0857362*ALOG(M(J))-DMAG(IFRM)
                  END IF
               END DO
               B=SQRT(NOBS(IMASTR)/(NOBS(IMASTR)-1.))
               A=1.2533*A/SUMWT                                          ! sigma(1 obs.)
               RR=RR/REAL(NOBS(IMASTR))
               A=1.0857362*A*B/MM(IMASTR)
               SUM=1.2533*B*SUM/SUMWT                                      ! CHI
               MM(IMASTR)=15.-1.0857362*ALOG(MM(IMASTR))
            ELSE
               SUM=0.
               SUMWT=0.
               SUMCHI=0.
               SUMSHP=0.
               RR=1.
               DO IFRM=1,NFRM
                  J=IWHICH(IMASTR,IFRM)
                  IF (J .GT. 0) THEN
                     A=M(J)+DMAG(IFRM)
                     WT=1./S(J)
                     SUM=SUM+A*WT
                     SUMCHI=SUMCHI+CHI(J)*WT
                     SUMSHP=SUMSHP+SHARP(J)*WT
                     SUMWT=SUMWT+WT
                  END IF
               END DO
               MM(IMASTR)=SUM/SUMWT
               SUMCHI=SUMCHI/SUMWT
               SUMSHP=SUMSHP/SUMWT
               A=0.
               SUM=0.
            END IF
C
C  1. Star ID
C  2. X coordinate, in system of master frame
C  3. Y coordinate, in system of master frame
C  4. Robust intensity-weighted mean instrumental magnitude
C  5. Standard error of the mean magnitude, based on individual frame sigma's
C  6. External standard error of one measurement, based on frame-to-frame repeatability
C  7. Number of frames in which the star appeared
C  8. Weighted average CHI value
C  9. Weighted average SHARP value
C 10. Variability index: ratio of external error to internal error
C 11. Blunder index: fraction of residuals that are positive
C
            IF ((XM(IMASTR) .LE. -999.9995) .OR. (YM(IMASTR) .LE.
     .           -999.9995) .OR. (XM(IMASTR) .GE. 9999.995) .OR.
     .           (YM(IMASTR) .GE. 9999.995)) THEN
               WRITE (1,110) IDMAST(IMASTR), XM(IMASTR), YM(IMASTR),
     .              MM(IMASTR), SQRT(SM(IMASTR)), A, NOBS(IMASTR), 
     .              SUMCHI, SUMSHP, SUM, RR
  110          FORMAT (I6, 2F9.2, F9.3, 2F9.4, I8, '.', 
     .              2F9.3, 2F9.2, F9.3)
            ELSE
               WRITE (1,1105) IDMAST(IMASTR), XM(IMASTR), YM(IMASTR),
     .              MM(IMASTR), SQRT(SM(IMASTR)), A, NOBS(IMASTR),
     .              SUMCHI, SUMSHP, SUM, RR
 1105          FORMAT (I6, 3F9.3, 2F9.4, I8, '.', 2F9.3, F9.2, F9.3)
            END IF
         END DO
         CALL CLFILE (1)
      END IF
C
C      CALL GETYN 
C     .     ('A file with corrected magnitudes and errors?', ANSWER)

      ANSWER = 'Y'
      WRITE (*,*) 'A file with corrected magnitudes and errors?: ',
     .      ANSWER

      IF (ANSWER .EQ. 'E') THEN 
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         OFILE=SWITCH(LFILE, CASE('.cor'))
         CALL GETFIL ('Output file name:', OFILE, 1, 'NEW')
         CALL INFILE (2, FILE(1), ISTAT)
         CALL RDCHAR (2, HEAD(1), K, ISTAT)
         CALL RDCHAR (2, HEAD(2), L, ISTAT)
         CALL CLFILE (2)
         IF ((K .GE. 3) .AND. (HEAD(1)(2:3) .EQ. 'NL')) THEN
            WRITE (1,115) HEAD(1)(1:K)
            IF (L .GE. 1) THEN
               WRITE (HEAD(2)(1:3),115) '  1'
               WRITE (1,115) HEAD(2)(1:L)
            ELSE
               WRITE (1,*)
            END IF
            WRITE (1,*)
         END IF
C
         DO IMASTR=1,NMASTR
            SUMCHI=0.
            SUMSHP=0.
            SUMWT=0.
            DO IFRM=1,NFRM
               J=IWHICH(IMASTR,IFRM)
               IF (J .GT. 0) THEN
                  DATUM(1,IFRM)=M(J)+DMAG(IFRM)
                  DATUM(2,IFRM)=AMIN1
     .                 (9.9999,SQRT(AMAX1(0., S(J)-SOFT)))
                  A=EXP(0.921034*(15.-DATUM(1,IFRM)))
                  B=S(J)*(0.921034*A)**2
                  WT = 1./B
                  SUMWT=SUMWT+WT
                  SUMCHI=SUMCHI+CHI(J)*WT
                  SUMSHP=SUMSHP+SHARP(J)*WT
               ELSE
                  DATUM(1,IFRM)=99.9999
                  DATUM(2,IFRM)=9.9999
               END IF
            END DO
            SUMCHI = SUMCHI/SUMWT
            SUMSHP = SUMSHP/SUMWT
            IF ((XM(IMASTR) .LE. -999.9995) .OR. (YM(IMASTR) .LE.
     .           -999.9995) .OR. (XM(IMASTR) .GE. 9999.995) .OR.
     .           (YM(IMASTR) .GE. 9999.995)) THEN
               WRITE (1,111) IDMAST(IMASTR),
     .               XM(IMASTR), YM(IMASTR), 
     .               ((DATUM(J,I), J=1,2), I=1,NFRM),
     .               SUMCHI, SUMSHP
  111          FORMAT (I6, 2F9.2, 12F9.4:/ (24X, 12F9.4))
            ELSE
               WRITE (1,1115) IDMAST(IMASTR),
     .               XM(IMASTR), YM(IMASTR), 
     .               ((DATUM(J,I), J=1,2), I=1,NFRM),
     .               SUMCHI, SUMSHP
 1115          FORMAT (I6, 2F9.3, 12F9.4:/ (24X, 12F9.4))
 1116          FORMAT (I6, 2F9.3, 8A, 4F9.4:/ (24X, 12F9.4))
            END IF
         END DO
         CALL CLFILE (1)
      END IF
C
C      CALL GETYN ('A file with raw magnitudes and errors?', ANSWER)
      
      IF (LFILE .EQ. 'ufilter.mch') THEN
      ANSWER = 'N'
      ELSE
          IF (LFILE .EQ. 'bfilter.mch') THEN
          ANSWER = 'N'
          ELSE
              IF (LFILE .EQ. 'vfilter.mch') THEN
              ANSWER = 'N'
              ELSE
                  IF (LFILE .EQ. 'ifilter.mch') THEN
                  ANSWER = 'N'
                  ELSE
                      IF (LFILE .EQ. 'daom.mch') THEN
                      ANSWER = 'Y'
                      ELSE
                          WRITE (*,*) 'Unknown error. Check code'
                      END IF
                  END IF
              END IF
          END IF
       END IF 
       
       WRITE (*,*) 'A file with raw magnitudes and errors?: ',ANSWER      
      
      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         OFILE=SWITCH(LFILE, CASE('.raw'))
         CALL GETFIL ('Output file name:', OFILE, 1, 'NEW')
         CALL INFILE (2, FILE(1), ISTAT)
         CALL RDCHAR (2, HEAD(1), K, ISTAT)
         CALL RDCHAR (2, HEAD(2), L, ISTAT)
         CALL CLFILE (2)
         IF ((K .GE. 3) .AND. (HEAD(1)(2:3) .EQ. 'NL')) THEN
            WRITE (1,115) HEAD(1)(1:K)
            IF (L .GE. 1) THEN
               WRITE (HEAD(2)(1:3),115) '  1'
               WRITE (1,115) HEAD(2)(1:L)
            ELSE
               WRITE (1,*)
            END IF
            WRITE (1,*)
         END IF
C
         DO IMASTR=1,NMASTR
            SUMCHI=0.
            SUMSHP=0.
            DO IFRM=1,NFRM
               J=IWHICH(IMASTR,IFRM)
               IF (J .GT. 0) THEN
                  DATUM(1,IFRM)=M(J)
                  DATUM(2,IFRM)=AMIN1
     .                 (9.9999,SQRT(AMAX1(0., S(J)-SOFT)))
      WRITE (DATUM2(1,IFRM), *) (DATUM(1,IFRM))
      WRITE (DATUM2(2,IFRM), *) (DATUM(2,IFRM))  
    
                  SUMCHI=SUMCHI+CHI(J)
                  SUMSHP=SUMSHP+SHARP(J)
               ELSE
                  DATUM(1,IFRM)=99.9999
                  DATUM(2,IFRM)=9.9999
      WRITE (DATUM2(1,IFRM), *) (DATUM(1,IFRM))
      WRITE (DATUM2(2,IFRM), *) (DATUM(2,IFRM))                   
                  DATUM2(1,IFRM)="  INDEF"
                  DATUM2(2,IFRM)="  INDEF"
               END IF
            END DO
            SUMCHI=SUMCHI/NOBS(IMASTR)
            SUMSHP=SUMSHP/NOBS(IMASTR)
            IF ((XM(IMASTR) .LE. -999.9995) .OR. (YM(IMASTR) .LE.
     .           -999.9995) .OR. (XM(IMASTR) .GE. 9999.995) .OR.
     .           (YM(IMASTR) .GE. 9999.995)) THEN
               WRITE (1,111) IDMAST(IMASTR),
     .               XM(IMASTR), YM(IMASTR), 
     .               ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP
            ELSE
               WRITE (1,1116) IDMAST(IMASTR),
     .               XM(IMASTR), YM(IMASTR), 
     .               ((DATUM2(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP
            END IF
         END DO
         CALL CLFILE (1)
      END IF
C
C      CALL GETYN ('A file with the new transformations?', ANSWER)

      ANSWER = 'N'
      WRITE (*,*) 'A file with the new transformations?: ',ANSWER  

      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         OFILE = SWITCH(LFILE, CASE('.mch'))
         CALL GETFIL ('Output file name:', OFILE, 1, 'NEW')
         DO I=1,NFRM
            IF (MODE .EQ. 4) THEN
               CON(6,I)=FLIP(I)*CON(3,I)
               CON(5,I)=-FLIP(I)*CON(4,I)
            END IF
            IF (MODE .LE. 6) THEN
               WRITE (1,109) FILE(I), (CON(J,I), J=1,6), DMAG(I),
     .              SQRT(SMAG(I))
  109          FORMAT (1X, '''', A, '''', 2F10.3, 4F10.5, F10.3, F10.4,
     .              15F10.4)
            ELSE
               WRITE (1,109) FILE(I), (CON(J,I), J=1,6), DMAG(I),
     .              SQRT(SMAG(I)), (CON(J,I), J=7,MODE)
            END IF
         END DO
         CALL CLFILE (1)
      END IF
C
C      CALL GETYN ('A file with the transfer table?', ANSWER)

      ANSWER = 'N'
      WRITE (*,*) 'A file with the transfer table?: ',ANSWER 

      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         OFILE = SWITCH(LFILE, CASE('.tfr'))
         CALL GETFIL ('Output file name:', OFILE, 1, 'NEW')
         WRITE (1,112) (FILE(I), 99.9999, 9.9999, I=1,NFRM)
  112    FORMAT (1X, A30, 2F9.4)
         WRITE (1,113)
  113    FORMAT (1X, 30('='))
         DO IMASTR=1,NMASTR
            IF ((XM(IMASTR) .LE. -999.9995) .OR. (YM(IMASTR) .LE.
     .           -999.9995) .OR. (XM(IMASTR) .GE. 9999.995) .OR.
     .           (YM(IMASTR) .GE. 9999.995)) THEN
               WRITE (1,114) IDMAST(IMASTR),
     .              XM(IMASTR), YM(IMASTR), 
     .              (MAX0(0,IWHICH(IMASTR,IFRM)-LASTAR(IFRM-1)), 
     .              IFRM=1,NFRM)
  114          FORMAT (I6, 2F9.2, 18I6 / (24X, 18I6))
            ELSE
               WRITE (1,1145) IDMAST(IMASTR),
     .              XM(IMASTR), YM(IMASTR), 
     .              (MAX0(0,IWHICH(IMASTR,IFRM)-LASTAR(IFRM-1)), 
     .              IFRM=1,NFRM)
 1145       FORMAT (I6, 2F9.3, 18I6 / (24X, 18I6))
         END IF
         END DO
         CALL CLFILE (1)
      END IF
C
C      CALL GETYN ('Individual .COO files?', ANSWER)

      ANSWER = 'N'
      WRITE (*,*) 'Individual .COO files?: ',ANSWER 

      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         CALL TBLANK
         DO IFRM=1,NFRM
            CALL INFILE (1, FILE(IFRM), ISTAT)
            CALL RDCHAR (1, HEAD(1), K, ISTAT)
            CALL RDCHAR (1, HEAD(2), L, ISTAT)
            CALL CLFILE (1)
            OFILE = SWITCH(FILE(IFRM), CASE('.coo'))
            CALL OUTFIL (1, OFILE, ISTAT)
            WRITE (TEXT,6) 0, OFILE
    6       FORMAT (I6, 3X, A)
            CALL OVRWRT (TEXT(1:38), 4)
            IF ((K .GE. 3) .AND. (HEAD(1)(2:3) .EQ. 'NL')
     .           .AND. (OFILE .NE. 'APPEND')) THEN
               WRITE (1,115) HEAD(1)(1:K)
  115          FORMAT (A)
               IF (L .GE. 1) THEN
                  WRITE (HEAD(2)(1:3),115) '  1'
                  WRITE (1,115) HEAD(2)(1:L)
               ELSE
                  WRITE (1,*)
               END IF
               WRITE (1,*)
            END IF
C
            RCOL(IFRM) = RCOL(IFRM)+0.5
            RROW(IFRM) = RROW(IFRM)+0.5
            K = 0
            DO IMASTR=1,NMASTR
C
               CALL TRFM (NOC(1,IFRM), RCOL(1), RROW(1), FLIP(IFRM), 
     .              MODE, XM(IMASTR), YM(IMASTR), XX, YY)
C
               IF (HEAD(1)(2:3) .EQ. 'NL') THEN
                  IF ((XX .GE. 0.5) .AND. (XX .LE. RCOL(IFRM)) .AND. 
     .                 (YY .GE. 0.5) .AND. (YY .LE. RROW(IFRM))) THEN
                     K = K+1
                     IF (MOD(K,200) .EQ. 0) THEN
                        WRITE (TEXT,6) K, OFILE
                        CALL OVRWRT (TEXT(1:38), 2)
                     END IF
                     WRITE (1, 116) IDMAST(IMASTR), XX, 
     .                    YY, MM(IMASTR)-DMAG(IFRM)
                  END IF
               ELSE
                  K = K+1
                  IF (MOD(K,200) .EQ. 0) THEN
                     WRITE (TEXT,6) K, OFILE
                     CALL OVRWRT (TEXT(1:38), 2)
                  END IF
                  WRITE (1, 116) IDMAST(IMASTR), XX, 
     .                 YY, MM(IMASTR)-DMAG(IFRM)
  116             FORMAT (I6, 2F9.3, F9.4)
               END IF
            END DO
            WRITE (TEXT,6) K, OFILE
            CALL OVRWRT (TEXT(1:38), 2)
            CALL CLFILE (1)
         END DO
         CALL TBLANK
      END IF
C
C      CALL GETYN ('Simply transfer star IDs?', ANSWER)

      ANSWER = 'N'
      WRITE (*,*) 'Simply transfer star IDs?: ',ANSWER 

      IF (ANSWER .EQ. 'E') THEN
         CALL BYEBYE
      ELSE IF (ANSWER .EQ. 'Y') THEN
         DO IFRM=1,NFRM
            CALL INFILE (1, FILE(IFRM), ISTAT)
            CALL RDCHAR (1, HEAD(1), K, ISTAT)
            CALL RDCHAR (1, HEAD(2), L, ISTAT)
            IF (HEAD(1)(2:3) .NE. 'NL') THEN 
               CALL CLFILE (1)
               CALL INFILE (1, FILE(IFRM), ISTAT)
               NL = 1
            ELSE
               READ (1,*)
               READ (HEAD(2), *) NL
            END IF
C
            I=0
 8101       I=I+1
 8102       CALL RDCHAR (1, TEXT, III, ISTAT)
            IF (ISTAT .GT. 0) GO TO 8103
            IF (ISTAT .LT. 0) GO TO 8102
            READ (TEXT,*,ERR=8102) J, XX, YY, ZZ
            NLINE(I) = MIN0(127,III-6)
            IF ((J .LE. 0) .OR. (ZZ .GT. 90.)) GO TO 8102
            READ (TEXT,11,END=8103) III, (LINE(J,I), J=1,NLINE(I))
   11       FORMAT (I6, 127A1)
            IF (NL .EQ. 2) READ (1,*)
            GO TO 8101
C
 8103       CALL CLFILE (1)
            OFILE=SWITCH(FILE(IFRM), CASE('.mtr'))
            CALL OUTFIL (1, OFILE, ISTAT)
            WRITE (TEXT,6) 0, OFILE
            CALL OVRWRT (TEXT(1:38), 4)
            IF (HEAD(1)(2:3) .EQ. 'NL') THEN
               WRITE (1,115) HEAD(1)(1:K)
               WRITE (1,115) HEAD(2)(1:L)
               WRITE (1,*)
            END IF
C
            K = 0
            DO IMASTR=1,NMASTR
               J = IWHICH(IMASTR,IFRM)
               IF (J .GT. 0) THEN 
                  J = J - LASTAR(IFRM-1)
                  WRITE (1,117) IDMAST(IMASTR),
     .                 (LINE(L,J), L=1,NLINE(J))
  117             FORMAT (I6, 127A1)
                  K = K+1
               END IF
               IF (MOD(K,1000) .EQ. 0) THEN
                  WRITE (TEXT,6) K, OFILE
                  CALL OVRWRT (TEXT(1:38), 2)
               END IF
            END DO
            WRITE (TEXT,6) K, OFILE
            CALL OVRWRT (TEXT(1:38), 2)
            CALL CLFILE (1)
         END DO
         CALL TBLANK
      END IF
C
 9900 CONTINUE
      RETURN
      CALL BYEBYE
      END!
C
C##########################################################################
C
      INTEGER*4 FUNCTION  IBNRY  (X, N, F)
C
C Locate the first item in a sorted list whose value is greater
C than fraction F of the largest value in the list.
C
      REAL X(*)
      T = F*X(N)                                          ! Target value
C
      IF (X(1) .GT. T) THEN
         IBNRY = 1
         RETURN
      ELSE IF (X(N) .LT. T) THEN
         IBNRY = 0
         RETURN
      END IF
C
      K = MAX0(1, NINT(0.5*N))                             ! Initial step
      IBNRY = K
C
 1000 K = MAX0(1, (K+1)/2)
      IF (X(IBNRY) .LT. T) THEN
         IBNRY = MIN0(N, IBNRY+K)
         GO TO 1000
      END IF
      IF (X(IBNRY-1) .GT. T) THEN
         IBNRY = MAX0(2, IBNRY-K)
         GO TO 1000
      END IF
      RETURN
      END!
C
C##########################################################################
C
      SUBROUTINE  RDFRMS  (FILE, CHIMAX, NFRM, MAXMTR, WORK, MAXSTR, 
     .     IDMAST, LASTAR, RCOL, RROW)
      IMPLICIT NONE
C
      REAL CHIMAX, DM, RCOL(*), RROW(*)
      INTEGER J, K, ID, IFRM, ISTAR, IMAX, MAXMTR, MAXSTR, 
     .     NL, NFRM, NMAX, ISTAT, NCOL, NROW
C
      CHARACTER*30 FILE(NFRM)
      CHARACTER*140 LINE
      REAL WORK(MAXSTR,*), SOFT
      INTEGER IDMAST(MAXMTR), LASTAR(0:NFRM)
      DATA SOFT /1.E-6/
C
C Read data from all frames into the working array.
C
      ISTAR = 0
      NMAX = 0
      CALL OVRWRT (' ', 4)
      DO IFRM=1,NFRM
         IMAX = 0
         CALL INFILE (1, FILE(IFRM), ISTAT)
         IF (ISTAT .NE. 0) THEN
            CALL STUPID ('Error opening input file '//FILE(IFRM))
            CALL BYEBYE
         END IF
         NL = -1
c         WRITE (LINE,6) IFRM, FILE(IFRM)
         WRITE (*,6) IFRM, FILE(IFRM) ! AWS
    6    FORMAT ('  Reading file', I3,': ', A)
c         CALL OVRWRT (LINE(1:47), 2)
         CALL RDHEAD (1, NL, NCOL, NROW, DM, DM, DM, DM, DM, DM, DM)
         RCOL(IFRM) = REAL(NCOL)
         RROW(IFRM) = REAL(NROW)
 1010    ISTAR = ISTAR+1
 1020    CONTINUE
         IF (NL .EQ. 2) THEN
            READ (1,100,END=1090,ERR=1020) ID, (WORK(ISTAR,J), J=1,3)
C           READ (1,100,END=1090,ERR=1020) ID, (WORK(ISTAR,J), J=1,4)
 100        FORMAT (I6, 3F9.3)
            IF (ID .LE. 0) GO TO 1020
            READ (1,101,ERR=1020) WORK(ISTAR,4)
 101        FORMAT (25X, F8.4)
            WORK(ISTAR,5) = 0.
            WORK(ISTAR,6) = 0.
         ELSE
            CALL RDCHAR (1, LINE, K, ISTAT)
            IF (ISTAT .GE. 1) GO TO 1090
            IF (K .LT. 9) GO TO 1020
            LINE = LINE(1:K)//' 0 0 0 0 0'
            READ (LINE,*,ERR=1020) ID, (WORK(ISTAR,J), J=1,4),
     .        DM, DM, (WORK(ISTAR,J), J=5,6)
            IF (ID .LE. 0) GO TO 1020
            IF (WORK(ISTAR,5) .GT. CHIMAX) GO TO 1020
         END IF
         IF (WORK(ISTAR,3) .GT. 90.) GO TO 1020
         IMAX = IMAX+1
         if (file(ifrm) .eq. 'save:NGC2419.mag') then
            work(istar,4) = 4.
         else
            WORK(ISTAR,4) = WORK(ISTAR,4)**2 + SOFT
         end if
         IF (IFRM .EQ. 1) IDMAST(ISTAR)=ID
         IF (ISTAR .LT. MAXSTR) THEN
            GO TO 1010
         ELSE
            WRITE (6,690) FILE(IFRM)
  690       FORMAT (/' Star limit reached in ', A/)
            CALL BYEBYE
         END IF
 1090    ISTAR = ISTAR-1
         print*,'  (read up to index',istar,')' ! AWS
         LASTAR(IFRM) = ISTAR
         IF (LASTAR(IFRM) .LE. LASTAR(IFRM-1)) THEN
            CALL STUPID ('No stars in '//FILE(IFRM))
            CALL CLFILE (1)
            CALL OOPS
         END IF
         IF (IMAX .GT. NMAX) NMAX=IMAX
         CALL CLFILE (1)
      END DO
      MAXSTR = NMAX
      RETURN
      END!
C
C#########################################################################
C
      SUBROUTINE  REPACK  (WORK, MAXSTR, NTOT)
      IMPLICIT NONE
C
      INTEGER I, J, K, M, N, MAXSTR, NTOT
C
      REAL WORK(*)
C
      CALL OVRWRT ('Repacking memory                              ', 2)
      DO J=2,6
         M = NTOT*(J-1)
         N = MAXSTR*(J-1)
         DO I=1,NTOT
            K = M+I
            WORK(K) = WORK(N+I)
         END DO
      END DO
C
      RETURN
      END!
C
C########################################################################
C
      INTEGER FUNCTION  ICLOSE  (XX, YY, RADIUS, RADSQ, XM, YM, 
     .     NOBS, LAST, IFIRST, RMIN)
      PARAMETER (MINY=-11 384, MAXY=15 000)
      REAL XM(*), YM(*), RMIN(*)
      INTEGER LAST(MINY:MAXY), NOBS(*)
C
      ICLOSE=0
      L=MAX0(MINY, MIN0(MAXY, INT(YY-RADIUS)-1))
      IF (L .EQ. MINY) THEN
         L=IFIRST
      ELSE
         L=LAST(L)+1
      END IF
      M=LAST(MAX0(MINY, MIN0(MAXY, INT(YY+RADIUS))))
      RR=RADSQ
      DO IMASTR=L,M
        DX=XX-XM(IMASTR)
        IF (ABS(DX) .LE. RADIUS) THEN
          RSQ=DX**2+(YY-YM(IMASTR))**2
          IF (RSQ .LT. RADSQ) THEN
            RSQ = RSQ/REAL(NOBS(IMASTR)) 
            IF ((RSQ .LT. RR) .AND. (RSQ .LT. RMIN(IMASTR))) THEN
              RR=RSQ
              ICLOSE=IMASTR
            END IF
          END IF
        END IF
      END DO
      IF (ICLOSE .NE. 0) THEN
        RMIN(ICLOSE)=RR
      END IF
      RETURN
      END
C
C########################################################################
C
      SUBROUTINE  TRFM  (CON, RCOL, RROW, FLIP, MODE, X, Y, XX, YY)
      IMPLICIT NONE
      INTEGER MODE
      REAL CON(MODE), X, Y, XX, YY, X2, Y2, FLIP, XS, XY, YS, RCOL, RROW
      IF (MODE .EQ. 2) then
         XX=X+CON(1)
         YY=Y+CON(2)
      ELSE IF (MODE .EQ. 4) THEN
         XX=CON(1)+CON(3)*X-FLIP*CON(4)*Y
         YY=CON(2)+CON(4)*X+FLIP*CON(3)*Y
      ELSE IF (MODE .GE. 6) THEN
         XX=CON(1)+CON(3)*X+CON(5)*Y
         YY=CON(2)+CON(4)*X+CON(6)*Y
         IF (MODE .GE. 7) THEN
            XS=2.*(X-1.)/(RCOL-1.)-1.
            YS=2.*(Y-1.)/(RROW-1.)-1.
            XY=XS*YS
            X2=1.5*XS**2-0.5
            Y2=1.5*YS**2-0.5
            XX=XX+
     .           CON(7)*X2+
     .           CON(9)*XY+
     .           CON(11)*Y2
            YY=YY+
     .           CON(8)*X2+
     .           CON(10)*XY+
     .           CON(12)*Y2
            IF (MODE .GE. 13) THEN
               XX = XX+
     .              (CON(13)*XS+CON(15)*YS)*X2+
     .              (CON(17)*XS+CON(19)*YS)*Y2
               YY = YY+
     .              (CON(14)*XS+CON(16)*YS)*X2+
     .              (CON(18)*XS+CON(20)*YS)*Y2
            END IF
         END IF
      ELSE
         CALL STUPID ('TRFM MODE ERROR')
         CALL BYEBYE
      END IF
      RETURN
      END
C
      SUBROUTINE  INVERT  (CON, COL, ROW, COLM, ROWM, FLIP, MODE)
      PARAMETER (MAXCON=20)
      DOUBLE PRECISION C(MAXCON/2,MAXCON/2), VX(MAXCON/2), 
     .     VY(MAXCON/2), D(MAXCON/2)
      REAL CON(MAXCON), PLACES(4)
      DATA PLACES / 0.125, 0.375, 0.625, 0.875 /
      M = MODE/2
      IF (MODE .GT. 7) THEN
         DO I=1,MAXCON/2
            VX(I) = 0.0D0
            VY(I) = 0.0D0
            DO J=1,MAXCON/2
               C(J,I) = 0.0D0
            END DO
         END DO
         DO J=1,4
            Y = PLACES(J)*(ROW-1.)+1.
            DO I=1,4
               X = PLACES(I)*(COL-1.)+1.
               CALL TRFM(CON, COL, ROW, 
     .              FLIP, MODE, X, Y, XX, YY)
               D(1) = 1.
               D(2) = XX
               D(3) = YY
               XX = 2.*(XX-1.)/(COLM-1.)-1.
               YY = 2.*(YY-1.)/(ROWM-1.)-1.
               D(4) = 1.5*XX**2-0.5
               D(5) = XX*YY
               D(6) = 1.5*YY**2-0.5
               IF (MODE .GT. 13) THEN
                  D(7) = XX*D(4)
                  D(8) = YY*D(4)
                  D(9) = XX*D(6)
                  D(10) = YY*D(6)
               END IF
               DO L=1,M
                  VX(L) = VX(L)+D(L)*X
                  VY(L) = VY(L)+D(L)*Y
                  DO K=1,M
                     C(K,L) = C(K,L)+D(K)*D(L)
                  END DO
               END DO
            END DO
         END DO
         CALL DINVRS (C, MAXCON/2, M, IFLAG)
         CALL DVMUL (C, MAXCON/2, M, VX, D)
         DO K=1,MODE-1,2
            CON(K) = D((K+1)/2)
         END DO
         CALL DVMUL (C, MAXCON/2, M, VY, D)
         DO K=2,MODE,2
            CON(K) = (D(K/2))
         END DO
      ELSE 
         CON(1)=-CON(1)
         CON(2)=-CON(2)
         CONSTANT = CON(3)*CON(6) -  CON(4)*CON(5)
         CON(4) = -CON(4)/CONSTANT
         CON(5) = -CON(5)/CONSTANT
         C2 = CON(3)/CONSTANT
         CON(3) = CON(6)/CONSTANT
         CON(6) = C2
         CONSTANT = CON(3)*CON(1) + CON(5)*CON(2)
         CON(2) = CON(4)*CON(1) + CON(6)*CON(2)
         CON(1) = CONSTANT
         DO I=7,MAXCON
            CON(I)=-CON(I)
         END DO
      END IF
      RETURN
      END
