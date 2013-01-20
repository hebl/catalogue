      PROGRAM u4dump 
C
C  g77 -o u4dump u4dump.f u4sub.f
C
C  - loop over range of zones
C  - read individual UCAC4 zone file (binary data) 
C  - output all data to ASCII file, no header lines,
C      at least 1 blank between columns, no blank columns
C  - optional output of column separator character
C  - for explanations of columns see "readme" file
C
C  111227 NZ create
C  111228 NZ use separate subroutines similar to UCAC2 release, bf
C  120117 NZ change output file names to append ".asc"
C  120711 NZ test run, no algorithm change

      IMPLICIT NONE
      INTEGER    dimc, dims
      PARAMETER (dimc = 53   ! number of data columns
     .          ,dims = 52)  ! number of separators

      INTEGER      dv(dimc), mi(dimc), ma(dimc), is(dims)
      CHARACTER*1  csc

      INTEGER      i,j,k, jo, uni,uno,zn1,zn2,zn, nst,nsz
      CHARACTER*40 pathi,patho, answer
      CHARACTER*45 fnin,fnout, line*235, a1*1
      LOGICAL      eozf, bf

      DATA is /11,21,27,33,36,38,41,45,49,52,55,58,64,70,77,84,88,92
     .       ,103,109,115,121,124,127,130,133,136,139,145,151,157,163
     .       ,169,173,177,181,185,189,191,193,195,197,199,201,203,205
     .       ,207,209,212,215,225,229/

* defaults
      pathi = '/d3/temp/side1/u4b/'
      patho = './'   
      zn1 = 1
      zn2 = 1
      bf  = .FALSE.   ! byte flip of binary data

* interactive
      WRITE (*,'(/a)') 'dump UCAC4 binary data to ASCII'
      WRITE (*,'(a)') 
     .   'hit "enter" to accept defaults or enter new values'

      WRITE (*,'(2a,$)') 'path input = ',pathi
      READ (*,'(a)') answer
      IF (answer.NE.' ') pathi = answer

      WRITE (*,'(a)') 'binary data are stored in PC-style sequence'
      WRITE (*,'(a)') 'some computers (UNIX-style) need a byteflip'
      WRITE (*,'(a,l1,a,$)') 'byte flip (true/false) ?  ',bf,'  '
      READ (*,'(a)') a1
      IF (a1.NE.' ') READ (a1,*) bf

      WRITE (*,'(2a,$)') 'path output= ',patho
      READ (*,'(a)') answer
      IF (answer.NE.' ') patho = answer

      WRITE (*,'(a,2i4,a,$)') 'from, to zone (1..900) = ',zn1,zn2,'  '
      READ (*,'(a)') answer
      IF (answer.NE.' ') READ (answer,*) zn1,zn2

      WRITE (*,'(a,$)') 
     .  'enter column separator character (default = blank) '
      READ (*,'(a)') csc 

* prepare
      jo = INDEX (patho,' ') - 1
      uni= 11
      uno= 20
      nst= 0

* loop zone files
      DO zn = zn1,zn2
        WRITE (fnout,'(a,a,i3.3,a)') patho(1:jo),'z',zn,'.asc'
        CALL open_zfile (pathi,uni,zn,fnin)
        OPEN (uno,FILE=fnout)
        eozf = .FALSE.

        WRITE (*,'(/a,a)') 'begin read file = ',fnin
        WRITE (*,'( a,a)') '... output to   = ',fnout

        DO nsz = 1,999000
          CALL getistar (uni,nsz,bf,eozf,dv,dimc)
          IF (eozf) GOTO 91

          WRITE (line,'(2i10,2i6,i3,i2,i3,2i4,3i3,2i6,2i7,2i4
     .                ,i11,3i6,6i3,5i6,5i4,i2,9i2,2i3,i10,i4,i7)') 
     .      (dv(j),j=1,dimc)

          IF (csc.NE.' ') THEN
            DO j=1,dims
              k = is(j)
              line(k:k) = csc
            ENDDO
          ENDIF

          WRITE (uno,'(a)') line
        ENDDO   ! loop stars on individ. zone files

  91    CLOSE (uni)
        CLOSE (uno)
        nsz = nsz - 1
        nst = nst + nsz
        WRITE (*,'(a,i7,i10)') 'numb.stars/zone, total = ',nsz,nst
      ENDDO     ! loop all zones

      END  ! main <u4dump>
 
