      PROGRAM u4access
C
C  g77 -o u4access u4access.f u4sub.f 
C
C  - top level program to access original UCAC4 binary data zone files
C  - extract all items for stars in selected areas (RA, Dec box, mag range)
C  - formatted (ASCII) output table, various format options
C  - officieal UCAC4 star designation = zone numb + run ID along zone
C
C  use Fortran unit 13 to auxiliary data and zone files
C  use Fortran unit 20 for output file (or screen dump)
C
C  111228/29 NZ taken from UCAC2, lots of updates
C  120111 NZ taken from u4access.f, fmt = 5 = MPOS for RORF, list of fields
C  120131 NZ option enter degree for list of targets, split subr.
C  120321 NZ u4sub.f includes sorti.f;  hpmd array with parameter (UCAC4 vers.4)
C  120711 NZ allow degree input also with interactive mode

      IMPLICIT NONE
      INTEGER    zmax,bmax, dimj, dima, dimh
      PARAMETER (zmax= 900       ! largest zone number
     .          ,bmax=1440       ! numb. bins along RA for index
     .          ,dimj= 59        ! all data items + upd. pos,err + ID numbers
     .          ,dima=  90000    ! max.numb. stars in a box
     .          ,dimh= 99)       ! max.numb. of HPM stars from external file

      INTEGER  alls(dima,dimj)   ! all data all selected star
      INTEGER  x0 (zmax,bmax)    ! index = accumul.numb.stars f(zone,RA bin)
      INTEGER  hpmd(dimh,5)      ! HPM data from extra table
      INTEGER  fmto

      INTEGER  i,j,jp,is,zn, nsr,nss,errc, uo, rah,ram,dcd,dcm
     .        ,un,ic1, nst, ntgt, nh
      REAL*8   ra1,ra2,dc1,dc2,ra0,dc0,ras,dcs,mag1,mag2, newep,wra,wdc
      CHARACTER*40 pathz,fnxu,fnhpm, fnout, answer, fnlst
      CHARACTER*1  a1, lstfmt
      LOGICAL      bf, was, negdec

* defaults
      pathz = '/d3/temp/side2/u4b/'  ! path to relase files on HDD
CC    pathz = '/mnt/dvd/u4z/'    ! CDROM at my Linux box
CC    pathz = 'D:\u4\'           ! sample path for Windows, CD-ROM drive

      fnxu  = '/d3/temp/side1/u4i/u4index.unf'
      fnhpm = '/d3/temp/side1/u4i/u4hpm.dat'
      fnout = 'u4.tab'    

      ra1   =  0.0d0
      ra2   = 24.0d0
      dc1   = -90.0d0
      dc2   = -89.9d0
      mag1  = -2.0d0
      mag2  = 21.0d0
      newep = 2000.0d0
      wra   = 1.0d0
      wdc   = 0.2d0
      fmto  = 1
      was   = .FALSE.    ! width in arcsec, else degree
      bf    = .FALSE.    ! byte flip needed
      un    = 13         ! Fortran unit number for zone files
      is    = 0          ! sort by column number (0 = no sort)
      ic1   = 1          ! first column to sort by
      nst   = 0          ! sum up total numb. stars output
      ntgt  = 0          ! number of targets from list file

* interactive: basic data, common for entire run of program
      WRITE (*,'(/a)') 'Program to access UCAC4 release data'
      WRITE (*,'( a)') '  hit "enter" for defaults'
      WRITE (*,'( a)') 'option: loop through several areas possible,' 
      WRITE (*,'( a)') '  however, output with single format to single'
      WRITE (*,'( a)') '  file per run of this program only'
      WRITE (*,'( a)') 'wrap around 24/0h RA allowed'

      WRITE (*,'(/a)') 'path for binary UCAC4 zone files = '
      WRITE (*,'(a,$)') pathz
      READ (*,'(a)') answer
      IF (answer.NE.' ') pathz = answer

      WRITE (*,'(/a)') 'binary data are stored in PC-style sequence'
      WRITE (*,'( a)') 'some computers (UNIX-style) need a byteflip'
      WRITE (*,'(a,l1,a,$)') 'byte flip (true/false) ?  ',bf,'  '
      READ (*,'(a)') a1
      IF (a1.NE.' ') READ (a1,*) bf

      WRITE (*,'(/a)') 'full path and file name of unform. index ='
      WRITE (*,'(a,$)') fnxu
      READ (*,'(a)') answer
      IF (answer.NE.' ') fnxu = answer

      WRITE (*,'(/a)') 'full path and file name for high PM stars='
      WRITE (*,'(a,$)') fnhpm
      READ (*,'(a)') answer
      IF (answer.NE.' ') fnhpm = answer

      WRITE (*,'(/a)') '  fmt = 1 = all original integers'
      WRITE (*,'( a)') '  fmt = 2 = updated RA, Dec (degree) PM..'
      WRITE (*,'( a)') '  fmt = 3 = updated RA, Dec (hms) photom..'
      WRITE (*,'( a)') '  fmt = 4 = updated RA, Dec (h,d) short'
      WRITE (*,' (a)') '  fmt = 5 = MPOS for RORF ref.stars'
      WRITE (*,' (a)') '  fmt = ... user defined in subr. put_stars'
      WRITE (*,'(a,i2,a,$)') 'format for output = ',fmto,'  '
      READ (*,'(a)') answer
      IF (answer.NE.' ') READ (answer,*) fmto

      WRITE (*,'(/a)') 'file name output table or "s" = screen dump '
      WRITE (*,'(a,$)') fnout
      READ (*,'(a)') answer
      IF (answer.NE.' ') fnout = answer

      WRITE (*,'(/a)')
     . 'select apert. magnitude range (UCAC mags, -2..21 = all):'
      WRITE (*,'(2f7.2,a,$)') mag1, mag2,'  '
      READ (*,'(a)') answer
      IF (answer.NE.' ') READ (answer,*) mag1,mag2

      IF (fmto.GT.1) THEN
        WRITE (*,'(/a)')
     .   'update positions for proper motions to new epoch'
        WRITE (*,'( a)') '2000.0 is default = original UCAC4'
        WRITE (*,'(f7.2,a,$)') newep,'  '
        READ (*,'(a)') answer
        IF (answer.NE.' ') READ (answer,*) newep
      ENDIF

      WRITE (*,'(/a)') 'see readme4.txt for all items'
      WRITE (*,'( a)') '  0 = no sort = default'
      WRITE (*,'( a)') '  1 = sort by col. 1 = RA '
      WRITE (*,'( a)') '  2 = sort by col. 2 = Dec'
      WRITE (*,'( a)') '  4 = sort by UCAC apert.mag'
      WRITE (*,'( a)') ' ...'
      WRITE (*,'( a)') ' 53 = sort by UCAC4 ID number'
      WRITE (*,'(a,i2,a,$)') 'sort by column number: ',is,'  '
      READ (*,'(a)') answer
      IF (answer.NE.' ') READ (answer,*) is

* read HPM star data
      CALL gethpm (fnhpm,hpmd,dimh,nh)

* read index file (900 x 1440 x 4 byte)
      WRITE (*,'(/a,a)') 'open index file = ',fnxu
      OPEN (un,FILE=fnxu,ACCESS='direct',RECL=5184000)
      READ (un,REC=1) x0  

C     OPEN (un,FILE=fnxu,ACCESS='direct',RECL=5760) ! 1440 * 4 byte
C     DO zn=1,zmax
C       READ (un,REC=zn,ERR=901) (x0(zn,j),j=1,bmax)
C     ENDDO

      CLOSE(un)
      WRITE (*,'(a)') 'index read successfully'
      DO zn=1,4
        WRITE (*,'(a,6i6)') 'zn, x0...=',zn,(x0(zn,j),j=1,5)
      ENDDO

      IF (bf) THEN
        CALL x0_byte_flip (x0,zmax,bmax) ! apply byte flip if needed
      ENDIF

* prepare table output
      IF (fnout.EQ.'s') THEN
        uo = 6
      ELSE
        uo = 20
        OPEN (uo,FILE=fnout)
      ENDIF

* input list or interactive
      WRITE (*,'(/a,$)') 'list of centers or blank for interactive: '
      READ (*,'(a)') fnlst

      WRITE (*,'(a,$)') '(h) = hms format,  (d) = decimal degree  '
      READ (*,'(a)') lstfmt

      IF (fnlst.NE.' ') THEN
        OPEN (12,FILE=fnlst)
        WRITE (*,'(/a,2f7.3,a,$)')
     .     'width of box RA*cosD, Dec [degree] =',wra,wdc,'  '
        READ (*,'(a)') answer
        IF (answer.NE.' ') READ (answer,*) wra,wdc
      ENDIF

* loop boxes 
 101  CONTINUE

      IF (fnlst.EQ.' ') THEN
        WRITE (*,'(/a,$)') '== new box: r=range, c=center, q=quit '
        READ (*,'(a)') a1

        IF (a1.EQ.'r'.OR.a1.EQ.'R') THEN
          WRITE (*,'(a,2f8.4,a,$)') 'RA (hour) = ',ra1,ra2,'  '
          READ (*,'(a)') answer
          IF (answer.NE.' ') READ (answer,*) ra1,ra2
          WRITE (*,'(a,2f8.3,a,$)') 'Dec (deg) = ',dc1,dc2,'  '
          READ (*,'(a)') answer
          IF (answer.NE.' ') READ (answer,*) dc1,dc2

        ELSEIF (a1.EQ.'c'.OR.a1.EQ.'C') THEN
          WRITE (*,'(2a)') 'width in degree or enter "a" for arcsec'
     .      ,' after values'
          WRITE (*,'(a,2f7.3,a,$)')
     .       'width of box RA*cosD, Dec =',wra,wdc,'  '
          READ (*,'(a)') answer
          j = INDEX(answer,'a')
          IF (j.GT.0) THEN
            was = .TRUE.
            answer (j:j) = ' '
          ELSE
            was = .FALSE.
          ENDIF

          IF (answer.NE.' ') THEN
            READ (answer,*) wra,wdc
            IF (was) THEN
              wra = wra / 3.6d3    ! arcsec to degree
              wdc = wdc / 3.6d3
            ENDIF
          ENDIF
          IF (wra.LE.0.0d0.OR.wdc.LE.0.0d0) GOTO 101

          IF (lstfmt.EQ.'h') THEN
            WRITE (*,'(a,$)') 'center RA (hms), Dec (dms) = '
            READ (*,'(a)') answer
            CALL str2deg  (answer,ra0,dc0)  ! [degree] both
            CALL calc_box (ra0,dc0,wra,wdc,ra1,ra2,dc1,dc2)
            WRITE (*,'(2f10.6,2f10.5)') ra1,ra2,dc1,dc2
          ELSE
            WRITE (*,'(a,$)') 'center RA (deg), Dec (deg) = '
            READ (*,*) ra0,dc0     ! degree
            WRITE (*,'(a,4f12.6)') 'ra0,dc0 = ',ra0,dc0, wra,wdc
            CALL calc_box (ra0,dc0,wra,wdc,ra1,ra2,dc1,dc2)
          ENDIF

        ELSE
          GOTO 99  ! exit interactive loop
        ENDIF
 
      ELSE         ! read centers from list file
        IF (lstfmt.EQ.'h') THEN
          READ (12,'(a)',END=99) answer  ! hms format
          CALL str2deg  (answer,ra0,dc0)
          CALL calc_box (ra0,dc0,wra,wdc,ra1,ra2,dc1,dc2)
        ELSE
          READ (12,*,END=99) ra0,dc0     ! degree
          WRITE (*,'(a,4f12.6)') 'ra0,dc0 = ',ra0,dc0, wra,wdc
          CALL calc_box (ra0,dc0,wra,wdc,ra1,ra2,dc1,dc2)
        ENDIF
        ntgt = ntgt + 1
        WRITE (*,'(/a,2f10.6,2f10.5)') 'new field:', ra1,ra2,dc1,dc2
      ENDIF

* output info
      IF (fmto.NE.5.AND.fnlst.EQ.' ') THEN
        WRITE (uo,'(/a, f10.3)') 'epoch      = ',newep
        WRITE (uo,'( a, i10  )') 'sort item  = ',is
        WRITE (uo,'( a, i10  )') 'output fmt = ',fmto
        WRITE (uo,'( a,2f10.2)') 'mag range  = ',mag1,mag2
        WRITE (uo,'( a,f10.6,f10.5)') 'min RA, Dec= ',ra1,dc1
        WRITE (uo,'( a,f10.6,f10.5)') 'max RA, Dec= ',ra2,dc2
      ENDIF

* get stars 
      CALL get_stars (bf,ra1,ra2,dc1,dc2,mag1,mag2, hpmd,dimh,nh
     .               ,newep,pathz,x0,zmax,bmax,dima,dimj
     .               ,alls,nss,nsr,errc)

      WRITE (*,'(a,3i7)')
     .   'stars read, selected, eofile = ',nsr,nss,errc

* sort
      IF (is.GE.1.AND.is.LE.dimj) THEN
        WRITE (*,'(a,i6)') 'start sort by column = ',is
C       WRITE (*,'(5i7)') is,ic1,nss,dima,dimj
        CALL sorti (alls,is,dimj,ic1,nss,dima,dimj)
      ENDIF

* output
      CALL put_stars (ntgt,newep,alls,dima,dimj,nss,fmto,uo)
      WRITE (*,'(a,i6,a)') 'output complete for ',nss,' stars'
      nst = nst + nss

      GOTO 101         ! more boxes

* the end
 901  WRITE (*,'(a)') 'error in reading index file'
      STOP

  99  CONTINUE
      WRITE (*,'(/a,i8)') 'total numb. stars output = ',nst
      WRITE (*,'(a/)') 'thanks for using <u4access>, have a nice day!'

      END  ! main <u4access>

************************************************************************

      SUBROUTINE str2deg (answer,ra0,dc0)
C
C  input : answer = string with RA, Dec in hms format
C  output: ra0,dc0= RA, Dec [degree]
C
C  120323 NZ fix ra0 [degree]

      IMPLICIT NONE
      CHARACTER*(*) answer
      REAL*8        ra0,dc0, ras,dcs
      INTEGER       rah,ram,dcd,dcm, j
      LOGICAL       negdec

      j = INDEX (answer,'-')

      IF (j.GT.0) THEN
        negdec = .TRUE.
        answer (j:j) = ' '
      ELSE
        negdec = .FALSE.
      ENDIF
      READ (answer,*) rah,ram,ras, dcd,dcm,dcs

      ra0 =(rah + ram/60.0d0 + ras/3.6d3) * 15.0d0 ! hour,min,sec -> degree
      dc0 = dcd + dcm/60.0d0 + dcs/3.6d3           !      -> decimal degree
      IF (negdec) dc0 = -dc0

      IF (ra0.LT.  0.0d0.OR.ra0.GT.360.0d0.OR.
     .    dc0.LT.-90.0d0.OR.dc0.GT. 90.0d0) THEN
        WRITE (*,'(/a)') 'error <str2deg>'
        WRITE (*,'(a)') answer
        WRITE (*,'(2f12.6)') ra0,dc0
        STOP
      ENDIF

      END  ! subr. <str2deg>

************************************************************************

      SUBROUTINE calc_box (rad,dcd,wra,wdc,ra1,ra2,dc1,dc2)
C
C  input : rad,dcd= RA, Dec [degree]
C          wra    = width of box [degree] along RA (*cosDec)
C          wdc    = width of box [degree] along Dec
C  output: ra1,ra2= range in RA [hour]
C          dc1,dc2= range in Declination [degree]

      IMPLICIT NONE
      REAL*8  rad,dcd, wra,wdc, ra1,ra2,dc1,dc2, rah
      REAL*8  degrad, cosdec, hwr,hwd, ras,dcs

      degrad = DATAN (1.0d0) / (4.5d1) ! degree to radian
      cosdec = DCOS (dcd * degrad)
      rah = rad / 15.0d0               ! degree -> hour

      hwr = wra / (cosdec * 30.0d0)    ! half width in hour for RA
      hwd = wdc / 2.0d0                ! half width in Dec

      ra1 = rah - hwr
      ra2 = rah + hwr
      IF (ra1.LT. 0.0d0) ra1 = ra1 + 24.0d0
      IF (ra2.GT.24.0d0) ra2 = ra2 - 24.0d0

      dc1 = dcd - hwd
      dc2 = dcd + hwd

      END  ! subr. <calc_box>

************************************************************************

      SUBROUTINE get_stars (bf,ra1,ra2,dc1,dc2,mag1,mag2, hpmd,dimh,nh
     .                     ,newep,pathz,x0,zmax,bmax,dima,dimj
     .                     ,alls,nss,nsr,errc)
C
C  get all data items for all stars within specified RA,Dec box and mag range
C
C  input : bf       = .TRUE. if flip of bytes is required
C          ra1,ra2  = RA range (hour, decimal)
C          dc1,dc2  = declination range (degree, decimal)
C          mag1,mag2= range aperture magnitude 
C          hpmd(dimh,5) = supplement data for high proper motion stars
C          dimh,nh  = dimension, actual number of HPM stars to use
C          newep    = requested epoch for positions (year)
C          pathz    = path name for UCAC4 zone files 
C          x0 (zmax,bmax) = index array = f(zone, RA bin)
C            = run along indiv. zone only = 1st star of bin - 1
C          zmax     = largest zone number  (0.2 deg steps from South Pole)
C          bmax     = largest bin along RA (0.25 deg steps = 1 arcmin)
C          dima     = max. number of stars to retrieve
C          dimj     = number of items per star incl. updated pos.+error
C
C  output: alls(dima,dimj) = all items for all retrieved stars
C          nss      = number of stars selected = within specified box
C          nsr      = number of star records read from zone files
C          errc     = number of errors while reading zone files

      IMPLICIT NONE
      INTEGER    dimc
      PARAMETER (dimc = 53)  ! numb. columns (data items from file)

      LOGICAL  bf
      REAL*8   ra1,ra2,dc1,dc2, mag1,mag2, newep
      CHARACTER*(*) pathz
      INTEGER  zmax,bmax,dima,dimj, nb, nsr,nss,errc, dimh,nh
      INTEGER  x0 (zmax,bmax)    ! accumul.numb.stars = f(zone,RA bin)
      INTEGER  hpmd(dimh,5)      ! supplem. data for very high PM stars
      INTEGER  alls(dima,dimj)   ! all integer data each star

      INTEGER  idat(dimc)        ! data items individual star

      INTEGER  i,j,zn,rr,uz, d1m,d2m,z1,z2,nz, mi1,mi2
     .        ,r1m(2),r2m(2),i1(2),i2(2),nr, j1,j2, recn
     .        ,ran,dcn, sxn,syn, sd1,sd2

      REAL*8        masrad
      CHARACTER*40  fnzone
      LOGICAL  errflg

* prepare
      IF (dimj.LT.59) THEN
        WRITE (*,'(a)') 'error <get_stars> need dimj >= 59'
        STOP
      ENDIF

      mi1 = IDNINT (mag1 * 1.0d3) ! from R*8 magnitudes
      mi2 = IDNINT (mag2 * 1.0d3) !   to I*4 [milli-mag]
      masrad = DATAN (1.0d0) / (4.5d1*3.6d6)  ! mas to radian

* calculate range 
      CALL get_zone_range (dc1,dc2,zmax,d1m,d2m,z1,z2,nz)

      WRITE (*,'(a,2i5,2i11)') '=== z1,2,d1,2 = ',z1,z2,d1m,d2m

      sd1 = d1m + 324000000       ! Dec -> SPD [mas]
      sd2 = d2m + 324000000

      IF (nz.LT.1) THEN
        nss = 0
        RETURN    ! exit, declination invalid
      ENDIF

      CALL get_ra_range (ra1,ra2,r1m,r2m,i1,i2,nr)

      WRITE (*,'(a,4i11)') '=== r1m,2 = ',r1m,r2m
      WRITE (*,'(a,4i11)') '=== i1,i2 = ',i1,i2

* initialize
      uz = 13     ! unit number for zone files
      nsr = 0
      nss = 0
      errc= 0

* loop all zones
      DO zn = z1,z2
        WRITE (*,'(a,i3)') 'open file for zone = ',zn
        CALL open_zfile (pathz,uz,zn,fnzone)  
        WRITE (*,'(a)') fnzone

* rr    = 1 or 2     = numb. of RA ranges, normal=1, else crossover 24/0 hour
* i1,i2 = 1,...,1440 = bin number (0.25 deg) along RA
* j1,j2 = 1,... many = record number on individ. zone file (1 record = 1 star)
* x0(..)= j - 1  for that zone and RA bin

        DO rr = 1,nr               ! 1 or 2 RA ranges possible
          errflg = .FALSE.         ! reset, no end of file encountered yet
          j1 = x0 (zn,i1(rr)) + 1  ! record numb. of 1st star along zone file
          IF (i2(rr).GE.1440) THEN ! last RA bin included for output
            j2 = 999000            ! read till end of zone file
          ELSE
            j2 = x0 (zn,i2(rr)+1)  ! index of next bin = last star of prev.
          ENDIF

          WRITE (*,'(a,i3,3i5,2i9)') 
     .     'zn, rr,i1,i2, j1,j2 = ',zn,rr,i1(rr),i2(rr),j1,j2

          DO recn= j1,j2
            CALL getistar (uz,recn,bf,errflg,idat,dimc)

            IF (errflg) THEN
              errc = errc + 1  ! count read errors = end of files encountered
              GOTO 91
            ELSE
              nsr = nsr + 1    ! count number of stars read
              IF (idat(4).GE.mi1.AND.idat(4).LE.mi2) THEN   ! apert.mag range
              IF (idat(1).GE.r1m(rr).AND.idat(1).LE.r2m(rr).AND.
     .            idat(2).GE.sd1.AND.idat(2).LE.sd2) THEN
                nss = nss + 1
                IF (nss.GT.dima) THEN
                  WRITE (*,'(/a,i8)')
     .             'WARNING: too many stars, take first dima =',dima
                  RETURN
                ENDIF

                IF (idat(15).EQ.32767) THEN         ! look up real PM for
                  DO j=1,nh                         !   those very HPM stars
                    IF (idat(51).EQ.hpmd(j,1)) THEN ! unique star ID number
                      idat(15) = hpmd(j,4)          ! PM RA*cosDec
                      idat(16) = hpmd(j,5)          ! PM Dec
                      GOTO 300
                    ENDIF
                  ENDDO   ! loop all supplem. data very high PM stars
                  WRITE (*,'(a,i9)') 'error HPM not found =',idat(51)
                  STOP
 300              CONTINUE
                ENDIF

                DO j=1,dimc   ! all data columns
                  alls (nss,j) = idat(j) 
                ENDDO

                CALL apply_pm (idat,newep,masrad, ran,dcn,sxn,syn)
                alls (nss,54) = ran
                alls (nss,55) = dcn
                alls (nss,56) = sxn
                alls (nss,57) = syn
                alls (nss,58) = zn      ! zone number
                alls (nss,59) = recn    ! running record numb. along zone

              ENDIF    ! within box
              ENDIF    ! within mag range
            ENDIF      ! case read o.k. or error

          ENDDO        ! all stars within a zone
  91      CONTINUE     ! end of zone file
          WRITE (*,'(a,i2,2i6)') '=== rr,nsr,nss = ',rr,nsr,nss
        ENDDO          ! 1 or 2 RA ranges 
      ENDDO            ! all zones

      END   ! subr. <get_stars>

************************************************************************

      SUBROUTINE put_stars (ntgt,newep,alls,dima,dimj,nss,fmto,uo)
C
C  input : ntgt  = number of target (from list) or 0
C          newep = epoch of position data 
C          alls  = array with selected stars
C          dima  = max.numb. of stars
C          dimj  = number of items per star
C          nss   = number of selected stars
C          fmto  = format modus for output (= 1 to 4 here)
C          uo    = Fortran unit number for output file
C  output : items to file with unit number uo (can be screen = std.out)
C           header for ntgt <= 1

      IMPLICIT NONE
      INTEGER  dima,dimj,nss, fmto,uo
      INTEGER  alls(dima,dimj)

      INTEGER  ntgt,i,j, nrj,spd,ep97, nmag,area,wdsf, dnn,hstf,badf
      REAL*8   newep,magm,maga, epr,epd, pmr,pmrc,pmd
     .        ,ran,dcn, magj, magh, magk, epmr,epmd
     .        ,magb,magv,magg,magr,magi
      CHARACTER*13 cra,cdc

* header line
      IF (ntgt.LE.1) THEN
      IF (fmto.EQ.1) THEN
        WRITE (uo,'(/2i10,2i6,i3,i2,i3,2i4,3i3,2i6,2i7,2i4
     .            ,i11,3i6,6i3,5i6,5i4,i3,9i2,2i3,i10,i4,i7)')
     .      (j,j=1,53)
        WRITE (uo,'(7a)') '        RA       SPD'
     .   ,' momag apmag sm o ds sra sdc nt ns nc'
     .   ,' cepra cepdc  pmrac   pmdc spr spd    2MASSid'
     .   ,'  Jmag  Hmag  Kmag Ji Hi Ki Je He Ke'
     .   ,'  Bmag  Vmag  gmag  rmag  imag Ber Ver ger rer ier'
     .   ,' gc T A B H Z U L N S le 2x  uniqueID U2z  U2rsn'
     .   ,' U4z  U4rsn'

      ELSEIF (fmto.EQ.2) THEN
        WRITE (uo,'(/3a)') '   RA (deg)    DE (deg)'
     .   ,' momag apmag eRA eDE    epRA    epDC'
     .   ,' pmRAcD   pmDE epmR epmD  uniqueID U4z  U4rsn'

      ELSEIF (fmto.EQ.3) THEN
        WRITE (uo,'(/3a)') '           RA            DE'
     .   ,' eRA eDE na nu nc  momag  mamag   Jmag   Hmag   Kmag'
     .   ,'   Bmag   Vmag   gmag   rmag   imag  uniqueID U4z  U4rsn'

      ELSEIF (fmto.EQ.4) THEN
        WRITE (uo,'(/2a)') '  RA (hour)    DE (deg)'
     .   ,' apmag eRA eDE U4z  U4rsn'

      ELSEIF (fmto.EQ.5) THEN
        WRITE (*,'(a,i6)') 'output nss = ',nss 

      ELSE
        WRITE (*,'(a,i3)') '** invalid format  fmto = ',fmto
        RETURN
      ENDIF
      ENDIF  ! header only for 1st target from list or interactive cases

* table with data
      DO i=1,nss
        IF (fmto.EQ.1) THEN               ! all original integer 
          WRITE (uo,'(2i10,2i6,i3,i2,i3,2i4,3i3,2i6,2i7,2i4,i11
     .              ,3i6,6i3,5i6,5i4,i3,9i2,2i3,i10,i4,i7,i4,i7)')
     .      (alls(i,j),j=1,53), alls(i,58),alls(i,59)

        ELSEIF (fmto.EQ.2) THEN
          ran = alls(i,54) / 3.6d6       ! updated mas -> degree
          dcn = alls(i,55) / 3.6d6       ! updated mas -> degree
          magm= alls(i, 3) * 1.0d-3
          maga= alls(i, 4) * 1.0d-3
          epr = alls(i,13) * 1.0d-2 + 1900.0d0
          epd = alls(i,14) * 1.0d-2 + 1900.0d0
          pmrc= alls(i,15) * 0.1d0
          pmd = alls(i,16) * 0.1d0
          epmr= alls(i,17) * 0.1d0
          epmd= alls(i,18) * 0.1d0

          WRITE (uo,'(f11.7,f12.7,2f6.2,2i4,2f8.2
     .               ,2f7.1,2f5.1,i10,i4,i7)')
     .      ran,dcn,magm,maga, alls(i,56),alls(i,57), epr,epd
     .     ,pmrc,pmd, epmr,epmd, alls(i,51),alls(i,58),alls(i,59)

        ELSEIF (fmto.EQ.3) THEN
          ran = alls(i,54) * 1.0d-3   ! updated pos. mas -> arcsec 
          dcn = alls(i,55) * 1.0d-3   ! updated pos. mas -> arcsec
          magm= alls(i, 3) * 1.0d-3
          maga= alls(i, 4) * 1.0d-3
          magj= alls(i,20) * 1.0d-3   ! from 2MASS
          magh= alls(i,21) * 1.0d-3
          magk= alls(i,22) * 1.0d-3
          magb= alls(i,29) * 1.0d-3   ! APASS mags
          magv= alls(i,30) * 1.0d-3
          magg= alls(i,31) * 1.0d-3
          magr= alls(i,32) * 1.0d-3
          magi= alls(i,33) * 1.0d-3

          CALL as2hms (ran,dcn,cra,cdc)

          WRITE (uo,'(a,1x,a,2i4,3i3,10f7.3,i10,i4,i7)')
     .    cra,cdc,alls(i,56),alls(i,57),alls(i,10),alls(i,11),alls(i,12)
     .   ,magm,maga,magj,magh,magk,magb,magv,magg,magr,magi, alls(i,51)
     .   ,alls(i,58),alls(i,59)

        ELSEIF (fmto.EQ.4) THEN
          ran = alls(i,54) / 5.4d7       ! updated mas -> hour
          dcn = alls(i,55) / 3.6d6       ! updated mas -> degree
          maga= alls(i, 4) * 1.0d-3
          WRITE (uo,'(f11.8,f12.7,f6.2,2i4,i4,i7)')
     .      ran,dcn,maga, alls(i,56),alls(i,57), alls(i,58),alls(i,59)

        ELSEIF (fmto.EQ.5) THEN
          nrj = 0
          spd = alls(i,55) + 324000000
          ep97= IDNINT ((newep-1997.0)*1.0d3)
          nmag= 1
          area= 9
          wdsf= 0
          dnn = 100
          hstf= 0
          badf= 0

          WRITE (uo,'(i10,1x,i9,4(1x,i4),4(1x,i3),1x,i5,2(1x,i5)
     .               ,2(1x,i4),2(1x,i2),1x,i4,2(1x,i2),1x,i5
     .               ,2(1x,i2),1x,i5)')
     .        alls(i,54),spd, alls(i,8),alls(i,9)
     .       ,alls(i,56),alls(i,57), alls(i,10),alls(i,11), nrj
     .       ,alls(i, 6),ep97, alls(i,3), alls(i,4),alls(i,5),alls(i,5)
     .       ,nmag,nmag,area, alls(i,7),wdsf,dnn,hstf,badf
     .       ,alls(i,3)-alls(i,4)

        ENDIF ! case fmto

      ENDDO   ! all stars

      END   ! subr. <put_stars>

************************************************************************

      SUBROUTINE pos_up (ra,spd,pmra,pmdc,newep, ran,dcn)
C
C  update position for proper motion with epoch difference
C
C  input : ra, spd   = RA, SPD (mas) at J2000 epoch assumed
C          pmra,pmdc = proper motion RA, Dec (0.1 mas/yr)
C                      no cos(Dec) for pmra, i.e. large values near pole
C          newep     = new epoch (year)
C  output: ran,dcn   = RA, Dec (mas) at new epoch

      IMPLICIT NONE
      INTEGER  ra,spd,pmra,pmdc,ran,dcn
      REAL*8   newep, dt, dra, ddc

      dt  = newep - 2000.0d0
      dra = DFLOAT(pmra) * 0.1d0 * dt
      ddc = DFLOAT(pmdc) * 0.1d0 * dt

      ran = ra + IDNINT (dra)
      IF (ran.GT.1296000000) ran = ran - 1296000000  ! 24 hour in mas
      IF (ran.LT.         0) ran = ran + 1296000000  ! restrict to 0..24 range

      dcn = spd - 324000000 + IDNINT (ddc) ! assume no jump over pole

CC    WRITE (20,'(a,f9.3,2f9.1)') 'pos_up: dt,dra,ddc = ',dt,dra,ddc

      END    ! subr. <pos_up>

************************************************************************

      SUBROUTINE pos_error (sigx,sigy,cepx,cepy,spmx,spmy,newep
     .                     ,sxn,syn)
C
C calculate position error at new epoch
C
C     x = RA * cosDec,  y = Dec  component
C
C input : sigx, sigy = position error (mas) at central epoch
C         cepx, cepy = central epoch (1/100 yr after 1900)
C         spmx, spmy = error of proper motion component (0.1 mas/yr)
C         newep      = requested epoch for position error (year)
C output: sxn, syn   = position error at new epoch (mas)

      IMPLICIT NONE
      INTEGER  sigx,sigy, spmx,spmy, cepx,cepy, sxn,syn
      REAL*8   newep, dtx,dty

      dtx = newep - DFLOAT(cepx) * 1.0d-2 - 1900.0d0
      dty = newep - DFLOAT(cepy) * 1.0d-2 - 1900.0d0

      sxn = IDNINT (DSQRT (sigx**2 + (spmx * 0.1d0 * dtx)**2))
      syn = IDNINT (DSQRT (sigy**2 + (spmy * 0.1d0 * dty)**2))

CC    WRITE (20,'(a,2f9.3)') 'pos_error: dtx,dty = ',dtx,dty

      END   ! subr. <pos_error>

************************************************************************

      SUBROUTINE apply_pm (idat,newep,masrad, ran,dcn,sxn,syn)
C
C input : idat = integer vector with original UCAC2 data for 1 star
C         newep= requested new epoch (year)
C         masrad = convert mas to radian
C output: ran, dcn = new position, updated for proper motion (mas)
C         sxn, syn = error of new position (at new epoch) (mas)

      IMPLICIT NONE
      INTEGER  idat(25), ran,dcn,sxn,syn, pmra
      REAL*8   newep, masrad, dec,cosd

      dec  = DFLOAT(idat(2)-324000000)     ! SPD [mas] -> Dec [mas]
      cosd = DCOS (dec*masrad)
      pmra = IDNINT (idat(15)/cosd)        ! PM RA (without cosDec factor)

      CALL pos_up (idat(1),idat(2),pmra,idat(16),newep, ran,dcn)

      CALL pos_error (idat( 8),idat( 9), idat(13),idat(14)
     .               ,idat(17),idat(18), newep, sxn,syn)

      END   ! subr. <apply_pm>

************************************************************************

