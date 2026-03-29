c   EXPAND - put "include" packages into fortran source
c
c   the files  MACROS1, MACROS2, & MACROS3  are read:
c     the first lines of MACROS1 contain 35 include followed by 10 exclude
c       strings in '':  'cf8lr' 'calpha'  ' '...   up to 8 chars each
c     the remaining lines of the files are include packages of the
c       form
c       <name(up to 20 characters)
c          lines of source
c       <newname
c
c   then unit 5 is read and copied to unit 6:
c     if a line begins with an include string, then the string
c        is removed and the line shifted left that many positions.
c     if the line begins with an exclude string, it is not copied
c     if a line is of the form
c         >name
c       then the include package is copied to the output.
c
c     In the macro files, lines starting with $ are ignored
c
c     the case of all names is ignored
c
c     4/3/91 - process includes in macro files, MACROS3, 15 strings
c     1/18/92 - expand tab keys to 7, 10, 13, ... (India)
c     9/8/92 - 25 include prefixes, 10 exclude prefixes
c     28/1/95 - 35 include prefixes, 10 exclude prefixes
c
c     NUMPREF is the number of include prefixes.  It must be .GE. NUMEXCL
c             which is the number of exclude
c
      parameter ( numpref = 35 ,  numexcl = 10,  nummacfile = 3 )
      character*20 mnames(200), mname
      integer mstrts(200), mends(200)
      character*80 line, line2, mlines(10000)
      character*8 prefs(numpref,2), pref
      integer prefln(numpref,2), numprf(2)
      character*24 tolowr
      character*30 fnames(nummacfile), fname, sname, sfname
      data fnames / 'macros1', 'macros2', 'macros3' /, sname /' '/
      logical badsw
      character*8 prefou
c
      parameter ( logun = 0 )
c
      nummac = 0
      mline = 0
      badsw = .false.
      lineno = 0
c
c   read macros1, macros2, & macros3
c
      do 199 ifile = 1, nummacfile
         call findfile ( fnames(ifile), fname )
         open ( unit=2, file=fname, status='old' )
         if ( ifile .eq. 1 )  then
            read (2, *) ( prefs(i,1), i = 1, numpref )
            read (2, *) ( prefs(i,2), i = 1, numexcl )
            do 79 j = 1, 2
               npref = numpref
               if ( j .eq. 2 )  npref = numexcl
               do 69  i = 1, npref
                  prefs(i,j) = tolowr(prefs(i,j))
                  k = index(prefs(i,j), ' ') - 1
                  if ( k .eq. -1 )  k = 8
                  prefln(i,j) = k
                  if ( k .eq. 0 )  go to 75
  69           continue
  70           i = npref+1
  75           numprf(j) = i-1
  79        continue
ccc        write (logun, *) numprf, prefln, prefs
         endif
c
c    read the macro lines
c
         assign 100 to iskipln
         assign 110 to ifound
         assign 120 to inopref
c
 100     read (2, 103, end=199) line2
         call extab ( line2, line )
 103     format ( a80 )
 110     if ( line(1:1) .ne. ' ' )  go to 800
 120     if ( line(1:1) .eq. '$' )  then
         else if ( line(1:1) .eq. '<' )  then
            nummac = nummac+1
            mnames(nummac) = tolowr( line(2:25) )
            mstrts(nummac) = mline+1
            mmline = 0
ccc         write (logun, *) ' macro ', mnames(nummac), ' starts at ',
ccc  1         mstrts(nummac)
         else
            mline = mline+1
            mmline = mmline+1
            call chklen ( lineno, line, badsw )
            if (mmline .eq. 1)  line(73:80) = mnames(nummac)(1:8)
            mlines(mline) = line
            mends(nummac) = mline
         endif
         go to 100
 199  continue
c
      if ( badsw ) stop 9876
c
c    now read the file to be expanded
c
      assign 300 to iskipln
      assign 310 to ifound
      assign 600 to inopref
c
c
 300  read ( 5, 103, end=900 ) line2
      call extab ( line2, line )
      lineno = lineno+1
      prefou = ' '
c
 310  if ( line(1:1) .eq. ' ' )  then
         go to 600
c
      else if ( line(1:1) .eq. '>' )  then
c
c    insert macro definitions
c
         mname = tolowr(line(2:25))
ccc          write (logun, *) ' looking', mname
         do 439 imac = 1, nummac
            if ( mnames(imac) .eq. mname ) then
ccc            write(logun, *) ' found' , imac, mstrts(imac), mends(imac)
               do 429 i = mstrts(imac), mends(imac)
                  call writel ( mlines(i) )
 429           continue
               go to 300
            endif
 439     continue
         write (logun,  *) mname, ' is not a defined macro'
         stop 1234
      else
         go to 800
      endif
c
      go to 300
c
 600  call chklen ( lineno, line, badsw )
      if ( mod(lineno,10) .eq. 0 )  then
         write (6, 603) line, lineno
 603     format ( a80, t73, i5 )
      else
         line(73:80) = prefou
         call writel (line)
      endif
      go to 300
c
c     internal subroutine to look for prefix strings
c
 800  pref = tolowr(line(1:8))
ccc      write (logun, *) ' possible pref ', pref
      do 549 j = 1, 2
         do 539  i = 1, numprf(j)
            if ( pref(1:prefln(i,j)) .eq. prefs(i,j) )  then
ccc               write (logun, *) 'found', j, i, pref(1:prefln(i,j))
               if ( j .eq. 1 ) then
                  line2 = line(prefln(i,j)+1:80)
                  line = line2
                  prefou = prefs(i,1)
                  go to ifound, ( 110, 310 )
               else
                  go to iskipln, ( 100, 300 )
               endif
            endif
 539     continue
 549  continue
ccc       write (logun, *) ' not found'
      go to inopref, ( 120, 600 )
c
c
 900  if ( badsw ) stop 9876
      stop
      end
c
c     convert string to lower case
c
      function tolowr ( in )
      character*24 tolowr
      character*(*) in
      character*24 cdum
      character*1 in1(24)
      equivalence ( cdum, in1(1) )
c
      cdum = in
      do 9 i = 1, len(in)
         if ( in1(i) .ge. 'A' .and. in1(i) .le. 'Z' ) then
            j = ichar(in1(i))
            j = j+32
            in1(i) = char(j)
         endif
 9    continue
      tolowr = cdum
      return
      end
c
c     find last non-blank character in line and output the line
c
      subroutine writel( line )
      character*1 line(80)
      do 90 i = 80, 1, -1
         if ( line(i) .ne. ' ' )  then
            go to 100
         endif
 90   continue
      i = 1
 100  write ( 6, 103 ) ( line(j), j = 1, i )
 103  format ( 80a1 )
      return
      end
c
c     search through parent directories for a file
c
      subroutine findfile ( fin, fout )
      character*(*) fin, fout
      character*50 fname, fname2
      logical ex
      parameter ( logun = 0 )
c
      fname = fin
      do 29 i = 1, 5
ccC         write (logun, *) "looking for ", fname
         inquire ( file=fname, exist = ex )
         if ( ex )  go to 100
         fname2 = '../' // fname
         fname = fname2
 29   continue
      write (logun, *) "Could not find ", fin
      stop 999
c
 100  fout = fname
      return
      end
c
c
c
      subroutine extab ( in, out )
c
c     expand tabs to   7, 10, 13, ...
c
      character*80 in, out
      character*1 tab /'	'/
      j = 1
      out = ' '
      do 99 i = 1,80
         if ( in(i:i) .eq. tab )  then
            if ( j .lt. 7 )  then
               j = 7
            else
               j = 3*( (j-7)/3 ) + 10
            endif
         else
            if ( j .lt. 81 ) then
               out(j:j) = in(i:i)
               j = j+1
            endif
         endif
 99   continue
      return
      end
c
c
c
c     make sure there is nothing in 73-80 that is significant
c
      subroutine chklen ( lineno, line, badsw )
      logical badsw
      character*80 line
      if ( line(1:1) .ne. 'c' .and. line(1:1) .ne. 'C'
     1   .and. line(1:1) .ne. '$' .and. line(1:1) .ne. '!' 
     2    )  then
         if ( line(73:80) .ne. ' ' )  then
            write (0,*) 'line ', lineno, ' is too long:'
            write (0,*) line
            badsw = .true.
         endif
      endif
      return
      end
