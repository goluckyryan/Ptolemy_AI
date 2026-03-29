C*ID* phiffer     PTOLEMY LINKULE
ccccni%'linkul' = 'phiffer'
      subroutine phiffer ( alias, savint, ipotyp, ireque, iretur,
     1   iunit, numout, l, jp, rstart, stepsz, numpts,
     2   array1, array2, flt, tempvs, param, intger, jblock,
     3   iswtch, intrnl, fintrn, cnstnt, wavbl, fkandm, id,
     4   alloc, illoc, facfr4, loc, length, nalloc, namloc, iparam )
 
!     linkule for phiffer-computed bound states
                                                                           10   
!     this linkule may be used in either of two ways:
 
!     to calculate both the wavefunction and the potential, the
!     minimum input is
!        "channel:  a + b = c"
!        (or an appropriate "reaction" and "projectile" or "target")
!        "r0 = 1   a = 0.5  E = "
!        "wavefunction = phiffer"
!        "l = xxx   nodes = ... ;"
!     then the linkule returns the wavefunction in array1, and             20   
!     (-1/v)*the potential in array2.
 
!  *** Following not available ***
!           to calculate only the effective potential, and let bound
!     calculate the wave function, use "realpotential" instead of
!     "wavefunction" in the input above".
!     then the linkule returns array2(1)*the potential in array1.
 
 
!           normalization:  the integral of the square of the wave         30   
!     function returned is 1.0.
 
!     2/13/07 - new linkule based on av18 somewhat
 
 
!     relevant arguments are:
 
!     alias = linkule name, used in error messages (input).
!     savint = a two-element real*4 array to store the integrals
!        of wf**2 and wf*r**(l+1).                                         40   
!     ipotyp = 6 (both wf and potential) or 1 (potential only) (input).
!     ireque = 1 for initialization, 2 for printing, 3 for calculation
!        (input).
!     iretur = error return code:  <0 = error (output ).
!     iunit = unit no. for printed output (input).
!     numout = no. of lines to be printed (output).
!     l = orbital angular momentum (input).
!     jp = 2*projectile "total" angular momentum; ?????? set to
!        1 or 3 for l = 0 or 2 respectively (output).
!     rstart = starting r value (fermis) (input).                          50   
!     stepsz = grid spacing (fermis) (input).
!     numpts = no. of grid points to be calculated (input).
!     array1 (length numpts) = primary output array, set to
!        array2(1)*the potential if ipotyp=1, or to the
!        wave function if ipotyp=6.
!     array2: if ipotyp=1, array2(1) is the factor by which the linkule
!        must multiply the potential.
!        if ipotyp=2, array2 (length numpts) is set to
!        (-1/v) * the potential (output).
!     flt(62) = v, input initial guess; set to depth found by phifer (ou   60.  
!     flt(43) = r, input
!     flt(1) = a, input
!     flt(12) = e,input
!     flt(24) = am = reduced mass in mev/c**2 (input).
!     flt(53) = spam = spectroscopic amplitude.
!     flt(54) = spamp = projectile spam. ???  set to spam if next = 1.
!     flt(55) = spamt = target spam.  ???  set to spam if next = 2.
!     param(1) = rho(wine bottle) input
!     param(2) = alpha(wine bottle) input
!     intger(12) = nodes, input                                            70   
!     intger(17) = iprint = print control.  dumps if the 1 digit
!        is >= 6.
!     jblock(1) = 2*total angular momentum,  ???  set to 2 (output).
!     jblock(8) = 2*projectile spin,  ???  set to 1 (output).
!     jblock(9) = 2*target spin,  ???  set to 1 (output).
!     iswtch(11) = next = 1 if this is the projectile, 2 if target.
!     intrnl(3) = notdef = constant for undefined integers.
!     fintrn(1) = undef = constant for undefined real numbers (input).
!     cnstnt(1) = pi (input).
!     cnstnt(2) = sqrt(4*pi) (input).                                      80   
!     cnstnt(6) = hbarc = hbar*c (mev fm) (input).
!     cnstnt(13) = biglog = natural log of a big number whose square
!        doesn's quite overflow (input).
 
 
      implicit real*8 ( a-h, o-z )                                      implicit
      character*8 alias
      character*1 sd(2) / 'S', 'D' /
      real*4  savint(2)
                                                                           90   
      dimension  array1(numpts), array2(numpts), flt(152),
     1   tempvs(6), param(20), intger(50), jblock(12), iswtch(23),
     2   intrnl(74), fintrn(37), cnstnt(13), wavbl(240), fkandm(26),
     3   alloc(1), illoc(1), loc(1), length(1), iparam(5)
ccc!      external nalloc, namloc
 
      logical wfsw
      real*8 params(14)
      real*8, allocatable :: rgrid(:)
                                                                          100   
!
 
!  for each case we have
!     ls(i) = l value
!     phis(,i) and vs(,i) on the grid; phis is normed to 1
!     spams(i) = spec amplitude with sign
!     ebounds(i) = (negative) energy of the state
!     phi_tails(,i) = A, kappa: tail of phi = A h_l( kappa r )
!     v_tails(,i) = A, kappa: nuclear tail of v = a exp(-kappa r)/(kappa
!     v_couls(i) - tail of V also has  v_couls(i)/r                       110   
 
!junk      parameter ( num_case = 2 )
 
!junk      integer  ls(num_case)
!junk      real*4 phis( 0:num_grid-1, num_case ), spams(num_case),
!junk     &   vs( 0:num_grid-1, num_case ), ebounds(num_case),
!junk     &   phi_tails(2, num_case), v_tails(2, num_case), v_couls(num_
 
 
      SAVE                                                                120   
 
!     skip to setup, printing, or calculation code.
 
      iretur = 0
      numout = 0
 
!   get the case number which points to the tables above
 
!junk      icase = l/2 + 1
      undef = fintrn(1)                                                   130   
      notdef = intrnl(3)
 
      if ( ireque - 2 ) 100, 200, 300
 
!     setup:  check for errors
 
 100  if ( ipotyp /= 6 ) then
         write ( iunit, * ) '**** bad ipotyp in ', alias, ipotyp
         numout = 1
         iretur = -1                                                      140   
         return
      endif
 
!     set jp, nodes, v, r, a, e, j, jsp, jst, spam
 
      if ( l == notdef .or. jp == notdef ) then
         write (iunit,*) 'L or JP not defined:', l, jp
         numout = 1
         iretur = -1
         return                                                           150   
      endif
 
 
      nodes = intger(12)
      v = flt(62)
      if ( v == undef )  v = 60
      rv = flt(43)
      rc = flt(45)
      av = flt(1)
      e = flt(12)                                                         160   
      am = flt(24)
      if ( e == undef .or. am == undef  ) then
         write (iunit,*) 'E or red. mass undefined:', e, am
         numout = 1
         iretur = -1
         return
      endif
      if ( rv == undef .or. av == undef ) then
         write (iunit,*) 'R or A undefined:',  rv, av
         numout = 1                                                       170   
         iretur = -1
         return
      endif
      if ( rc == undef ) rc = rv
      rhowb = param(1)
      alphawb = param(2)
      if ( alphawb == undef ) alphawb = 0
      if ( alphawb == 0 ) rhowb = 1
      if ( rhowb == undef ) then
         write (iunit,*) 'rho(wb) undefined, rho, alpha:', param(1:2)     180   
         numout = 1
         iretur = -1
         return
      endif
 
      izp = intger(24)
      izt = intger(25)
      if ( izp == notdef .or. izt == notdef ) then
         write (iunit,*) 'zp or zt not defined', izp,izt
         numout = 1                                                       190   
         iretur = -1
         return
      endif
 
      vso = flt(64)
      rso = flt(49)
      aso = flt(7)
      if ( vso == undef ) vso = 0
      if ( vso /= 0 ) then
         if ( rso == undef .or. aso == undef ) then                       200   
            write (iunit,*) 'rso or vso undef" vso, rso, aso =',
     &         vso, rso, aso
            numout = 1
            iretur = -1
            return
         endif
      else
         if ( rso == undef ) rso = 1
         if ( aso == undef ) aso = 1
      endif                                                               210   
 
 
!junk       jblock(1) = 2
!junk       jblock(8) = 1
!junk      jblock(9) = 1
!     the default spec. amplitudes
!junk      if ( flt(53) .eq. undef ) flt(53) = spams(icase)
!     "spamp" or "spamt", as appropriate, is always set to "spam".
!junk      i = iswtch(11)
!junk      if ( i .eq. 1  .or.  i .eq. 2 )  flt(53+i) = flt(53)           220   
      return
 
 
!     ireque = 2:  print stuff about the linkule.
 
 200  if ( ipotyp == 6 ) then
 
!     printout for wave function calculation.
 
         write (iunit,*) "Phiffer calculation of wave function"           230   
         write (iunit,'(a,3i4,a)') "L, nodes, jp =", l, nodes, jp, "/2"
         write (iunit,'(a,3f10.3)') "E, mu =", e, am
         write (iunit,'(a,3f10.3)') "V(guess), V(convrg) =",
     &      v, -params(10)
         write (iunit,'(a,3f10.3)') "R, A =", rv, av
         write (iunit,'(a,3f10.3)') "w.b. rho, alpha =", rhowb, alphawb
         write (iunit,'(a,3f10.3)') "Vso, Rso, Aso =", vso, rso, aso
         write (iunit,'(a,2i4,3f10.3)') "Zp, Zt, Rc =", izp, izt, rc
 
         numout = 7                                                       240   
 
      else
 
!     printout for potential calculation
         write ( iunit, * ) "******* we should not be here****"
         numout = 1
 
      endif
 
!     printout for both wf and potential calculations.                    250   
 
!junk      write (iunit, 273 ) flt(53)
!junk 273  format ( 10x, 'spectroscopic amplitude ("SPAM") =', f8.5 )
!junk      numout = numout + 1
      return
 
 
!     set up for phifer call
 
 300  am = flt(24)                                                        260   
      e = flt(12)
      rt4pi = cnstnt(2)
      hbarc = cnstnt(6)
 
      hb2o2m = hbarc**2 / (2*am)
      afine = cnstnt(8)
      alpha_coul = izp*izt*hbarc/afine
 
      wfsw = ipotyp .eq. 6
                                                                          270   
      rt = rstart
      rend = rstart + (numpts-1)*stepsz
 
      if ( rt /= 0 ) then
         write (6,*) "rstart not 0:", rt
         iretur = -1
         return
      endif
 
      params(:) = 0                                                       280   
      params(1) = rv
      params(2) = av
      params(3) = e
      params(4) = rhowb
      params(5) = alphawb
      params(6) = rso
      params(7) = aso
      params(8) = vso
      params(10) = v
      params(11) = alpha_coul                                             290   
      params(12) = rc
 
      allocate ( rgrid(numpts+1) )
 
      itype = 11
      ievopt = 1
 
!debug      write (6,*) 'calling phifer:', itype, ievopt, nodes,
!debug     &   l, jp, hb2o2m
!debug      write (6,*) 'Zp,...', izp, izt, hbarc, afine                  300   
!debug      write (6,*) 'params', params(1:12)
!debug      write (6,*) 'numpts...', numpts, rend
 
      CALL PHIFER ( ITYPE, IEVOPT, NODES, L, jp, hb2o2m,
     &   params, .true., numpts, rend, rgrid,
     &   array1, numpts, array2, x )
 
!  phifer returns  phi/r^L .  we want u
 
!debug      write (6,*) 'converting r**(L+1)', l, numpts                  310   
!debug      write (6,*) 'rgrid', rgrid(1:5)
!debug      write (6,*) 'phi ', array1(1:5)
      array1(1:numpts) =
     &      rgrid(1:numpts)**(L+1) * array1(1:numpts)
!debug      write (6,*) 'done', array1(1:5)
 
      flt(62) = -params(10)
!debug      write (6,*) 'v returned =', flt(62)
 
      deallocate ( rgrid )                                                320   
      return
      end
 
 
c*id* phifer
!!!!
!!!!  call stops commented out for speakeas !!!!!
!!!!
 
!!!!!  PTOLEMY version returns potential in  vongrid                      330   
 
      subroutine phifer ( itype, ievopt, nodes, l, jtwo,
     &   h2o2mu, params, prntsw, npts, rmax, rgrid, phi_by_rl,
     &   nvgrid, vongrid, vrgrid )
!
!    generates  single-particle radial function.
!
!    returns phi_by_rl = phi / r**L
!
!    integral r^2 phi^2 dr = 1                                            340   
!
!     itype (decimal format htu) determins the form of correlation:
 
!  u defines the nuclear potential being solved:
!     1 = woods-saxon with wine bottle and/or fermi :
!         v(r) = wsv*{ [1/(1+exp((r-wsr)/wsa))]*[1 + alpha_fermi*(r/wsr)
!                     - alpha*exp(-(r/rho)**2) }
!     2 = shifted gaussian:
!         v(r) = wsv*exp(-((r-wsr)/wsa)**2)
!     8 = output is a normalized h.o. solution for 0 nodes:               350   
!        phi = N r^l exp( -(r/a)^2 )
!        N = sqrt(2**(2*l+3)/((2*l+1)!! * a**(2*l+3) * sqrt(pi/2)) )
!     9 = potential shape input on grid:   wsv*vongrid
 
!  t controls the Coulomb potential
!     0 = no Coulomb
!     1 = The Coulomb pot is
!              alpha_coul * (3 - (r/Rc)**2)/(2Rc)    r < Rc
!              alpha_coul * 1/r                      r > Rc
!     2 = Coulomb is folding of dipole formfactors for both               360   
!         nuclei, with rms radii ra & rb.  strength is alpha_coul
!
!   alpha_coul must be product of atomic numbers times proton charge**2
!              = Z_a*Z_b * hbarc * alpha(fine structure)
!              = Z_a*Z_b * 197.327/137.036
 
!  h controls the boundary conditions
!     0 = Phi goes to exp(-kappa*r)/(kappa*r) at Rmax
!     1 = Phi goes to 0 at Rmax (Infinite well solution)
!     2 = Phi goes to Coulomb (not implimented)                           370   
!     3 = Scattering state solution, E > 0.  params(14) gives the
!         log-deriv of Phi = u/r at RMAX.  if Params(14)=0, then
!         phi(rmax) = 0.
 
!  In addition,
 
!   if L > 0, then a centrifugal potential is added
 
!   if jtwo > 0,  the spin-orbit potential is
!         2*Vso * L.sigma * 1/r * df/dr                                   380   
!         f = 1 / ( 1 + exp( (r - Rso)/aso ) )
 
!    wsv is adjusted to produce given single particle binding wse.
!    phi is provided on a grid with npts # of points up to
!    distance rmax, and is passed with grid values rgrid
!
!   The energy is WSE = params(3).  This is both input and output.
!   For bound states (h=0,1,2), E_sep = WSE > 0; E = -E_sep < 0.
!   for scattering states (h=3), E = WSE > 0.
!                                                                         390   
!     param      value
!      (1)        wsr
!      (2)        wsa
!      (3)        wse
!      (4)        rho
!      (5)        alpha
!      (6)        rso
!      (7)        aso
!      (8)        vso
!      (9)        alpha_ferm                                              400   
!     (10)        wsv
!     (11)        alpha_coulomb
!     (12)        Rc or ra
!     (13)        rb
!     (14)        log deriv or 0 for scattering (h=3)
!
!     ievopt = 0 means vary wse to fit the given wsv
!            = 1 means vary wsv to fit the given wse
!              9 means no solution is found; the complete potential
!                (nuclear + Coulomb + Centrifugal) is returned in PHI.    410   
!                The input value of wsv is used.
!     the modified wse or wsv is returned.
!
!     nodes - number of nodes
!     l - l value
!     jtwo - 2*j for s.o. pot or -1 if no s. o.
!
!     params - the array of parameter values described above.
!
!     prntsw - LOGICAL - if true printing occurs (error messages are      420   
!          always printed)
!
!     npts - number of points to make on the grid.
!     rmax - last r-point; first is 0; stepsize = rmax/(npts-1)
!
!     rgrid is returned as  0, stepsize, ..., rmax
!     phi_by_rl is returned as    phi(r)/r**l
!
!     phi is normalized so that
!         integral(0 rmax) dr r**2 psi(r)**2 = 1                          430   
!
!     11/29/90 - add s.o. pot
!     12/6/90 - add cases 40&46
!     12/7/93 - add 2'd W.S.
!     11/5/99 - massive simplification
!     10/10/00 - use allocate to have no built in dimensions
!     4/19/01 - scattering states
!     4/19/04 - h.o. solutions for 0 nodes
!
      IMPLICIT real*8 ( a-h, o-z )                                        440   
      logical prntsw, infwellsw, varye, scatsw,
     &   eminsw, emaxsw, gridsw
      dimension rgrid(npts), phi_by_rl(npts), params(14)
      real*8, dimension(nvgrid) :: vongrid, vrgrid
!
      dimension nodval(2), xval(2)
!
!
      real*8, allocatable, dimension(:) :: ws, psir, centrg, splb,
     &   splc, spld                                                       450   
!
!
 
      if (prntsw)  write(6, 7) itype, ievopt, nodes, l,
     1     jtwo, npts, rmax, h2o2mu, params(1:13)
 7    format ('phifer input:', i5, 3i3, i4, '/2', i6, f6.1, f10.4 /
     &    ( 1x, 8f10.4 ) )
 
      maxgrid = npts
      if ( gridsw ) maxgrid = max( maxgrid, nvgrid )                      460   
      allocate ( ws(maxgrid), psir(maxgrid), centrg(maxgrid),
     &   splb(maxgrid), splc(maxgrid), spld(maxgrid) )
 
      iu = mod( itype, 10 )
      it = mod( itype/10, 10 )
      ih = mod( itype/100, 10 )
      gridsw = iu == 9
      infwellsw = ih == 1
      scatsw = ih == 3
      varye = ievopt .eq. 0                                               470   
!
      sml=1.e-10
      dr=rmax/(npts-1)
      dx=dr*dr/(12.*h2o2mu)
      nptsm=npts-1
      nptsm2=npts-2
      do i = 1, npts
         rgrid(i) = dr*(i-1)
      enddo
                                                                          480   
      wsr = params(1)
      wsa = params(2)
      wse = params(3)
      rho = params(4)
      alpha = params(5)
      rso = params(6)
      aso = params(7)
      vso = params(8)
      alpha_fermi = params(9)
      wsv = params(10)                                                    490   
      alpha_coul = params(11)
      rc = params(12)
      rb = params(13)
      xa = rc/sqrt(12.)
      xb = rb/sqrt(12.)
      if ( iu == 8 ) go to 800
      if (scatsw) deriv = params(14)
      if ( rho .eq. 0 )  rho = 1
      dl=l*(l+1)*dr*dr/12
      consso = 0                                                          500   
      if ( jtwo .ne. -1 )  consso =
     1   -2.*vso*(jtwo*(jtwo+2)/4. - l*(l+1) - .75)/aso
!
      if ( gridsw ) then
         call splncb ( nvgrid, vrgrid, vongrid, splb, splc, spld )
         call intrpc ( nvgrid, vrgrid, vongrid, splb, splc, spld,
     1      npts, rgrid, ws )
      endif
!
      do  i = 1, npts                                                     510   
         r = dr*(i-1)
 
         select case (iu)
         case (1)
            ws(i) = 1/(1.+exp((r-wsr)/wsa))
            if ( wsr > 0 ) ws(i) = ws(i)*
     &         ( 1 + (alpha_fermi/wsr**2)*r**2 )
            ws(i) = ws(i) - alpha*exp(-(r/rho)**2)
         case (2)
            ws(i) = exp(-((r-wsr)/wsa)**2)                                520   
         case (9)
         case default
            call stop ( 'PHIFER: invalid itype', itype )
         end select
 
         ws(i) = dx*ws(i)
         centrg(i)=dl/(r+1.e-20)**2
         if ( consso .ne. 0 )  then
            x = exp((r-rso)/aso)
            so = consso*x/((r+1.e-20)*(1+x)**2)                           530   
            centrg(i) = centrg(i) + dx*so
         endif
 
         select case (it)
         case (0)
            vc = 0
         case (1)
            if ( r .gt. rc )  then
               vc = alpha_coul/r
            else                                                          540   
               vc = (alpha_coul/(2*rc)) * (3 - (r/rc)**2)
            endif
         case (2)
            if (r > 1e-3) then
               vc = alpha_coul*( 1 - ( .5/(1-(xa/xb)**2)**2 )
     &              *exp(-r/xb)*( 2 + (r/xb) + 4/(1-(xb/xa)**2) )
     &            - ( .5/(1-(xb/xa)**2)**2 )
     &              *exp(-r/xa)*( 2 + (r/xa) + 4/(1-(xa/xb)**2) ) )/r
            else
               vc = (alpha_coul/(2*(xb**2-xa**2)**3))                     550   
     &           *( xb**3*(xb**2-5*xa**2) + xa**3*(5*xb**2-xa**2) )
            end if
         case default
            call stop ( 'PHIFER: invalid itype', itype )
         end select
         centrg(i) = centrg(i) + dx*vc
 
      enddo
 
      if ( ievopt == 9 ) then                                             560   
         do i = 1, npts
            phi_by_rl(i) = ( wsv*ws(i) + centrg(i) )/dx
         enddo
         go to 900
      endif
!
      nodlop = 1
      efac0 = 1.3
      efac = efac0
      nodelast = 9999                                                     570   
      nodhave = 0
      eminsw = .false.
      emaxsw = .false.
      if ( .not. infwellsw ) then
         if ( .not. gridsw ) wsv = -abs(wsv)
         wse = abs(wse)
      endif
 
!   loop back
                                                                          580   
 15   estrt = wse
      vstrt = wsv
      vme = ( wsv*ws(1)/dx + wse )/h2o2mu
      match=int( sqrt( -l*(l+1)/vme )/dr + .5 )
      match=max(match,npts/20)
      matchp=match+1
!
!     converge on the eigenvalue
!
      do 80 loop = 1, 100                                                 590   
!
!     integrate from outside in
!     construct function at last two points from boundary condition
!
         ak2 = -dx*wse
         akap=sqrt(+wse/h2o2mu)
         if ( infwellsw ) then
            psir(npts) = 0
            psir(nptsm) = .01
         else if ( scatsw ) then                                          600   
            ak2 = -ak2
            if ( deriv == 0 ) then
               psir(npts) = 0
               psir(nptsm) = .01
            else
!  log deriv is for  phi = u/r.  2nd deriv assumed to be just
!  energy and angular momentum, potential = 0
               psir(npts) = 1
               d1 = deriv + 1./rmax
               d2 = -akap**2 + l*(l+1)/rmax**2                            610   
               psir(nptsm) = 1 - d1*dr + .5*d2*dr**2
            endif
         else
            psir(npts)=exp(-akap*rgrid(npts))
            psir(nptsm)=exp(-akap*rgrid(nptsm))
         endif
         xa=psir(npts)
         xb=psir(nptsm)
         xd=wsv*ws(npts)+centrg(npts)-ak2
         xe=wsv*ws(nptsm)+centrg(nptsm)-ak2                               620   
         do  j = nptsm2, match, -1
            xf=wsv*ws(j)+centrg(j)-ak2
            psir(j)=((2.+10.*xe)*xb-(1.-xd)*xa)/(1.-xf)
            xa=xb
            xb=psir(j)
            xd=xe
            xe=xf
         enddo
         gg1=xb
         gg2=xa                                                           630   
!
!      integrate from inside out
!
         psir(1)=0.
         psir(2)=(1.+12.*(wsv*ws(2)-ak2)/(4*l+6))*dr**(l+1)
         xb=psir(2)
         xe=wsv*ws(2)+centrg(2)-ak2
         xf=wsv*ws(3)+centrg(3)-ak2
         psir(3)=(2.+10.*xe)*xb/(1.-xf)
         if (l.eq.1) psir(3)=psir(3)+2.*h2o2mu*dx/(1.-xf)                 640   
         xa=xb
         xb=psir(3)
         xd=xe
         xe=xf
         do  j = 4, matchp
            xf=wsv*ws(j)+centrg(j)-ak2
            psir(j)=((2.+10.*xe)*xb-(1.-xd)*xa)/(1.-xf)
            xa=xb
            xb=psir(j)
            xd=xe                                                         650   
            xe=xf
         enddo
         ff1=xa
         ff2=xb
         dldif=(ff2*gg1-ff1*gg2)/(abs(ff1*gg1)+abs(ff2*gg2))
 
         if ( loop .eq. 1 ) then
            dldif2=dldif
            wse2=wse
            wsv2=wsv                                                      660   
            if (varye) then
               wse=wse2*1.2
            else
               wsv=wsv2*1.2
            endif
            go to 80
         else if ( loop .eq. 2 ) then
            dldif1=dldif
            wse1=wse
            wsv1=wsv                                                      670   
         else
            if (abs(dldif) .le. sml) go to 90
            if (dldif*dldif1) 54,54,56
 54         if (dldif*dldif2) 58,58,60
 56         if (dldif*dldif2) 65,65,58
 58         if (abs(dldif1)-abs(dldif2)) 60,60,65
!
 60         dldif2=dldif1
            wse2=wse1
            wsv2=wsv1                                                     680   
            dldif1=dldif
            wse1=wse
            wsv1=wsv
            go to 70
 65         dldif1=dldif
            wse1=wse
            wsv1=wsv
         endif
!
 70      if (varye) then                                                  690   
            wse=(dldif2*wse1-dldif1*wse2)/(dldif2-dldif1)
            if ( wse .lt. 0. .and. .not. infwellsw )  then
               wse = .5*min(wse1, wse2)
               if ( wse .lt. .001 ) then
                  if (prntsw)
     &               write (6,*) ' wse tried for 0, start=', estrt
                  wse = 2*estrt
                  go to 15
               endif
            endif                                                         700   
            if ( min(abs(wse-wse1), abs(wse-wse2))/abs(wse1-wse2)
     &         .lt. .1 ) wse = .5*(wse1+wse2)
            if ( eminsw ) then
               if ( wse .lt. emin )  then
                  wse = .5*(emin+min(wse1,wse2))
               endif
            endif
            if ( emaxsw ) then
               if ( wse .lt. emax )  then
                  wse = .5*(emax+max(wse1,wse2))                          710   
               endif
            endif
         else
            wsv = (dldif2*wsv1-dldif1*wsv2)/(dldif2-dldif1)
!  don't let wsv change sign, but it could be either sign
            if ( wsv*vstrt .le. 0.)  wsv =
     &         sign( .5*min( abs(wsv1), abs(wsv2) ), vstrt )
            if ( min(abs(wsv-wsv1), abs(wsv-wsv2))/abs(wsv1-wsv2)
     &         .lt. .1 ) wsv = .5*(wsv1+wsv2)
         endif                                                            720   
 80   continue
!
      if ( abs(dldif) .gt. 10000*sml ) then
         write(6, 88) wse, wsv, wse1, wsv1, dldif1,
     1      wse2, wsv2, dldif2, ff1, gg1, ff2, gg2, sml
 88      format ( '0*********** did not converge on wavefunction ',
     1      20('*') / 2g25.15 / 3g25.15 / 3g25.15 / 4g25.15 /  4g25.15)
!!!!         call stop ( 'phifer - did not converge on wavefunction', 0
         write (6,*) '******** phiffer could not converge ********'
      else                                                                730   
         if (prntsw) write ( 6, '("phifer accepts poor converg",
     &      2g15.5)' ) dldif, sml
      endif
!
!     normalize inside r(matchp)
!
 90   fac=gg2/psir(matchp)
      do  i = 2, matchp
         psir(i)=fac*psir(i)
      enddo                                                               740   
!
!     normalize and count nodes
!
      sum1=0.
      sum2=0.
      f = psir(4)
      nod = 0
      do  i = 2,npts
         sum1=sum1+psir(i)**2
         sum2=sum2+psir(i)**2*rgrid(i)**2                                 750   
         phi_by_rl(i)=psir(i)/rgrid(i)**(l+1)
         if ( psir(i)*f <= 0  .and.  i > 4 .and.  i < npts ) then
            nod = nod+1
            f = psir(i+1)
         endif
      enddo
      if ( nod .eq. nodes )  go to 130
!
!     wrong number of nodes; try again
!                                                                         760   
      if (prntsw)  write(6, 113) nod, estrt, wse, efac, vstrt, wsv, loop
 113  format ( ' got', i3, ' nodes:  e =', 2g15.5,' efac=', g15.5,
     1   5x, 'v =', 2f10.5, i5, ' matching loops' )
      nodlop = nodlop + 1
      if ( nodlop .gt. 40 ) then
         write (6, *) ' *** node loop limit reached:',
     &      nod, estrt, wse, efac, vstrt, wsv
         go to 130
!!!!!         call stop ( 'could not get number of nodes', nod )
      endif                                                               770   
!
      if ( varye )  then
!     varing e
         if ( nod .eq. nodelast ) then
            efac = 1.1*efac
         else
            efac = efac0
         endif
         nodelast = nod
!                                                                         780   
!     set barriers on e based on what we have found
!
         xdel = .1*abs(wse)
         if ( nod .gt. nodes ) then
            if ( .not. eminsw  .or.
     &         ( eminsw .and. wse-xdel .gt. emin ) )  then
               eminsw = .true.
               emin = wse - xdel
            endif
         else                                                             790   
            if ( .not. emaxsw  .or.
     &         ( emaxsw .and. wse+xdel .lt. emax ) )  then
               emaxsw = .true.
               emax = wse + xdel
            endif
         endif
!
!     accumulate eigenvalues and number of nodes
!
         if ( nodhave .eq. 0 ) then                                       800   
            nodhave = 1
            nodval(1) = nod
            xval(1) = wse
            go to 120
         else if ( nodhave .eq. 1 ) then
            if ( nod .eq. nodval(1) ) go to 120
            nodhave = 2
            nodval(2) = nod
            xval(2) = wse
         else                                                             810   
!
!     if we got a different number of nodes, we can use it
!
            if ( nod .eq. nodval(1) .or. nod .eq. nodval(2) )
     &         go to 120
!
!     replace the worse number of nodes, but only with a better value
!
            if ( abs(nodval(1)-nodes) .gt. abs(nodval(2)-nodes) ) then
               i = 1                                                      820   
            else
               i = 2
            endif
            if ( abs(nodval(i)-nodes) .gt. abs(nod-nodes) ) then
               nodval(i) = nod
               xval(i) = wse
            else
               go to 120
            endif
         endif                                                            830   
!
!     predict the new eigenval based on number of nodes found
!
         x = xval(1) + ((xval(2)-xval(1))/(nodval(2)-nodval(1)))
     &      * (nodes-nodval(1))
         wse = x
         go to 125
!
!     either we don't have two different solutions, or the prediction
!     based on two leads back to one of them.                             840   
!
 120     if ( nod .lt. nodes ) then
            wse = min( wse, estrt )
            if ( wse .gt. 0 ) then
               wse = wse/efac
            else
               wse = wse*efac
            endif
         else
            wse = max( wse, estrt )                                       850   
            if ( wse .gt. 0 ) then
               wse = wse*efac
            else
               wse = wse/efac
            endif
         endif
      else
!     varing v
         if ( nod .lt. nodes )  then
            wsv = sign( 1.3*max( abs(wsv), abs(vstrt) ), vstrt )          860   
         else
            wsv = sign( min( abs(wsv), abs(vstrt) )/1.2, vstrt )
         endif
      endif
 125  go to 15
!
!     have correct number of nodes
!
 130  phi_by_rl(1) = (rgrid(3)**2*phi_by_rl(2)
     &   - rgrid(2)**2*phi_by_rl(3)) / (rgrid(3)**2-rgrid(2)**2)          870   
      anorm = sqrt(sum1*dr)
      phi_by_rl(1:npts) = phi_by_rl(1:npts)/anorm
      rms=sqrt(sum2/sum1)
      if ( prntsw ) write(6, 503) wsv, wse, rms, loop, dldif,
     1  ( phi_by_rl(i), i = 1, 5 )
  503 format(' woods-saxon strength v =', f15.10, 5x,
     1   'e(sep) =', f15.10, ';   rms radius =', f8.5,
     2    i5, ' iters;' /
     2    ' convrg =', g12.3,
     3    5x, 'phi/r**l:', 5g18.9 )                                       880   
!
!     if we completely failed to converge, stop
!
      IF ( ABS(DLDIF) .GT. 1000*SML ) then
         write (6,*) 'could not converge:', dldif, 1000*sml
!!!         call stop (  'PHIFER: could not converge:', 553 )
      endif
      params(3) = wse
      params(10) = wsv
                                                                          890   
!  PTOLEMY
 
      vongrid(1:nvgrid) = ws(1:nvgrid)/dx
 
      go to 900
 
!  return a h.o. solution
 
 800  if ( nodes /= 0 ) then
!!!         call stop ('H.O. must have 0 nodes', nodes)                   900   
         stop 9876
      endif
      pi = acos(-1.d0)
      anorm = (2./wsa)**(2*l+3) * sqrt(2/pi)
      do ll = 1, 2*l+1, 2
         anorm = anorm/ll
      enddo
      anorm = sqrt(anorm)
      write (6,*) 'H.O.; L=', L, ' A=', wsa, ' Norm=', anorm
      do i = 1, npts                                                      910   
         x = rgrid(i)/wsa
!  we leave out the r**l because phifer returns  phi/r**l
         phi_by_rl(i) = anorm * exp( -x**2 )
      enddo
      return
 
 900  deallocate ( ws, psir, centrg, splb, splc, spld )
      return
      end
c*id* stop                                                                920   
      subroutine stop ( msg, icode )
c
c     print an error message and stop.
c     If compiled for MPI, abort the whole job.
c
cmpi      include 'mpif.h'
      character*(*) msg
c
      write ( 6, 13 )  msg, icode
      write ( 0, 23 )                                                     930   
      write ( 0, 13 )  msg, icode
      write ( 0, 23 )
 13   format ( / 1x, 78('*') /
     &   ' *',20x, 'Stopping because of an error  !!', t79, '*' /
     &   ' * ', a, ' : ', i10,  t79, '*' /
     &   1x, 78('*') / )
 23   format ( '' )
c
cmpi      call MPI_abort ( MPI_COMM_WORLD, icode, ierr )
cmpi      write (0,*) ' MPI_Abort is done', ierr                          940   
c
      stop 9999
      end
