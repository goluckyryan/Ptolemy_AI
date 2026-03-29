C*ID* AV18     PTOLEMY LINKULE
ccccni%'linkul' = 'av18'
      subroutine av18 ( alias, savint, ipotyp, ireque, iretur,
     1   iunit, numout, l, jp, rstart, stepsz, numpts,
     2   array1, array2, flt, tempvs, param, intger, jblock,
     3   iswtch, intrnl, fintrn, cnstnt, wavbl, fkandm, id,
     4   alloc, illoc, facfr4, loc, length, nalloc, namloc, iparam )
 
!     linkule for Argonne v18 deuteron wave function
                                                                           10   
!     this linkule may be used in either of two ways:
 
!     to calculate both the wavefunction and the potential, the
!     minimum input is
!        "channel:  n + p = d"
!        (or an appropriate "reaction" and "projectile" or "target")
!        "r0 = 1   a = 0.5"      (to avoid complaints)
!        "wavefunction = av18"
!        "l = 0;"  (or "l = 2;")
!     then the linkule returns the wavefunction in array1, and             20   
!     (-1/v)*the potential in array2.
!           to calculate only the effective potential, and let bound
!     calculate the wave function, use "realpotential" instead of
!     "wavefunction" in the input above, and add "nodes = 0".
!     then the linkule returns array2(1)*the potential in array1.
!     (note:  the present version of  bound  (10/79) will not calculate
!     the wave fuction from this potential because it chokes on the
!     core.)
 
                                                                           30   
!           normalization:  the integral of the square of the wave
!     function returned is 1.0.  spectroscopic information is entered
!     using the keyword "spam", which is not actually used until much
!     later in the dwba calculation.  if "spam" is not defined, this
!     linkule supplies default values from the actual wave functions:
!     spam(l=0) = sqrt(1-pd), spam(l=2) = sqrt(pd),
 
!           note:  in dwba calculations, the default spectroscopic
!     amplitude is always used unless the keyword "spam" is used in
!     the description of this channel.  setting "spamp" or "spamt" will    40   
!     not work.  in stand-alone calculations, you must "undefine spam"
!     to force use of the default in the second and subsequent
!     calculations.  f "spam" is used, remember that a prolate
!     deuteron requires a negative d-state amplitude.
 
!     method:
!           the wave functions (both s and d) are calculated by
!     interpolating between tabulated values which is good to ~.2%.
!     The potentials are also computed by interpolation
!                                                                          50   
 
!     3/16/03 - new linkule based on Reid somewhat
 
 
!     relevant arguments are:
 
!     alias = linkule name, used in error messages (input).
!     savint = a two-element real*4 array to store the integrals
!        of wf**2 and wf*r**(l+1).
!     ipotyp = 6 (both wf and potential) or 1 (potential only) (input).    60   
!     ireque = 1 for initialization, 2 for printing, 3 for calculation
!        (input).
!     iretur = error return code:  <0 = error (output ).
!     iunit = unit no. for printed output (input).
!     numout = no. of lines to be printed (output).
!     l = orbital angular momentum, must = 0 or 2 (input).
!     jp = 2*projectile "total" angular momentum; set to
!        1 or 3 for l = 0 or 2 respectively (output).
!     rstart = starting r value (fermis) (input).
!     stepsz = grid spacing (fermis) (input).                              70   
!     numpts = no. of grid points to be calculated (input).
!     array1 (length numpts) = primary output array, set to
!        array2(1)*the potential if ipotyp=1, or to the
!        wave function if ipotyp=6.
!     array2: if ipotyp=1, array2(1) is the factor by which the linkule
!        must multiply the potential.
!        if ipotyp=2, array2 (length numpts) is set to
!        (-1/v) * the potential (output).
!     flt(62) = v, set to 1.0 (output).
!     flt(43) = r, set to 1.0 if undefined.                                80   
!     flt(1) = a, set to 0.4 if undefined.
!     flt(12) = e, set to -2.224644 if undefined.
!     flt(24) = am = reduced mass in mev/c**2 (input).
!     flt(53) = spam = spectroscopic amplitude.
!        set to sqrt(1-pd) or sqrt(pd) if l = 0 or 2, respectively,
!        if initially undefined, where pd = 0.064696.
!     flt(54) = spamp = projectile spam.  set to spam if next = 1.
!     flt(55) = spamt = target spam.  set to spam if next = 2.
!     param(1)
!     param(2)                                                             90   
!     intger(12) = nodes, set to 0 (output).
!     intger(17) = iprint = print control.  dumps if the 1 digit
!        is >= 6.
!     jblock(1) = 2*total angular momentum, set to 2 (output).
!     jblock(8) = 2*projectile spin, set to 1 (output).
!     jblock(9) = 2*target spin, set to 1 (output).
!     iswtch(11) = next = 1 if this is the projectile, 2 if target.
!     intrnl(3) = notdef = constant for undefined integers.
!     fintrn(1) = undef = constant for undefined real numbers (input).
!     cnstnt(1) = pi (input).                                             100   
!     cnstnt(2) = sqrt(4*pi) (input).
!     cnstnt(6) = hbarc = hbar*c (mev fm) (input).
!     cnstnt(13) = biglog = natural log of a big number whose square
!        doesn's quite overflow (input).
 
 
      implicit real*8 ( a-h, o-z )                                      implicit
      character*8 alias
      character*1 sd(2) / 'S', 'D' /
      real*4  savint(2)                                                   110   
 
      dimension  array1(numpts), array2(numpts), flt(152),
     1   tempvs(6), param(20), intger(50), jblock(12), iswtch(23),
     2   intrnl(74), fintrn(37), cnstnt(13), wavbl(240), fkandm(26),
     3   alloc(1), illoc(1), loc(1), length(1), iparam(5)
ccc!      external nalloc, namloc
 
      logical wfsw
 
!                                                                         120   
!     the grid defining phi's is r = 0(.1)12.
 
      parameter ( num_grid = 241, grid_step = .05d0 )
      parameter ( grid_max = (num_grid-1)*grid_step )
 
!  for each case we have
!     ls(i) = l value
!     phis(,i) and vs(,i) on the grid; phis is normed to 1
!     spams(i) = spec amplitude with sign
!     ebounds(i) = (negative) energy of the state                         130   
!     phi_tails(,i) = A, kappa: tail of phi = A h_l( kappa r )
!     v_tails(,i) = A, kappa: nuclear tail of v = a exp(-kappa r)/(kappa
!     v_couls(i) - tail of V also has  v_couls(i)/r
 
      parameter ( num_case = 2 )
 
      integer  ls(num_case)
      real*4 phis( 0:num_grid-1, num_case ), spams(num_case),
     &   vs( 0:num_grid-1, num_case ), ebounds(num_case),
     &   phi_tails(2, num_case), v_tails(2, num_case), v_couls(num_case)  140   
 
!  av18 deuteron - S & D
 
      data ls(1:2) / 0, 2 /
 
      data phis(:,1) /
     & .08152, .083498, .089482, .09958, .1139, .13244, .15504,
     & .18127, .21044, .24155, .2734, .30473, .33428, .36098, .38402,
     & .40291, .41745, .42774, .43406, .43685, .43658, .43378,
     & .42891, .42243, .4147, .40605, .39675, .38701, .37699, .36684,     150   
     & .35665, .34652, .3365, .32664, .31698, .30754, .29833, .28938,
     & .28069, .27226, .2641, .2562, .24856, .24118, .23405, .22716,
     & .22051, .21409, .2079, .20192, .19615, .19058, .18521, .18002,
     & .17501, .17018, .16551, .161, .15665, .15244, .14837, .14444,
     & .14064, .13697, .13341, .12997, .12664, .12342, .1203, .11728,
     & .11436, .11152, .10878, .10611, .10353, .10103, .098599,
     & .096243, .093957, .091737, .089583, .087491, .085459, .083485,
     & .081567, .079703, .077891, .07613, .074417, .072751, .071131,
     & .069554, .06802, .066527, .065073, .063658, .06228, .060938,
     & .059631, .058357, .057116, .055907, .054728, .053579, .052459,     160   
     & .051367, .050301, .049262, .048249, .047259, .046295, .045353,
     & .044434, .043537, .042661, .041806, .040971, .040155, .039359,
     & .038581, .037821, .037078, .036353, .035644, .03495, .034273,
     & .033611, .032964, .032331, .031712, .031107, .030515, .029936,
     & .029369, .028815, .028273, .027743, .027224, .026716, .026219,
     & .025733, .025257, .024791, .024334, .023887, .02345, .023022,
     & .022602, .022191, .021789, .021395, .021008, .02063, .02026,
     & .019897, .019541, .019192, .018851, .018516, .018188, .017866,
     & .017551, .017242, .016939, .016642, .016351, .016065, .015785,
     & .015511, .015241, .014977, .014718, .014465, .014215, .013971,     170   
     & .013731, .013496, .013265, .013039, .012817, .012599, .012385,
     & .012175, .01197, .011767, .011569, .011375, .011183, .010996,
     & .010812, .010631, .010454, .01028, .010109, .0099407,
     & .0097759, .0096141, .0094552, .0092991, .0091459, .0089954,
     & .0088476, .0087024, .0085599, .0084199, .0082823, .0081472,
     & .0080145, .0078842, .0077561, .0076303, .0075068, .0073853,
     & .007266, .0071488, .0070336, .0069205, .0068093, .0067,
     & .0065926, .0064871, .0063834, .0062815, .0061813, .0060829,
     & .0059861, .005891, .0057976, .0057057, .0056154, .0055266,
     & .0054393, .0053535, .0052692, .0051862, .0051047, .0050246,        180   
     & .0049457, .0048683, .0047921, .0047172 /
 
      data phis(:,2) /
     & .0, .0052663, .019706, .042113, .072017, .10918, .15328,
     & .20356, .25871, .31689, .37582, .43311, .48642, .53383,
     & .57393, .6059, .62952, .64505, .65313, .6546, .65043, .64161,
     & .62907, .61368, .59616, .57715, .55717, .53665, .51593,
     & .49526, .47487, .45488, .43543, .41658, .39838, .38088,
     & .36408, .34798, .33259, .31789, .30386, .29049, .27774,
     & .26561, .25405, .24305, .23257, .22261, .21312, .20409, .1955,     190   
     & .18732, .17953, .17211, .16504, .15831, .15189, .14578,
     & .13995, .13439, .12909, .12403, .1192, .11459, .11019, .10599,
     & .10197, .098132, .094463, .090953, .087597, .084385, .081312,
     & .07837, .075552, .072854, .070268, .06779, .065415, .063137,
     & .060953, .058857, .056846, .054916, .053063, .051283, .049574,
     & .047931, .046353, .044835, .043376, .041973, .040623, .039324,
     & .038074, .036871, .035712, .034596, .033521, .032485, .031488,
     & .030526, .029598, .028704, .027841, .027009, .026206, .025431,
     & .024683, .023961, .023264, .02259, .021939, .021311, .020703,
     & .020116, .019548, .018999, .018468, .017955, .017458, .016977,     200   
     & .016512, .016062, .015626, .015204, .014795, .014399, .014015,
     & .013644, .013284, .012935, .012596, .012268, .01195, .011642,
     & .011343, .011052, .010771, .010498, .010233, .0099752,
     & .0097254, .0094828, .0092473, .0090185, .0087964, .0085805,
     & .0083709, .0081671, .0079691, .0077767, .0075896, .0074077,
     & .0072309, .007059, .0068918, .0067291, .0065709, .006417,
     & .0062672, .0061215, .0059797, .0058416, .0057072, .0055764,
     & .005449, .005325, .0052042, .0050866, .004972, .0048603,
     & .0047516, .0046456, .0045424, .0044417, .0043436, .004248,
     & .0041548, .0040639, .0039754, .003889, .0038047, .0037225,         210   
     & .0036424, .0035642, .0034879, .0034135, .0033409, .00327,
     & .0032009, .0031334, .0030675, .0030032, .0029404, .0028791,
     & .0028192, .0027608, .0027037, .002648, .0025935, .0025404,
     & .0024884, .0024377, .0023881, .0023396, .0022923, .002246,
     & .0022008, .0021566, .0021134, .0020712, .0020299, .0019896,
     & .0019501, .0019116, .0018738, .001837, .0018009, .0017656,
     & .0017311, .0016973, .0016643, .001632, .0016004, .0015695,
     & .0015392, .0015096, .0014806, .0014523, .0014245, .0013974,
     & .0013708, .0013448, .0013193, .0012943, .0012699, .001246,
     & .0012226, .0011997, .0011772 /                                     220   
 
      data vs(:,1) /
     & 2408., 2368.4, 2247.6, 2055.6, 1812.6, 1541.7, 1263.9, 995.83,
     & 749.81, 533.59, 351.2, 203.56, 89.222, 5.0502, -53.16,
     & -90.118, -110.55, -118.84, -118.74, -113.3, -104.84, -95.015,
     & -84.918, -75.227, -66.31, -58.323, -51.294, -45.172, -39.872,
     & -35.292, -31.332, -27.901, -24.917, -22.312, -20.028, -18.018,
     & -16.24, -14.664, -13.262, -12.011, -10.892, -9.8895, -8.99,
     & -8.1816, -7.454, -6.7984, -6.2072, -5.6734, -5.191, -4.7547,
     & -4.3597, -4.0017, -3.677, -3.3823, -3.1144, -2.8707, -2.6489,      230   
     & -2.4467, -2.2622, -2.0937, -1.9397, -1.7987, -1.6695, -1.5511,
     & -1.4423, -1.3424, -1.2504, -1.1657, -1.0877, -1.0156, -.94907,
     & -.88754, -.8306, -.77785, -.72895, -.68357, -.64142, -.60224,
     & -.56579, -.53185, -.50022, -.47073, -.44322, -.41752, -.3935,
     & -.37104, -.35002, -.33034, -.31191, -.29462, -.27841, -.26319,
     & -.24889, -.23547, -.22284, -.21097, -.19979, -.18927, -.17936,
     & -.17003, -.16122, -.15292, -.14509, -.1377, -.13071, -.12412,
     & -.11789, -.11199, -.10642, -.10115, -.096165, -.091444,
     & -.086974, -.082739, -.078728, -.074926, -.071322, -.067904,
     & -.064662, -.061586, -.058667, -.055896, -.053266, -.050767,        240   
     & -.048394, -.046139, -.043995, -.041958, -.040021, -.038179,
     & -.036428, -.034761, -.033175, -.031666, -.030229, -.028861,
     & -.027559, -.026318, -.025137, -.024011, -.022939, -.021916,
     & -.020942, -.020013, -.019127, -.018283, -.017478, -.016709,
     & -.015976, -.015277, -.01461, -.013973, -.013365, -.012785,
     & -.012231, -.011702, -.011197, -.010714, -.010253, -.0098131,
     & -.0093925, -.0089906, -.0086066, -.0082396, -.0078888,
     & -.0075535, -.007233, -.0069265, -.0066335, -.0063533,
     & -.0060853, -.005829, -.0055839, -.0053494, -.0051251,
     & -.0049104, -.0047051, -.0045086, -.0043205, -.0041405,             250   
     & -.0039683, -.0038034, -.0036455, -.0034944, -.0033497,
     & -.0032112, -.0030786, -.0029516, -.0028299, -.0027134,
     & -.0026018, -.002495, -.0023926, -.0022945, -.0022005,
     & -.0021105, -.0020242, -.0019416, -.0018624, -.0017865,
     & -.0017137, -.001644, -.0015772, -.0015131, -.0014517,
     & -.0013928, -.0013364, -.0012823, -.0012304, -.0011807,
     & -.001133, -.0010873, -.0010434, -.0010014, -9.6106E-4,
     & -9.2238E-4, -8.8528E-4, -8.497E-4, -8.1557E-4, -7.8283E-4,
     & -7.5142E-4, -7.2129E-4, -6.9238E-4, -6.6465E-4, -6.3804E-4,
     & -6.1251E-4, -5.8802E-4, -5.6451E-4, -5.4196E-4, -5.2032E-4,        260   
     & -4.9955E-4, -4.7961E-4, -4.6048E-4, -4.4212E-4, -4.245E-4,
     & -4.0759E-4, -3.9136E-4, -3.7577E-4, -3.6082E-4, -3.4646E-4,
     & -3.3268E-4 /
 
      data vs(:,2) /
     & 1763.5, -7088.1, -2610.3, -1239.7, -678.63, -451.01, -384.97,
     & -398.84, -447.45, -504.36, -554.5, -590.14, -608.53, -610.1,
     & -597.15, -572.89, -540.68, -503.63, -464.35, -424.85, -386.57,
     & -350.43, -316.96, -286.41, -258.8, -234.02, -211.89, -192.19,
     & -174.67, -159.09, -145.25, -132.92, -121.94, -112.12, -103.34,     270   
     & -95.457, -88.365, -81.968, -76.182, -70.934, -66.161, -61.809,
     & -57.831, -54.185, -50.837, -47.755, -44.911, -42.283, -39.85,
     & -37.592, -35.494, -33.542, -31.723, -30.024, -28.437, -26.953,
     & -25.562, -24.257, -23.032, -21.882, -20.799, -19.781, -18.821,
     & -17.915, -17.061, -16.254, -15.492, -14.771, -14.088, -13.442,
     & -12.83, -12.249, -11.699, -11.176, -10.68, -10.209, -9.7612,
     & -9.3353, -8.9301, -8.5446, -8.1775, -7.8279, -7.4947, -7.1772,
     & -6.8745, -6.5858, -6.3103, -6.0474, -5.7963, -5.5566, -5.3277,
     & -5.1089, -4.8998, -4.6999, -4.5088, -4.3259, -4.1511, -3.9837,
     & -3.8236, -3.6702, -3.5234, -3.3829, -3.2482, -3.1192, -2.9957,     280   
     & -2.8772, -2.7637, -2.6549, -2.5506, -2.4506, -2.3546, -2.2626,
     & -2.1743, -2.0896, -2.0084, -1.9304, -1.8556, -1.7837, -1.7148,
     & -1.6486, -1.585, -1.524, -1.4654, -1.4091, -1.355, -1.3031,
     & -1.2532, -1.2053, -1.1592, -1.115, -1.0725, -1.0316, -.99234,
     & -.95459, -.91831, -.88344, -.84991, -.81768, -.7867, -.75691,
     & -.72826, -.70072, -.67423, -.64876, -.62426, -.6007, -.57805,
     & -.55625, -.53529, -.51513, -.49573, -.47708, -.45913, -.44186,
     & -.42525, -.40926, -.39389, -.37909, -.36485, -.35116, -.33798,
     & -.32529, -.31309, -.30134, -.29004, -.27916, -.26869, -.25862,
     & -.24892, -.23959, -.23061, -.22196, -.21364, -.20564, -.19793,     290   
     & -.19051, -.18337, -.17649, -.16987, -.1635, -.15737, -.15147,
     & -.14579, -.14032, -.13505, -.12998, -.1251, -.1204, -.11588,
     & -.11153, -.10733, -.1033, -.099412, -.095671, -.09207,
     & -.088603, -.085264, -.08205, -.078956, -.075976, -.073108,
     & -.070346, -.067687, -.065126, -.062661, -.060288, -.058003,
     & -.055803, -.053684, -.051645, -.049681, -.04779, -.045969,
     & -.044216, -.042528, -.040903, -.039339, -.037832, -.036382,
     & -.034985, -.03364, -.032346, -.031099, -.029899, -.028743,
     & -.027631, -.026559, -.025528, -.024535, -.023579, -.022658,
     & -.021772, -.020919, -.020097, -.019307, -.018545, -.017812,        300   
     & -.017106, -.016427, -.015773, -.015143 /
 
      data spams(1:2) / .97069 , -.24034 /
      data ebounds(1:2) / -2.22457 , -2.22457  /
      data phi_tails(:,1) / .21114 ,  .2316 /
      data phi_tails(:,2) / .021156 , .23113 /
      data v_tails(:,1) /  -18.019 , .72772 /
      data v_tails(:,2) /  -656.41 , .71085 /
      data v_couls(1:2) / 0., 0. /
                                                                          310   
      SAVE
 
!     skip to setup, printing, or calculation code.
 
      iretur = 0
      numout = 0
 
!   get the case number which points to the tables above
 
      icase = l/2 + 1                                                     320   
      undef = fintrn(1)
 
      if ( ireque - 2 ) 100, 200, 300
 
!     setup:  check for error in l (only likely error).
 
 100  if ( l /= 0  .and.  l /= 2 ) then
         write ( iunit, 103 ) alias, l
 103     format ( '0****** error in linkule ', a8 /
     1      ' ******  L =', i5, '   L must = 0 or 2' )                    330   
         numout = 2
cib      rewind iunit
         iretur = -1
         return
      endif
 
!     set jp, nodes, v, r, a, e, j, jsp, jst, spam
 
      jp = l + 1
      intger(12) = 0                                                      340   
!     the potential depth is always set to 1.0
      flt(62) = 1.
      if ( flt(43) .eq. undef ) flt(43) = 1.
      if ( flt(1) .eq. undef ) flt(1) = 0.4
      if ( flt(12) .eq. undef ) flt(12) = ebounds(icase)
      jblock(1) = 2
      jblock(8) = 1
      jblock(9) = 1
!     the default spec. amplitudes
      if ( flt(53) .eq. undef ) flt(53) = spams(icase)                    350   
!     "spamp" or "spamt", as appropriate, is always set to "spam".
      i = iswtch(11)
      if ( i .eq. 1  .or.  i .eq. 2 )  flt(53+i) = flt(53)
      return
 
 
!     ireque = 2:  print stuff about the linkule.
 
 200  if ( ipotyp == 6 ) then
                                                                          360   
!     printout for wave function calculation.
 
         write ( iunit, 203 ) sd(icase), spams(icase)**2
 203     format ( 10x, 'Argonne v18 deuteron wave function,',
     1      5x, a1, '-state probability =', f7.4 )
         numout = 1
 
      else
 
!     printout for potential calculation                                  370   
         write ( iunit, 253 ) flt(62)
 253     format ( 10x, 'Effective potential derived from Argonne v18',
     1      ' deuteron wavefunction,' /
     2      5x, 'multiplier =', f9.6 )
         numout = 1
 
      endif
 
!     printout for both wf and potential calculations.
                                                                          380   
      write (iunit, 273 ) flt(53)
 273  format ( 10x, 'spectroscopic amplitude ("SPAM") =', f8.5 )
      numout = numout + 1
cib      rewind iunit
      return
 
 
!     set up for the main calculation loop.
 
 300  am = flt(24)                                                        390   
      e = flt(12)
      rt4pi = cnstnt(2)
      hbarc = cnstnt(6)
      alsq = -2.*am*e/hbarc**2
      alpha = dsqrt( alsq )
 
!     the "param1" keyword gives the radial compression factor
!     (2 means factor-of-2 compression, i.e. small deuteron)
 
      wfsw = ipotyp .eq. 6                                                400   
 
      rt = rstart
      rend = rstart + (numpts-1)*stepsz
      iend = ( min( rend, grid_max ) - rstart )/stepsz
      hi = 1./grid_step
 
!     loop through radii, up to the highest tabulated value.
 
      do  i = 1, iend
                                                                          410   
c
c three-point interpolation for fp correlations
c
         frac = rt*hi
         j = frac+0.5
         frac = frac-j
         fras = frac**2
         if (j == 0) then
            j=1
            x0=2.*frac-fras                                               420   
            xpl=-.5*(frac-fras)
            xmn=1.-1.5*frac+.5*fras
         else
            x0=1.-fras
            xpl=.5*(fras+frac)
            xmn=.5*(fras-frac)
         end if
 
         if (wfsw) then
            phi = x0*phis(j,icase) + xpl*phis(j+1,icase)                  430   
     &         + xmn*phis(j-1,icase)
            array1(i) = rt*phi
         endif
 
         v = x0*vs(j,icase)+xpl*vs(j+1,icase)+xmn*vs(j-1,icase)
         array2(i) = -v
 
         rt = rt + stepsz
 
      enddo                                                               440   
 
!   add the tails
 
      ist = iend+1
      do i = ist, numpts
 
         if (wfsw) then
            x = rt*phi_tails(2,icase)
            xi = 1/x
            f = exp(-x)*xi                                                450   
            if ( l == 1 )  f = ( 1 + xi ) * f
            if ( l == 2 )  f = ( 1 + 3*xi + 3*xi**2 ) * f
            array1(i) = rt * phi_tails(1,icase)*f
         endif
 
         x = rt*v_tails(2,icase)
         v = v_tails(1,icase) * exp(-x)/x
         v = v + v_couls(icase)/rt
         array2(i) = -v
                                                                          460   
         rt = rt + stepsz
 
      enddo
 
cib      rewind iunit
      return
      end
