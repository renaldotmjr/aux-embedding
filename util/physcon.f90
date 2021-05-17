      SUBROUTINE PHYSCON
      use global_variables
!
!     Purpose: Initialization of physical constants.
!
!     Lit.: CODATA constants from http://www.codata.org (1986)
!
!     History: - From deMon Physcon subroutine
!            : - Modification: Atomic names with lengh four.
!
!     ******************************************************************
!
!     List of variables
!
!     COVRAD  : Atomic covalent radii [Angstrom] -> [Bohr].
!     ELSYM   : Element symbols.
!     EMASS   : Electron mass [kg].
!
!     ------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER LM,Z
!
!     ------------------------------------------------------------------
!
!     *** Definition of standard covalent radii [Angstrom] ***
!     *** Lit.: R.T. Sanderson, Inorganic Chemistry, Reinhold 1967 ***
!
      COVRAD(0)  = 0.00
      COVRAD(1)  = 0.32
      COVRAD(2)  = 0.93
      COVRAD(3)  = 1.23
      COVRAD(4)  = 0.90
      COVRAD(5)  = 0.82
      COVRAD(6)  = 0.77
      COVRAD(7)  = 0.75
      COVRAD(8)  = 0.73
      COVRAD(9)  = 0.72
      COVRAD(10) = 0.71
      COVRAD(11) = 1.54
      COVRAD(12) = 1.36
      COVRAD(13) = 1.18
      COVRAD(14) = 1.11
      COVRAD(15) = 1.06
      COVRAD(16) = 1.02
      COVRAD(17) = 0.99
      COVRAD(18) = 0.98
      COVRAD(19) = 2.03
      COVRAD(20) = 1.74
      COVRAD(21) = 1.44
      COVRAD(22) = 1.32
      COVRAD(23) = 1.22
      COVRAD(24) = 1.18
      COVRAD(25) = 1.17
      COVRAD(26) = 1.17
      COVRAD(27) = 1.16
      COVRAD(28) = 1.15
      COVRAD(29) = 1.17
      COVRAD(30) = 1.25
      COVRAD(31) = 1.26
      COVRAD(32) = 1.22
      COVRAD(33) = 1.20
      COVRAD(34) = 1.16
      COVRAD(35) = 1.14
      COVRAD(36) = 1.12
      COVRAD(37) = 2.16
      COVRAD(38) = 1.91
      COVRAD(39) = 1.62
      COVRAD(40) = 1.45
      COVRAD(41) = 1.34
      COVRAD(42) = 1.30
      COVRAD(43) = 1.27
      COVRAD(44) = 1.25
      COVRAD(45) = 1.25
      COVRAD(46) = 1.28
      COVRAD(47) = 1.34
      COVRAD(48) = 1.48
      COVRAD(49) = 1.44
      COVRAD(50) = 1.41
      COVRAD(51) = 1.40
      COVRAD(52) = 1.36
      COVRAD(53) = 1.33
      COVRAD(54) = 1.31
      COVRAD(55) = 2.35
      COVRAD(56) = 1.98
      COVRAD(57) = 1.69
      COVRAD(58) = 1.65
      COVRAD(59) = 1.65
      COVRAD(60) = 1.64
      COVRAD(61) = 1.63
      COVRAD(62) = 1.62
      COVRAD(63) = 1.85
      COVRAD(64) = 1.61
      COVRAD(65) = 1.59
      COVRAD(66) = 1.59
      COVRAD(67) = 1.58
      COVRAD(68) = 1.57
      COVRAD(69) = 1.56
      COVRAD(70) = 1.74
      COVRAD(71) = 1.56
      COVRAD(72) = 1.44
      COVRAD(73) = 1.34
      COVRAD(74) = 1.30
      COVRAD(75) = 1.28
      COVRAD(76) = 1.26
      COVRAD(77) = 1.27
      COVRAD(78) = 1.30
      COVRAD(79) = 1.34
      COVRAD(80) = 1.49
      COVRAD(81) = 1.48
      COVRAD(82) = 1.47
      COVRAD(83) = 1.46
      COVRAD(84) = 1.46
      COVRAD(85) = 1.45
      COVRAD(86) = 1.90
      COVRAD(90) = 1.65
      COVRAD(92) = 1.42
      COVRAD(93) = 1.34
      COVRAD(94) = 1.55
      COVRAD(95) = 1.89
      COVRAD(96) = 2.00
      COVRAD(97) = 2.00
      COVRAD(98) = 2.00
      COVRAD(99) = 2.00
      COVRAD(100) = 2.00
      COVRAD(101) = 2.00
      COVRAD(102) = 2.00
      COVRAD(103) = 2.00
!
!     *** Definition of element symbols ***
!
!     MODIFICATION: Definition of atomic names with lengh four.
      ELSYM(0)  = 'X   '
      ELSYM(1)  = 'H   '
      ELSYM(2)  = 'He  '
      ELSYM(3)  = 'Li  '
      ELSYM(4)  = 'Be  '
      ELSYM(5)  = 'B   '
      ELSYM(6)  = 'C   '
      ELSYM(7)  = 'N   '
      ELSYM(8)  = 'O   '
      ELSYM(9)  = 'F   '
      ELSYM(10) = 'Ne  '
      ELSYM(11) = 'Na  '
      ELSYM(12) = 'Mg  '
      ELSYM(13) = 'Al  '
      ELSYM(14) = 'Si  '
      ELSYM(15) = 'P   '
      ELSYM(16) = 'S   '
      ELSYM(17) = 'Cl  '
      ELSYM(18) = 'Ar  '
      ELSYM(19) = 'K   '
      ELSYM(20) = 'Ca  '
      ELSYM(21) = 'Sc  '
      ELSYM(22) = 'Ti  '
      ELSYM(23) = 'V   '
      ELSYM(24) = 'Cr  '
      ELSYM(25) = 'Mn  '
      ELSYM(26) = 'Fe  '
      ELSYM(27) = 'Co  '
      ELSYM(28) = 'Ni  '
      ELSYM(29) = 'Cu  '
      ELSYM(30) = 'Zn  '
      ELSYM(31) = 'Ga  '
      ELSYM(32) = 'Ge  '
      ELSYM(33) = 'As  '
      ELSYM(34) = 'Se  '
      ELSYM(35) = 'Br  '
      ELSYM(36) = 'Kr  '
      ELSYM(37) = 'Rb  '
      ELSYM(38) = 'Sr  '
      ELSYM(39) = 'Y   '
      ELSYM(40) = 'Zr  '
      ELSYM(41) = 'Nb  '
      ELSYM(42) = 'Mo  '
      ELSYM(43) = 'Tc  '
      ELSYM(44) = 'Ru  '
      ELSYM(45) = 'Rh  '
      ELSYM(46) = 'Pd  '
      ELSYM(47) = 'Ag  '
      ELSYM(48) = 'Cd  '
      ELSYM(49) = 'In  '
      ELSYM(50) = 'Sn  '
      ELSYM(51) = 'Sb  '
      ELSYM(52) = 'Te  '
      ELSYM(53) = 'I   '
      ELSYM(54) = 'Xe  '
      ELSYM(55) = 'Cs  '
      ELSYM(56) = 'Ba  '
      ELSYM(57) = 'La  '
      ELSYM(58) = 'Ce  '
      ELSYM(59) = 'Pr  '
      ELSYM(60) = 'Nd  '
      ELSYM(61) = 'Pm  '
      ELSYM(62) = 'Sm  '
      ELSYM(63) = 'Eu  '
      ELSYM(64) = 'Gd  '
      ELSYM(65) = 'Tb  '
      ELSYM(66) = 'Dy  '
      ELSYM(67) = 'Ho  '
      ELSYM(68) = 'Er  '
      ELSYM(69) = 'Tm  '
      ELSYM(70) = 'Yb  '
      ELSYM(71) = 'Lu  '
      ELSYM(72) = 'Hf  '
      ELSYM(73) = 'Ta  '
      ELSYM(74) = 'W   '
      ELSYM(75) = 'Re  '
      ELSYM(76) = 'Os  '
      ELSYM(77) = 'Ir  '
      ELSYM(78) = 'Pt  '
      ELSYM(79) = 'Au  '
      ELSYM(80) = 'Hg  '
      ELSYM(81) = 'Tl  '
      ELSYM(82) = 'Pb  '
      ELSYM(83) = 'Bi  '
      ELSYM(84) = 'Po  '
      ELSYM(85) = 'At  '
      ELSYM(86) = 'Rn  '
      ELSYM(87) = 'Fr  '
      ELSYM(88) = 'Ra  '
      ELSYM(89) = 'Ac  '
      ELSYM(90) = 'Th  '
      ELSYM(91) = 'Pa  '
      ELSYM(92) = 'U   '
      ELSYM(93) = 'Np  '
      ELSYM(94) = 'Pu  '
      ELSYM(95) = 'Am  '
      ELSYM(96) = 'Cm  '
      ELSYM(97) = 'Bk  '
      ELSYM(98) = 'Cf  '
      ELSYM(99) = 'Es  '
      ELSYM(100) = 'Fm  '
      ELSYM(101) = 'Md  '
      ELSYM(102) = 'No  '
      ELSYM(103) = 'Lr  '
!
!     *** Definition of standard isotopic masses ***
!     *** Lit.: CRC Handbook of Chemistry and Physics, 1989 ***
!
      STDMATOM(0)  =  0.000000
      STDMATOM(1)  =  1.007940
      STDMATOM(2)  =  4.002602
      STDMATOM(3)  =  6.941000
      STDMATOM(4)  =  9.012182
      STDMATOM(5)  = 10.811000
      STDMATOM(6)  = 12.011000
      STDMATOM(7)  = 14.006740
      STDMATOM(8)  = 15.999400
      STDMATOM(9)  = 18.998400
      STDMATOM(10) = 20.179700
      STDMATOM(11) = 22.989768
      STDMATOM(12) = 24.305000
      STDMATOM(13) = 26.981539
      STDMATOM(14) = 28.085500
      STDMATOM(15) = 30.973762
      STDMATOM(16) = 32.066000
      STDMATOM(17) = 35.452700
      STDMATOM(18) = 39.948000
      STDMATOM(19) = 39.098300
      STDMATOM(20) = 40.078000
      STDMATOM(21) = 44.955910
      STDMATOM(22) = 47.880000
      STDMATOM(23) = 50.941500
      STDMATOM(24) = 51.996100
      STDMATOM(25) = 54.938050
      STDMATOM(26) = 55.847000
      STDMATOM(27) = 58.933200
      STDMATOM(28) = 58.693400
      STDMATOM(29) = 63.546000
      STDMATOM(30) = 65.390000
      STDMATOM(31) = 69.723000
      STDMATOM(32) = 72.610000
      STDMATOM(33) = 74.921590
      STDMATOM(34) = 78.960000
      STDMATOM(35) = 79.904000
      STDMATOM(36) = 83.800000
      STDMATOM(37) = 85.467800
      STDMATOM(38) = 87.620000
      STDMATOM(39) = 88.905850
      STDMATOM(40) = 91.224000
      STDMATOM(41) = 92.906380
      STDMATOM(42) = 95.940000
      STDMATOM(43) = 98.000000
      STDMATOM(44) = 101.070000
      STDMATOM(45) = 102.905500
      STDMATOM(46) = 106.420000
      STDMATOM(47) = 107.868200
      STDMATOM(48) = 112.411000
      STDMATOM(49) = 114.820000
      STDMATOM(50) = 118.710000
      STDMATOM(51) = 121.757000
      STDMATOM(52) = 127.600000
      STDMATOM(53) = 126.904470
      STDMATOM(54) = 131.290000
      STDMATOM(55) = 132.905430
      STDMATOM(56) = 137.327000
      STDMATOM(57) = 138.905500
      STDMATOM(58) = 140.115000
      STDMATOM(59) = 140.907650
      STDMATOM(60) = 144.240000
      STDMATOM(61) = 145.000000
      STDMATOM(62) = 150.360000
      STDMATOM(63) = 151.965000
      STDMATOM(64) = 157.250000
      STDMATOM(65) = 158.925340
      STDMATOM(66) = 162.500000
      STDMATOM(67) = 164.930320
      STDMATOM(68) = 167.260000
      STDMATOM(69) = 168.934210
      STDMATOM(70) = 173.040000
      STDMATOM(71) = 174.967000
      STDMATOM(72) = 178.490000
      STDMATOM(73) = 180.947900
      STDMATOM(74) = 183.850000
      STDMATOM(75) = 186.207000
      STDMATOM(76) = 190.200000
      STDMATOM(77) = 192.220000
      STDMATOM(78) = 195.080000
      STDMATOM(79) = 196.966540
      STDMATOM(80) = 200.590000
      STDMATOM(81) = 204.383300
      STDMATOM(82) = 207.200000
      STDMATOM(83) = 208.980370
      STDMATOM(84) = 209.000000
      STDMATOM(85) = 210.000000
      STDMATOM(86) = 222.000000
      STDMATOM(87) = 223.000000
      STDMATOM(88) = 226.000000
      STDMATOM(89) = 227.000000
      STDMATOM(90) = 232.038100
      STDMATOM(91) = 231.035880
      STDMATOM(92) = 238.028900
      STDMATOM(93) = 237.000000
      STDMATOM(94) = 244.000000
      STDMATOM(95) = 243.000000
      STDMATOM(96) = 247.000000
      STDMATOM(97) = 247.000000
      STDMATOM(98) = 251.000000
      STDMATOM(99) = 252.000000
      STDMATOM(100) = 257.000000
      STDMATOM(101) = 258.000000
      STDMATOM(102) = 259.000000
      STDMATOM(103) = 262.000000
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE PHYSCON ***
!
      END SUBROUTINE PHYSCON
