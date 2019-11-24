Takeda Model 1 Benchmark XS
 2 4
 1.0E+07 0.68256
!
!Data here is derived from the OECD/NEA report NEA/NEACRP/L(1990)330.
!Titled "Final Report on the 3-D Neutron Transport Benchmarks" by T. Takeda
!and H. Ikeda. 
!
!These cross sections refer specifically to Model 1 which is a small LWR.
!

!Cross sections for the Core
! absorption(g) nu*fission(g) kappa*fission(g) chi(g)
! scattering matrix
XSMACRO CORE 0
  8.52709E-03 9.09319E-03 9.09319E-03 1.0
  1.58196E-01 2.90183E-01 2.90183E-01 0.0 
  1.92423E-01 0.00000E+00
  2.28253E-02 8.80439E-01
  
!Cross sections for the Reflector
XSMACRO REFL 0
  4.16392E-04 0.00000E+00 0.00000E+00 0.0
  2.02999E-02 0.00000E+00 0.00000E+00 0.0
  1.93446E-01 0.00000E+00
  5.65042E-02 1.62452E+00
  
!Cross sections for the Control Rod
XSMACRO CROD 0
  1.74439E-02 0.00000E+00 0.00000E+00 0.0
  1.82224E-01 0.00000E+00 0.00000E+00 0.0
  6.77241E-02 0.00000E+00
  6.45461E-05 3.52358E-02

!Cross sections for the Empty (Void)
XSMACRO VOID 0
  4.65132E-05 0.00000E+00 0.00000E+00 0.0
  1.32890E-03 0.00000E+00 0.00000E+00 0.0
  1.27700E-02 0.00000E+00
  2.40997E-05 1.07387E-02
  
