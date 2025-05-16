       ! sse2.f
      subroutine sse_evolve(mass, z, tphysf, 
     &                      neta_in, bwind_in, hewind_in, sigma_in,
     &                      ifflag_in, wdflag_in, bhflag_in, nsflag_in,
     &                      mxns_in, pts1_in, pts2_in, pts3_in,
     &                      scm_out, phase_out, time_out)
Cf2py intent(in) mass, z, tphysf
Cf2py neta_in, bwind_in, hewind_in, sigma_in
Cf2py intent(in) ifflag_in, wdflag_in, bhflag_in, nsflag_in
Cf2py intent(in) mxns_in, pts1_in, pts2_in, pts3_in
Cf2py intent(out) type_out, Mo_out, Mt_out, L_out, R_out, T_out, Mc_out
Cf2py intent(out) Menv_out, epoch_out, spin_out, rc_out, Renv_out
Cf2py intent(out) scm_out, phase_out, time_out

c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
c      INCLUDE 'const_bse.h'
*
      integer kw,j,k
*
      real*8 mass,mt,z,zpars(20)
      real*8 neta_in, bwind_in, hewind_in, sigma_in, mxns_in
      real*8 type_out, Mo_out, Mt_out, L_out, R_out, T_out, Mc_out
      real*8 rc_out, Renv_out, Menv_out, epoch_out, spin_out
      real*8 epoch,tms,tphys,tphysf,dtp
      real*8 r,lum,ospin
      real*8 mc,rc,menv,renv
      character*50 text1,text2,text3
      character*30 label(16)

      data label /' Low Mass MS Star ',' Main sequence Star ',
     &            ' Hertzsprung Gap ',' Giant Branch ',
     &            ' Core Helium Burning ',
     &            ' First AGB ',' Second AGB ',
     &            ' Naked Helium MS ',' Naked Helium HG ',
     &            ' Naked Helium GB ',' Helium WD ',
     &            ' Carbon/Oxygen WD ',' Oxygen/Neon WD ',
     &            ' Neutron Star ',' Black Hole ',
     &            ' Massless Supernova '/

      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag,bhflag
      INTEGER ceflag_in,tflag_in,ifflag_in,nsflag_in,wdflag_in,bhflag_in
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 sigma,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      REAL*8 pts1_in,pts2_in,pts3_in
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      REAL scm_out(50000,13)
      real phase_out(20), time_out(20)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp
*
*
************************************************************************
* Input:
*
* mass is in solar units.
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* tphysf is the maximum evolution time in Myr.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally). 
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
*
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used in the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05)
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* If you enter a negative mass then parameters for an evolved star are
* required in the order of:
* initial mass, current mass, type, current time & epoch,
* otherwise the star will start on the ZAMS.
*
      neta = neta_in
      bwind = bwind_in
      hewind = hewind_in
      sigma = sigma_in
      ifflag = ifflag_in
      wdflag = wdflag_in
      bhflag = bhflag_in
      nsflag = nsflag_in
      mxns = mxns_in
      pts1 = pts1_in
      pts2 = pts2_in
      pts3 = pts3_in

*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
      if(idum.gt.0) idum = -idum
*
      if(mass.gt.0.0)then
*
* Initialize the star to begin on the ZAMS.
*
         mt = mass
         kw = 1
         tphys = 0.d0
         epoch = 0.d0
      else
*
* Obtain parameters of evolved star.
*
         READ(22,*)mass,mt,kw,tphys,epoch
      endif
      CLOSE(22)
c      WRITE(*,*)
*
* Set the initial spin of the star. If ospin is less than or equal to 
* zero at time zero then evolv1 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
*
      ospin = 0.d0
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the scm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0.d0
* 
      CALL evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &            epoch,tms,tphys,tphysf,dtp,z,zpars)
*
************************************************************************
* Output:
*
      j = 0
      if(scm(1,1).lt.0.0) goto 50
*
* The scm array stores the stellar parameters at the specified output 
* times. The parameters are (in order of storage):
*    
*    Time, stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, epoch and spin.
*
c      OPEN(23,file='evolve2.dat',status='unknown')
c      text1 = ' Tev(Myr)    type      Mo        Mt      log10(L) '
c      text2 = ' log10(R) log10(Teff)  Mc        Menv     ' 
c      text3 = ' epoch      spin' 
c      WRITE(23,'(a,a,a)')text1,text2,text3
 30   j = j + 1
      if(scm(j,1).lt.0.0)then
         scm_out(1:j-1, 1:13) = scm(1:j-1,1:13)
         scm(j-1,1) = scm(j,1)
         j = j - 1
      endif
      WRITE(23,99)(scm(j,k),k=1,8),scm(j,10),scm(j,12),scm(j,13)
      if(scm(j,1).ge.0.0) goto 30
      CLOSE(23)
 99   FORMAT(8f10.4,1p,e12.4,0p,f12.4,1p,e12.4)
*
* The spp array acts as a log, storing the time and mass at each change
* of evolution stage.
*
      j = 0
 50   j = j + 1
      if(spp(j,1).lt.0.0) goto 60
      kw = INT(spp(j,2))
      phase_out(j) = kw
      time_out(j) = spp(j,1)

c      WRITE(*,100)label(kw+1),spp(j,1),spp(j,3)
      goto 50
 60   continue
c 100  format(a30,' Time [Myr] ',f10.1,' Mass [M_sun]',f7.3)
c      WRITE(*,*)
*
************************************************************************
*
      END
