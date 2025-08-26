****
*
*     Cas neutre
*
*     1. spfa est un facteur qui permet de corriger la constante 
*     de normalisation asymptotique
*
*     2. la fonction d'onde libre est multipliee par un facteur qui
*     la fait correspondre a u-tilde a basse energie
*
*     3. Un calcul de la fonction d'onde libre a energie 0 est
*     inclus mais sa normalisation est arbitraire
*
*****
* nuclear potential
      function vnuc(x,l,i,j,ityp)
      implicit real*8(a-h,o-z)
      common/parapot/para(20),npara
      npara=4
      z=j*(j+2)-i*(i+2)-4*l*(l+1)
      if(ityp.eq.2)then
c  Woods-Saxon potential
c     p2=para(2)/para(4)
         e=exp(-(x-para(3))/para(4))
c     vnuc=(-para(1)-p2*z/(8*x)*e/(1+e))/(1+e)
c         vnuc=(-para(1)-para(2)*z/8)*e/(1+e)
         vnuc=(-para(1)-para(2)/para(4)*z/8/x/(1+e))*e/(1+e)

      else
c  Gaussian potential
c  Changed from original Racam to fit the definition of the EFT potential
c  which are of Gaussian shape.
*      p2=para(2)/para(3)**2
*      vnuc=exp(-(x/para(3))**2)*(-para(1)-p2*z)
         e=(x-para(3))/para(4)
         e=exp(-e**2/2)     !Gaussian form factor
         vnuc=-(para(1)+para(2)*(x-para(3))**2)*e
      endif
      return
      end

* main program
c My version of Racam from Daniel Baye.
c The program has been changed in the definition of the potential,
c to include the definition of the EFT Gaussian potential.

c New self-standing version of the code
      program RACAM
      implicit real*8(a-h,o-z)
      parameter(nmax=125000)
      character hpar(2)*2,noy(0:20)*2
      character*80 nom
      character*85 fichfo,fichsca
      character*2 typ
      dimension potf(nmax),wff(nmax),poti(nmax),wfi(nmax),nlam(0:4)
     1 ,wf0(nmax)
      common/parapot/para(20),npara
c  fsc=inverse of the fine-structure constant
c  amn=nucleon mass (in MeV)
c  hc=hbar*c  (in MeV*fm)
      data fsc,amn,hc/137.0359895d0,931.4943228d0,197.32705359d0/
      data hpar/'/2','  '/,nlam/1,3,15,105,945/ !840???
      data pi/3.14159265359d0/ ,eps/1d-6/
      data noy/' n',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     1 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca'/
      c=2.99792458d23
      e2=hc/fsc
      fsc=hc/e2
      hm=hc*hc/2/amn
      print*,' hbar^2/2m',hm
*
*      print*,' data file?'
*      read1007,nom
*      print1007,nom
cc      read(*,*)nom
      read(*,'(a80)')nom
*      open(unit=1,file=nom)
c     print*,' file for figure ?'
c     read1007,nom
c     print1007,nom
c     open(unit=2,file=nom)
*      read (1,2000)a1,a2,kz1,kz2
      read (*,*)a1,a2,kz1,kz2
cc      print*,'A1=',a1,'  A2=',a2,'  Z1=',kz1,'  Z2=',kz2
      write(*,1010)a1,a2,kz1,kz2
*      read (1,2001)lam,i1,i2,if,li,ji,lf,jf,n0
      read (*,*)lam,i1,i2,if,li,ji,lf,jf,n0
cc      print*,' lam=', lam,' I1=',i1,' I2=',i2,' I=',if
cc     1,' li=',li,' ji=',ji,' lf=',lf,' jf=',jf,' nr=',n0
      write(*,1011)lam,i1,i2,if,li,ji,lf,jf,n0
*     read (1,2002)rmax,n
      read (*,*)rmax,n
cc      print*,'  rmax=',rmax,'  N=',n
      write(*,1012)rmax,n
*     read (1,2004)typ
      read (*,*)typ
      if(typ.eq.'G ')then
      print*,' Gaussian potential'
      ityp=1
c Original version of Racam
*      read (*,*)v0i,vsoi,rni,v0f,vsof,rnf,rc
*      print*,' V0i=',v0i,' Vsoi=',vsoi,' rni=',rni
*      print*,' V0f=',v0f,' Vsof=',vsof,' rnf=',rnf
c New version that suits the definition of the EFT Gaussian potentials
      read (*,*)v0i,vsoi,rni,ai,v0f,vsof,rnf,af,rc
cc      print*,' V0i=',v0i,' Vsoi=',vsoi,' rni=',rni,' ai=',ai
cc      print*,' V0f=',v0f,' Vsof=',vsof,' rnf=',rnf,' af=',af
      elseif(typ.eq.'ws'.or.typ.eq.'WS'
     a     .or.typ.eq.'Ws'.or.typ.eq.'wS')then
      print*,' Woods-Saxon potential'
      ityp=2
*      read (1,2003)v0i,vsoi,rni,ai,v0f,vsof,rnf,af,rc
      read (*,*)v0i,vsoi,rni,ai,v0f,vsof,rnf,af,rc
cc      print*,' V0i=',v0i,' Vsoi=',vsoi,' rni=',rni,' ai=',ai
cc      print*,' V0f=',v0f,' Vsof=',vsof,' rnf=',rnf,' af=',af
      else
         print*,' Undefined potential'
         STOP
      endif
      close(1)
      if(ji.lt.0.and.vsoi.ne.0)then
      print*
      print*,' ji > 0 (Vsoi ne 0)'
      stop
      endif
c         lam=multipole order
c         n=number of discretization points (typically  100 ~ 1000)
c         h=step in fm  (typically 0.02 ~ 0.1)
      h=rmax/n
cc      print*,'lam=',lam,'  n=',n,'  h=',h
      if(lam.lt.0.or.lam.gt.4.or.n.le.0.or.n.gt.nmax.or.h.le.0)then
      print*,'Input data are incorrect !!!'
      stop
      endif
cc      print*,' Last integration point=',n*h,' fm'
      lam2=2*lam+1

c         a1,a2=masses (real)
c         kz1,kz2=charges (integer)
c         i1,i2=spins (exact if integer, twice if half-integer)
c         rc=Coulomb radius (in fm)
      n1=a1+0.5d0
      n2=a2+0.5d0
      lz=2-mod(n1+n2,2)
      lz1=2-mod(n1,2)
      lz2=2-mod(n2,2)
      print1001,noy(kz1),n1,noy(kz2),n2,noy(kz1+kz2),n1+n2
      print1004,a1,a2,kz1,kz2,i1,hpar(lz1),i2,hpar(lz2),rc
      i1i2=(lz1*i1+1)*(lz2*i2+1)
      ama=a1+a2
      rmu=a1*a2/ama
      rm=rmu/hm
      print*,' hbar^2/2mu',1/rm
      ze=kz1*kz2*e2
      eta0=kz1*kz2*sqrt(rmu*amn/2)/fsc
      if(lam.ne.0)fac0=8*pi*lam2*(lam+1)*((kz1*(a2/ama)**lam+
     1 kz2*(-a1/ama)**lam)/nlam(lam))**2/(lam*i1i2)

c         lf=orbital momentum of the final state (integer)
c         if=channel spin (arising from the coupling of i1 and i2)
c         jf=total spin
c         n0=number of nodes in the final wave function
c         (para(i),i=1,npar)= parameters of the final potential
c         spfa=spectroscopic factor
    7 continue
      spfa=1
      npara=4
      para(1)=v0f
      para(2)=vsof
      para(3)=rnf
      para(4)=af
      if(lf.lt.0)stop
      jf=jf*lz
      if=if*lz
      if(jf.lt.abs(2*lf-if).or.jf.gt.2*lf+if.or.n0.lt.0)then
      print*,'Input data are incorrect !!!'
      stop
      endif
      s2=0 
      e=0
c computing the final potential (potf)
      umax=0
      do2j=1,n
      x=j*h
      vn=vnuc(x,lf,if,jf,ityp)
      vc=ze/x
      if(x.le.rc)vc=ze*(3-(x/rc)**2)/2/rc
      vce=lf*(lf+1)/x/x
      potf(j)=(vn+vc)*rm+vce
c     if(lam.ne.0.or.j.le.2.or.potf(j-1).le.max(potf(j-2),potf(j)))goto2
c     umax=potf(j-1)
c     imax=j-1
*      write(*,*)x,potf(j)
    2 continue
      print*,' '
      print1002,lf,if/lz,hpar(lz),jf/lz,hpar(lz),n0,spfa,
     1 (para(i),i=1,npara)
c     print*,'Final potential (in MeV) given by step of',10*h,' fm'
c     print1003,(potf(j)/rm,j=10,n,10)
      umax=0
c     if(lam.eq.0)then
c     do j=1,imax
c     potf(j)=potf(j)-umax
c     enddo
c     do j=imax+1,nmax
c     potf(j)=0
c     enddo
c     endif
c computing the final energy (ebound) and wave function (wff)
cc      print*,umax
c     print1003,(potf(j),j=1,n,100)
      call num1l(n,h,e,s2,potf,wff,n0,eps)
c     print1000,(wff(j),j=10,3000,100)
c1000 format(1h ,5e13.4)
cc      print*,e,1/rm
      ebound=(e+umax)/rm
      if(n0.eq.-1)then
      print*,'no bound state with this potential'
      goto7
      endif
      print*,'Bound-state energy=',ebound,' MeV'
      if(lam.eq.0)goto7
c         li=orbital momentum of the initial state (integer)
c         ji=total spin
c         (para(i),i=1,npar)= parameters of the initial potential
      para(1)=v0i
      para(2)=vsoi
      para(3)=rni
      para(4)=ai
      ji=ji*lz
      if(ji.lt.abs(2*li-if).or.ji.gt.2*li+if.or.lam*2.lt.abs(ji-jf).
     1 or.lam*2.gt.ji+jf.or.lam.lt.abs(li-lf).or.lam.gt.li+lf)then
      print*,'Input data are incorrect !!!'
      stop
      endif
c computing the initial potential (poti)
      do3j=1,n
      x=j*h
      vn=vnuc(x,li,if,ji,ityp)
      vc=ze/x
      if(x.le.rc)vc=ze*(3-(x/rc)**2)/2/rc
      vce=li*(li+1)/x/x
      poti(j)=(vn+vc)*rm+vce
*      write(*,*)x,poti(j),vc*rm
    3 continue
      print*,' '
      print1006,li,ji/lz,hpar(lz),(para(i),i=1,npara)
c     print*,'Initial potential (in MeV) given by step of',10*h,' fm'
c     print1003,(poti(j)/rm,j=10,n,10)
*
      call num0(n,h,s2,poti,wf0)
*
      call clebs(2*li,2*lam,2*lf,0,0,0,c1)
      call sixj(jf,ji,2*lam,2*li,2*lf,if,s1)
      fac1=fac0*(jf+1)*(ji+1)*(2*li+1)*(c1*s1)**2
cc ajout pour le calcul 14C(n,g)15C
cc      ji1=1
cc      call sixj(jf,ji1,2*lam,2*li,2*lf,if,s11)
cc      fac11=fac0*(jf+1)*(ji1+1)*(2*li+1)*(c1*s11)**2
cc      print*,' N_E = ',fac1
cc      print*,'  Energy ','     k   ','   eta   ',' phase shift '
cc     1 ,'cross section ',' S factor '
cc      print*,'  (MeV)  ','(fm^{-1}) ','       ','   (degrees) '
cc     1 ,'    (fm^2)     ','(fm^(2l+1))'

c reading ne=number of energies
c         e0=first energy (in MeV)
c         ep=energy step (in MeV)
*    1 print*,' ne?, e0?, ep?'
      read(*,*)ne,e0,ep
      if(ne.eq.0)stop
      read(*,*)lec

      if(lec==0)then !simply plot the wf for the bound and the 1st cont. states
         fichfo=trim(nom)//'.dwf'
         eta=eta0/sqrt(e0)
         qk=sqrt(e0*rm)
         call dephase0(n,h,poti,wfi,eps,dep,li,eta,qk,rfin)
         open(unit=68,status='unknown',file=fichfo)
         write(68,*)'#  r    ','    wff      ','     wfi     '
         write(68,*)'# (fm)  ',' (fm^{-1/2}) ',' (fm^{-1/2}) '
         do ir=1,n
            write(68,1013)ir*h,wff(ir),wfi(ir)
         enddo
         close(68)
         stop
      endif

      fichsca=trim(nom)//'.sca'      
      open(unit=69,status='unknown',file=fichsca)
      if(kz1*kz2.ne.0)then
      write(69,*)'#  Energy ','     k   ','     eta   ',' phase shift '
     1 ,'cross section ',' S factor '
      write(69,*)'#  (MeV)  ','  (fm^{-1})','         ','  (degrees)  '
     1 ,'     (b)      ',' (MeV b)  '
      else
      write(69,*)'#  Energy ','     k   ','     eta   ',' phase shift '
     1 ,'cross section ',' S factor '
      write(69,*)'#  (MeV)  ','  (fm^{-1})','         ','   (degrees) '
     1 ,'    (fm^2)     ','(fm^(2l+1))'
      endif         
c     do 5 ie=1,13
c     if(ie.le.1)e0=.01
c     if(ie.gt.1)e0=(ie-1)*.05
c     do 5 ie=1,28
c     if(ie.le.9)e0=ie*.001
c     if(ie.gt.9)e0=(ie-9)*.01
c     if(ie.ge.20)e0=(ie-18)*.1
      do5 ie=1,ne
      eta=eta0/sqrt(e0)
      qk=sqrt(e0*rm)
      vit=sqrt(2*e0/(rmu*amn))*c
c computing the initial wave function (wfi) and phase shift (dep)
      call dephase0(n,h,poti,wfi,eps,dep,li,eta,qk,rfin) 
c     print1000,(wfi(j),j=10,3000,100)
*
      if(kz1*kz2.ne.0)then
      fac2=fac1*((e0-ebound)/hc)**lam2/fsc*hc/200
cc      fac2=(fac1+fac11)*((e0-ebound)/hc)**lam2/fsc*hc/200
      fac3=e0*exp(2*pi*eta)/(qk*qk)*sqrt(rmu*amn/e0/2)*2/hc
      else
      fac2=fac1*((e0-ebound)/hc)**lam2/fsc*rmu*amn/hc
cc      fac2=(fac1+fac11)*((e0-ebound)/hc)**lam2/fsc*rmu*amn/hc
      fac3=qk**(-2*li-2)
      endif
*
      do j=1,n
      wff(j)=wff(j)*spfa
      enddo
      do j=1,n
c     wfi(j)=wfi(j)
      wfi(j)=wfi(j)*sqrt(fac3)
      enddo
*
      ra=wfi(20)/wf0(20)
cc      lec=1
cc      if(lec.gt.1)then
cc      print*,'   r         wff          wfi          wf0'
cc      do j=20,300,20
cc      print2222,j*h,wff(j),wfi(j),wf0(j)*ra,wff(j)*wfi(j)*(j*h)**lam
cc      enddo
cc      endif
c     do j=10,n,10
c     write(2,2222)j*h,wfi(j)
c     enddo
c     write(2,2222)
c     do j=2,n,2
c     write(2,2222)j*h,wff(j)
c     enddo
c     write(2,2222)
c     do j=10,n,10
c     write(2,2222)j*h,wff(j)*wfi(j)*(j*h)**lam
c     enddo
c     write(2,2222)
 2222 format(1h ,f8.2,1p,5e13.4)
*
      rmax=0
      zmax=0
      elmat=0
      elmat0=0
      do6 j=1,n
      z=wff(j)*wfi(j)*(j*h)**lam
      if(z.gt.zmax)then
      rmax=j*h
      zmax=z
      endif
      elmat=elmat+z
c     if(mod(j,100).eq.10.and.j.le.3000)
c    1 print5000,j*h,wfi(j),wff(j),z,elmat*h
c5000 format(1h ,f10.4,5e15.4)
      elmat0=elmat0+wff(j)*wf0(j)*(j*h)**lam*ra
    6 continue
cc      if(lec.gt.0)then
cc      print*,' rmax = ',rmax,' zmax = ',zmax
cc      print*,elmat*h,elmat0*h
cc      endif
      elmat=(elmat*h)**2*fac2
*      print1005,e0,qk,eta,dep*180/pi,elmat*qk**(2*li-1),elmat/qk*vit
      if(kz1*kz2.ne.0)then
      write(69,1005)e0,qk,eta,dep*180/pi,elmat/e0/exp(2*pi*eta),elmat
      else
      write(69,1005)e0,qk,eta,dep*180/pi,elmat*qk**(2*li-1),elmat/qk*vit
      endif
      e0=e0+ep
    5 continue
      close(unit=69)
*      goto1
 1001 format(/,' Reaction: ',a2,i2,'(',a2,i2,',gamma)',a2,i2,/) 
 1002 format(' lf=',i1,' If=',i2,a2,' Jf=',i2,a2,' n0=',i1,
     1 ' spectroscopic factor=',f10.6,/,' potential parameters=',6f9.4)
 1003 format(1h ,15f8.3) 
 1004 format(/,' Masses=',2f10.5,'  charges=',2i3,'  spins=',
     1 i2,a2,2x,i2,a2,'  Rc=',f6.3,' fm')
 1005 format(f12.7,2f9.5,1p,4e13.5)
 1006 format(' li=',i1,' Ji=',i2,a2,' potential parameters=',6f8.3)
 1007 format(a80)
 1010 format('A1=',f10.5,'  A2=',f10.5,'  Z1=',i3,'  Z2=',i3)
 1011 format(' lam=',i2,' I1=',i2,' I2=',i2,' I=',i2,
     a ' li=',i2,' ji=',i2,' lf=',i2,' jf=',i2,' nr=',i2)
 1012 format('  rmax=',f8.2,'  Nr=',i6)
 1013 format(f8.2,2e13.5)
 2000 format(2(f13.10,/),2(i13,/))
cc 2001 format(9(i13,/))
 2002 format(1(f13.10,/),1(i13,/))
 2003 format(9(f13.10,/))
 2004 format(1(a2))
      end 
*num0
      subroutine num0 (n,h,s2,u,s) 
c*****integration of Schrodinger equation by the Numerov method
c*****at E=0
      implicit real*8 (a-h,o-z) 
      dimension u(1),s(1)
      h12=h*h/12 
      e=0
      s(1)=1.0d-10 
      b0=0
      aa=h12*u(1) 
      if (s2) 16,18,16
   16 b0=-s(1)*aa 
   18 b1=s(1)*(1-aa) 
      do 38 k=2,n 
      b2=12*s(k-1)-10*b1-b0 
      if (abs(b2).lt.1.e+10) go to 22 
      b2=b2*1.0d-20
      b1=b1*1.0d-20
   22 aa=h12*u(k) 
      s(k)=b2/(1-aa) 
      b0=b1 
   38 b1=b2 
      return
      end 

*NUM1L
      SUBROUTINE NUM1L(N,H,E,S2,U,S,NO,EPS) 
CC  VERSION CORRIGEE LE 21 NOV 72 
C*****INTEGRATION DE L"EQUATION DE SCHROEDINGER PAR LA METHODE DE NUMERO
C*****POUR E NEGATIF  
C*****RECHERCHE DE L"ENERGIE PROPRE PAR LA METHODE DE RAPHSON-NEWTON
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION U(1),S(1) 
      DATA RAP1,RAP2/0,0/ 
      H12=H*H/12 
C*****CONTROLE DES CONDITIONS ASYMPTOTIQUES 
      IF(E.GT.0) E=0
      DEI=0
      EPSS=.1d-10 
      IF(U(N-1).GT.EPSS) GO TO 10 
      DEI=U(N-1)-EPSS 
      DO 8 K=1,N
    8 U(K)=U(K)-DEI 
   10 U(N)=U(N-1) 
C*****CALCUL DU NOMBRE D"ETATS LIES PAR INTEGRATION A ENERGIE NULLE 
      S(1)=1.d-10 
      B0=0 
      AA=H12*U(1) 
      IF (S2) 16,18,16
   16 B0=-S(1)*AA 
   18 B1=S(1)*(1-AA) 
      DO 38 K=2,N 
      B2=12*S(K-1)-10*B1-B0 
      IF (ABS(B2).LT.1.d+10) GO TO 22 
      B2=B2*1.d-20
      B1=B1*1.d-20
   22 AA=H12*U(K) 
      S(K)=B2/(1-AA) 
      B0=B1 
   38 B1=B2 
      DO 42 K=5,N 
      N0=K  
      IF(U(K).LT.0) GO TO 44 
   42 CONTINUE
   44 NEL=0 
      DO 52 K=N0,N
      IF (S(K-1)*S(K)) 46,50,52 
   46 NEL=NEL+2 
      GO TO 52
   50 NEL=NEL+1 
   52 CONTINUE
      NEL=NEL/2 
      IF(NEL.GT.NO) GO TO 64
      IF(NEL.EQ.NO) GO TO 60
   62 NO=-1 
      RETURN
   60 RAP1=S(N-1)/S(N)
      RAP2=EXP(H*SQRT(U(N-1)-E))
      IF(RAP1.LT.RAP2) GO TO 62 
C*****CALCUL DE EMIN ET EMAX ENTRE LESQUELLES SE TROUVE L"ENERGIE PROPRE
   64 UMIN=U(1) 
      DO 70 K=2,N 
      IF(U(K).LT.UMIN) UMIN=U(K)
   70 CONTINUE
      EMIN=UMIN 
      EMAX=0 
C*****DEBUT DE LA RECHERCHE DE L"ENERGIE PROPRE DANS L"INTERVALLE MAXIMU
      TE=EMAX-EMIN
C*****REJET DE L"ENERGIE D"ESSAI E PROPOSEE SI ELLE EST A L"EXTERIEUR DE
C*****BORNES (EMIN,EMAX)
      IF((E.LT.EMIN).OR.(E.GT.EMAX)) E=EMIN+TE/2 
      E1=EMIN 
      E2=EMAX 
      J=2 
      I=1 
      GO TO 102 
C*****REDUCTION DES BORNES EMIN ET EMAX 
   90 EMIN=E1 
      EMAX=E2 
      TE=EMAX-EMIN
      J=2 
   98 I=1 
  100 E=EMIN+TE*I/J 
  102 DE=0 
  104 E=E+DE
      IF(E.GT.0) GO TO 204 
      S(N)=1.d-10 
      N1=N-1
      RAP2=EXP(H*SQRT((U(N-1)+U(N))/2-E))
      S(N1)=S(N)*RAP2 
      AA=H12*(U(N1)-E)
      B0=S(N)*(1-AA) 
      B1=S(N1)*(1-AA)
      N1=N-2
      DO 138 KAUX=1,N1
      K=N1-KAUX+1 
      B2=12*S(K+1)-10*B1-B0 
      AA=H12*(U(K)-E) 
      S(K)=B2/(1-AA) 
      B0=B1 
      B1=B2 
      IF(U(K).LT.E) GO TO 140 
  138 CONTINUE
  140 N1=K  
C*****NORMALISATION DE LA FONCTION D"ONDE A S(N1) 
      DO 146 KAUX=N1,N
      K=N-KAUX+N1 
  146 S(K)=S(K)/S(N1) 
C*****DEBUT DE L"INTEGRATION VERS L"EXTERIEUR JUSQU"A N1
      S(1)=1.d-10 
      B0=0 
      AA=H12*(U(1)-E) 
      IF(S2) 156,158,156
  156 B0=-S(1)*AA 
  158 B1=S(1)*(1-AA) 
      DO 170 K=2,N1 
      B2=12*S(K-1)-10*B1-B0 
      AA=H12*(U(K)-E) 
      S(K)=B2/(1-AA) 
      B0=B1 
  170 B1=B2 
C*****NORMALISATION DE LA FONCTION A S(N1)
      DO 174 K=1,N1 
  174 S(K)=S(K)/S(N1) 
C*****CALCUL DE LA CORECTION D"ENERGIE
      SOM=0
      DO 180 K=1,N
  180 SOM=SOM+S(K)*S(K) 
      DE=((-S(N1-1)+2-S(N1+1))/(H*H)+U(N1)-E)/SOM  
      IF(ABS(DE).GT.EPS) GO TO 104
C*****CALCUL DU NOMBRE DE NOEUDS DE L"ETAT PROPRE TROUVE
      DO 182 K=5,N
      IF(U(K).LT.E) GO TO 184 
  182 CONTINUE
  184 N0=K  
      NEL=0 
      DO 192 K=N0,N1  
      IF(S(K-1)*S(K)) 186,190,192 
  186 NEL=NEL+2 
      GO TO 192 
  190 NEL=NEL+1 
  192 CONTINUE
      NEL=NEL/2 
C*****L"ETAT PROPRE TROUVE EST-IL LE BON  
      IF(NEL-NO) 198,214,202
  198 IF(E.GT.E1) E1=E
      GO TO 204 
  202 IF(E.LT.E2) E2=E
  204 I=I+2 
      IF (I.LE.J)  GO TO 100
      J=2*J 
      IF(ABS(E1-EMIN).GT.EPS.OR.ABS(EMAX-E2).GT.EPS) GO TO 90 
      GO TO 98
C*****NORMALISATION DE LA FONCTION PROPRE 
  214 SOM=1/SQRT(SOM*H)
      DO 218 K=1,N
  218 S(K)=S(K)*SOM 
      E=E+DEI 
      RETURN
C*****DEBUT FORMATS 
 2000 FORMAT(/,35X,56HL"ETAT DEMANDE N"EST PAS LIE. RETOUR DE NUM1L AVEC
     1 NO=-1,/) 
C*****FIN FORMATS 
      END 

*DEPHASE
      SUBROUTINE DEPHASE0(MAX,H,W,Y,EPS,DELTA,L,ETA,QK,RFIN) 
c my version of dephase with a call of coufra(*,*,0,L etc)
c this helps getting a better convergence at low energy
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION FG(2),DFC(101),GC(101),DGC(101),FC(101),W(1),Y(1) 
      DD=0 
      PAS=H*EPS 
      QQ=QK*QK
      LL=L+1
      L1=L*LL 
      R=2*ETA*QK
      RMAX=MAX*H
      DO20 K=2,MAX
      RMAX=RMAX-H 
      IF(ABS(W(MAX-K+1)-(L1/RMAX+R)/RMAX).LT.EPS)GOTO20 
      IF(K.GT.3)GOTO21
      PRINT2001 
      STOP  
   20 CONTINUE
      RMAX=H*(MAX-2)  
   21 H2 = H**2 
      H212 = H2/12 
      AA=H212*(QQ-W(1)) 
      Y(1)=H**LL
      B0=0  
      B1=Y(1)*(1+AA)  
      DO22 K=2,MAX
      AA=H212*(QQ-W(K)) 
      B2=12*Y(K-1)-10*B1-B0 
      Y(K)=B2/(1+AA)  
      B0=B1 
      B1=B2 
   22 CONTINUE
      R=RMAX
      N1=INT((RMAX+1d-5)/H)+1 
      CALL COUFRA(QK*R,ETA,0,L,FC,DFC,GC,DGC) 
      FG(1)=FC(LL)
      FG(2)=GC(LL)
      DO23 N=N1,MAX 
      R=R+H 
      CALL COUFRA(QK*R,ETA,0,L,FC,DFC,GC,DGC) 
      C=FG(2)*Y(N)-GC(LL)*Y(N-1)
      D=FC(LL)*Y(N-1)-FG(1)*Y(N)
      Z=FG(1)*GC(LL)-FG(2)*FC(LL) 
      DEL=ABS(Z)/SQRT(C*C+D*D)  
      IF(ABS(1-DD/DEL).LT.PAS) GO TO 24 
      DD = DEL
      FG(1)=FC(LL)
      FG(2)=GC(LL)
   23 CONTINUE
      PRINT2000,L,QK  
   24 RFIN=R
      D=D/C 
      IF(Y(N)*(FC(LL)+GC(LL)*D).LT.0.)DEL=-DEL
      DELTA=ATAN(D) 
      DO25 K=1,MAX
   25 Y(K)=Y(K)*DEL 
      RETURN
* 2000 FORMAT(/27H PAS DE CONVERGENCE POUR L=,I3,5H   K=,F10.6/) 
 2000 FORMAT(/22H NO CONVERGENCE FOR L=,I3,5H   K=,F10.6/) 
 2001 FORMAT(/34H POTENTIAL NEVER PURELY COULOMBIAN)
      END

*COUFRA 
      SUBROUTINE COUFRA(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP)
C*** FONCTIONS COULOMBIENNES CALCULEES EN R = RHO PAR LA METHODE DES FRA
C*** CONTINUES DE STEED. MINL ET MAXL CORRESPONDENT AUX VRAIES VALEURS D
C*** VOIR BARNETT, FENG, STEED ET GOLDFARB, CPC 1974 *******************
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(MAXL),FCP(MAXL),GC(MAXL),GCP(MAXL)
      SAVE  
      DATA ACCUR,STEP/1.D-7,100.0D0/
      PACE = STEP 
      ACC = ACCUR 
      R = RHO 
      KTR = 1 
      LMAX = MAXL 
      LMIN1 = MINL+1  
      XLL1 = MINL*LMIN1 
      ETA2 = ETA**2 
      TURN = ETA+SQRT(ETA2+XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.D-6) KTR = -1  
      KTRP = KTR
      GO TO 2 
    1 R = TURN
      TF = F
      TFP = FP
      LMAX = MINL 
      KTRP = 1
    2 ETAR = ETA*R
   21 RHO2=R*R
      PL = LMAX+1 
      PMX = PL+0.5D0  
C** FRACTION CONTINUE POUR FP(MAXL)/F(MAXL) ; XL=F ; XLPRIME=FP ********
      FP = ETA/PL+PL/R
      DK = ETAR+ETAR  
      DEL = 0
      D = 0
      F = 1 
      K = (PL*PL-PL+ETAR)*(PL+PL-1) 
      IF(PL*PL+PL+ETAR.NE.0.) GO TO 3 
      R = R*1.0000001D0 
      GO TO 2 
    3 H = (PL*PL+ETA2)*(1-PL*PL)*RHO2 
      K = K+DK+PL*PL*6
      D = 1/(D*H+K) 
      DEL = DEL*(D*K-1) 
      IF(PL.LT.PMX) DEL = -R*(PL*PL+ETA2)*(PL+1)*D/PL 
      PL = PL+1 
      FP = FP+DEL 
      IF(D.LT.0) F = -F
      IF(PL.GT.20000.0D0) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP = F*FP 
      IF(LMAX.EQ.MINL) GO TO 5  
      FC(LMAX+1) = F  
      FCP(LMAX+1) = FP
C*** RECURRENCE ARRIERE POUR F ET FP ; GC,GCP UTILISES POUR STOCKAGE ***
      L = LMAX
      DO 4 LP=LMIN1,LMAX
      PL = L
      GC(L+1) = ETA/PL+PL/R 
      GCP(L+1) = SQRT(ETA2+PL*PL)/PL
      FC(L) =(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)
      FCP(L) = GC(L+1)*FC(L)-GCP(L+1)*FC(L+1) 
    4 L = L-1 
      F = FC(LMIN1) 
      FP = FCP(LMIN1) 
    5 IF(KTRP.EQ.-1) GO TO 1
C*** MEME CALCUL POUR R = TURN SI RHO.LT.TURN 
C*** P + I.Q CALCULE EN MINL , EQUATION (32)
      P = 0
      Q = R-ETA 
      PL = 0 
      AR = -(ETA2+XLL1) 
      AI = ETA
      BR = Q+Q
      BI = 2
      WI = ETA+ETA
      DR = BR/(BR*BR+BI*BI) 
      DI = -BI/(BR*BR+BI*BI)
      DP = -(AR*DI+AI*DR) 
      DQ = AR*DR-AI*DI
    6 P = P+DP
      Q = Q+DQ
      PL = PL+2 
      AR = AR+PL
      AI = AI+WI
      BI = BI+2 
      D = AR*DR-AI*DI+BR
      DI = AI*DR+AR*DI+BI 
      T = 1/(D*D+DI*DI) 
      DR = T*D
      DI = -T*DI
      H = BR*DR-BI*DI-1
      K = BI*DR+BR*DI 
      T = DP*H-DQ*K 
      DQ = DP*K+DQ*H  
      DP = T
      IF(PL.GT.46000.0D0) GO TO 11
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P = P/R 
      Q = Q/R 
C*** CALCUL DE FP,G,GP, ET NORMALISATION DE F EN L = MINL **************
      G = (FP-P*F)/Q  
      GP = P*G-Q*F
      W = 1/SQRT(FP*G-F*GP) 
      G = W*G 
      GP = W*GP 
      IF(KTR.EQ.1) GO TO 8
      F = TF
      FP = TFP
      LMAX = MAXL 
C*** CALCUL DE G(MINL) ET GP(MINL) PAR INTEGRATION RUNGE-KUTTA A PARTIR 
C***         VOIR FOX ET MAYERS(1968) PG 202
      IF(RHO.LT.0.2D0*TURN) PACE = 999.0D0
      R3=1.0D0/3.0D0  
      H = (RHO-TURN)/(PACE+1) 
      H2 = H/2
      I2 = INT(PACE+0.001D0)
      ETAH = ETA*H
      H2LL = H2*XLL1  
      S = (ETAH+H2LL/R)/R-H2
    7 RH2 = R+H2
      T = (ETAH+H2LL/RH2)/RH2-H2
      K1 = H2*GP
      M1 = S*G
      K2 = H2*(GP+M1) 
      M2 = T*(G+K1) 
      K3 = H*(GP+M2)  
      M3 = T*(G+K2) 
      M3 = M3+M3
      K4 = H2*(GP+M3) 
      RH = R+H
      S = (ETAH+H2LL/RH)/RH-H2  
      M4 = S*(G+K3) 
      G = G+(K1+K2+K2+K3+K4)*R3 
      GP = GP+(M1+M2+M2+M3+M4)*R3 
      R = RH
      I2 = I2-1 
      IF(ABS(GP).GT.1.D300) GO TO 11
      IF(I2.GE.0) GO TO 7 
      W = 1/(FP*G-F*GP) 
C*** RECURRENCE AVANT A PARTIR DE GC(MINL) ET GCP(MINL) 
C*** RENORMALISATION DE FC ET FCP POUR CHAQUE VALEUR DE L **************
    8 GC(LMIN1) = G 
      GCP(LMIN1) = GP 
      IF(LMAX.EQ.MINL) GO TO 10 
      DO 9 L=LMIN1,LMAX 
      T = GC(L+1) 
      GC(L+1) = (GC(L)*GC(L+1)-GCP(L))/GCP(L+1) 
      GCP(L+1) = GC(L)*GCP(L+1)-GC(L+1)*T 
      FC(L+1) = W*FC(L+1) 
    9 FCP(L+1) = W*FCP(L+1) 
      FC(LMIN1) = W*FC(LMIN1) 
      FCP(LMIN1) = W*FCP(LMIN1) 
      RETURN
   10 FC(LMIN1) = W*F 
      FCP(LMIN1) = W*FP 
      RETURN
   11 W = 0
      G = 0
      GP = 0 
      GO TO 8 
      END 

*CLEBS
      SUBROUTINE CLEBS (L1,L2,L3,M1,M2,M3,Q)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION FT(0:1000)
      SAVE  
      DATALMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      IS=(L1+L2+M1-M2)/2
      K3=-M3
      Q1=L3+1 
   12 Q=0  
      I1=L1 
      I2=L2 
      I3=L3 
      K1=M1 
      K2=M2 
      IF(K1+K2+K3.NE.0)RETURN 
      L=I1+I2+I3
      IF(MOD(L,2).EQ.1)GOTO8
      L=L/2 
      IF(L.LE.LMEM)GOTO6
      DO10I=LMEM,L
   10 FT(I+1)=(I+1)*FT(I) 
      LMEM=L
    6 J1=ABS(K1) 
      J2=ABS(K2) 
      J3=ABS(K3) 
      IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3)RETURN
      IF(J1+J2.EQ.0)GOTO11
      J1=I1-J1
      J2=I2-J2
      J3=I3-J3
      IF(J1)8,2,1 
    2 IF(I1.NE.0)GOTO1
      IF(J2.LT.0)GOTO8
   13 IF(J3.LT.0)GOTO8
    4 Q=SQRT(Q1/(I2+1)) 
      IS=IS+(I2-K2)/2 
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    1 IF(J2.GT.J1)GOTO3 
      IF(J2.LT.0)GOTO8
      IS=IS+L 
      J1=J2 
      K1=K2 
      K2=M1 
      I1=I2 
      I2=L1 
      IF(I1.EQ.0)GOTO13 
    3 IF(J3.GT.J1)GOTO5 
      IF(J3.LT.0)GOTO8
      IS=IS+L 
      J1=K3 
      K3=K1 
      K1=J1 
      I3=I1 
      I1=L3 
      J1=J3 
      IF(I1.EQ.0)GOTO4
    5 IF(K1.GE.0)GOTO9
      K1=-K1
      K2=-K2
      K3=-K3
      IS=IS+L 
    9 CONTINUE
      Q1=Q1*FT(L-I3)/FT(L-I1)/FT(L-I2)/FT(L+1)
      I2=(I2+K2)/2
      I3=(I3+K3)/2
      K2=I2-K2
      K3=I3-K3
      J1=J1/2 
      I1=J1+K1
      J2=I3-K2
      J3=MAX(J2,0) 
      IS=IS+I1+K2 
      X=0  
      DO7I=J3,J1
    7 X=-X+FT(I1+I)*FT(I2+I3-I)/FT(J1-I)/FT(I3-I)/FT(I-J2)/FT(I)
      Q=X*       SQRT(Q1*FT(J1)*FT(K2)*FT(K3)*FT(I3)/FT(I1)/FT(I2)) 
      IF(MOD(IS,2).EQ.1)Q=-Q
      RETURN
    8 Q=0  
      PRINT1010,L1,L2,L3,M1,M2,M3 
 1010 FORMAT(10H ERREUR 3J,2(3X,3I3)) 
      RETURN
   11 IF(MOD(L,2).EQ.1)RETURN 
      I1=L-I1 
      I2=L-L2 
      I3=L-L3 
      Q=SQRT(FT(I1)*FT(I2)*FT(I3)/FT(L+1)*Q1) 
      I1=I1/2 
      I2=I2/2 
      I3=I3/2 
      L =L/2
      Q=Q*FT(L )/FT(I1)/FT(I2)/FT(I3) 
      IF(MOD(L +IS,2).EQ.1)Q=-Q 
      RETURN
      ENTRY TROISJI (L1,L2,L3,M1,M2,M3,Q) 
      ENTRY TROISJ (L1,L2,L3,M1,M2,M3,Q)  
      K3=M3 
      Q1=1 
      IS=0  
      GOTO12
      END 
*SIXJ 
      SUBROUTINE SIXJ(J1,J2,J3,L1,L2,L3,Q)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION M(7),M1(4),M2(4),M3(4),FT(0:1000)
      LOGICAL ICE 
      COMMON/TEUFJ/ICE
      SAVE  
      DATA LMEM,(FT(I),I=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/  
      ENTRYSIXJI(J1,J2,J3,L1,L2,L3,Q) 
      I1=J1 
      I2=J2 
      I3=J3 
      K1=L1 
      K2=L2 
      K3=L3 
      IS=0  
   24 Q=0  
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      DO 17 I=1,4 
      IF(MOD(M(I),2).EQ.1) GO TO 8
   17 CONTINUE
      L=MAX(I1+I2+K1+K2,I1+I3+K1+K3,I2+I3+K2+K3) 
      L=L/2 
      IF(L.LE.LMEM) GO TO 6 
      DO 10 I=LMEM,L  
   10 FT(I+1)=FT(I)*(I+1) 
      LMEM=L
    6 IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3) RETURN 
      IF(I1.LT.ABS(K2-K3).OR.I1.GT.K2+K3) RETURN 
      IF(K1.LT.ABS(I2-K3).OR.K1.GT.I2+K3) RETURN 
      IF(K1.LT.ABS(K2-I3).OR.K1.GT.K2+I3) RETURN 
      IF(I1) 8,2,1
    2 IF(I2.LT.0) GO TO 8 
    9 IF(I3.LT.0) GO TO 8 
   14 IF(K1.LT.0) GO TO 8 
   19 IF(K2.LT.0) GO TO 8 
   23 IF(K3.LT.0) GO TO 8 
   27 Q=SQRT(1.0d0/(I2+1)/(K2+1))  
      IS=(I2+K2+K1)/2+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    1 IF(I1.GT.1) GO TO 3 
      IF(I2.LT.0) RETURN
   12 IF(I3.LT.0) RETURN
   16 IF(K1.LT.0) RETURN
   21 IF(K2.LT.0) RETURN
   25 IF(K3.LT.0) RETURN
   28 IF(I2.LT.I3) GO TO 4
      IC=I2 
      I2=I3 
      I3=IC 
      IC=K2 
      K2=K3 
      K3=IC 
    4 IF(K2.GT.K3) GO TO 5
      I11=I1+K1+I2-K2 
      I11=I11/2 
      I12=I11-I2+K2 
      Q=SQRT(I11*I12*1./I3/(I3+1)/K3/(K3+1))
      IS =I11+K2+IS 
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    5 I11=K3-K1+I2
      I11=I11/2+1 
      I12=I11+K1+1
      Q=SQRT(I11*I12*1.0d0/I3/(I3+1)/K2/(K2+1))
      IS =I12-1+IS
      IF(MOD(IS ,2).EQ.1) Q=-Q  
      RETURN
    3 IF(I2.GE.I1) GO TO 7
      IF(I2.LT.0) GO TO 8 
      IC=I2 
      I2=I1 
      I1=IC 
      IC=K1 
      K1=K2 
      K2=IC 
      IF(I1.EQ.0) GO TO 9 
      IF(I1.EQ.1) GO TO 12
    7 IF(I3.GE.I1) GO TO 13 
      IF(I3.LT.0) GO TO 8 
      IC=I3 
      I3=I1 
      I1=IC 
      IC=K3 
      K3=K1 
      K1=IC 
      IF(I1.EQ.0) GO TO 14
      IF(I1.EQ.1) GO TO 16
   13 IF(K1.GE.I1) GO TO 18 
      IF(K1.LT.0) GO TO 8 
      IC=K1 
      K1=I1 
      I1=IC 
      IC=K2 
      K2=I2 
      I2=IC 
      IF(I1.EQ.0) GO TO 19
      IF(I1.EQ.1) GO TO 21
   18 IF(K2.GE.I1) GO TO 22 
      IF(K2.LT.0) GO TO 8 
      IC=K2 
      K2=I1 
      I1=IC 
      IC=K1 
      K1=I2 
      I2=IC 
      IF (I1.EQ.0) GO TO 23 
      IF(I1.EQ.1) GO TO 25
   22 IF(K3.GE.I1) GO TO 26 
      IF(K3.LT.0) GO TO 8 
      IC=K3 
      K3=I1 
      I1=IC 
      IC=K1 
      K1=I3 
      I3=IC 
      IF(I1.EQ.0) GO TO 27
      IF(I1.EQ.1) GO TO 28
   26 M1(4)=I3
      M1(1)=I3
      M1(3)=K3
      M1(2)=K3
      M2(2)=I1
      M2(1)=I1
      M2(4)=K1
      M2(3)=K1
      M3(3)=I2
      M3(1)=I2
      M3(4)=K2
      M3(2)=K2
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      Q1=1 
      DO 11 I=1,4 
      M(I)=M(I)/2 
   11 Q1=FT(M(I)-M1(I))*FT(M(I)-M2(I))*FT(M(I)-M3(I))*Q1/FT(M(I)+1) 
      Q1=SQRT(Q1) 
      M1(1)=I1+K1 
      M1(2)=I2+K2 
      M1(3)=I3+K3 
      IC=M1(1)+M1(2)  
      M(5)=IC/2 
      IC=M1(2)+M1(3)  
      M(6)=IC/2 
      IC=M1(1)+M1(3)  
      M(7)=IC/2 
      MAXZ=MIN(M(5),M(6),M(7)) 
      MINZ=MAX(M(1),M(2),M(3),M(4))
      X=0  
      DO 15 I=MINZ,MAXZ 
      Q2=1 
      DO 20 J=1,7 
      IJ=I-M(J) 
      IF(J.GT.4) IJ=-IJ 
   20 Q2=Q2*FT(IJ)
      Q2=FT(I+1)/Q2 
   15 X=-X+Q2 
      Q=X*Q1
      IS=MAXZ+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    8 PRINT 1010,J1,J2,J3,L1,L2,L3
      ICE=.FALSE. 
 1010 FORMAT(10H ERREUR 6J,2(3X,3I3)) 
      RETURN
      ENTRY RACAH(J1,J2,J3,L1,L2,L3,Q)
      IS=(J1+J2+J3+L1)/2
      I1=J1 
      I2=J2 
      I3=L2 
      K1=L1 
      K2=J3 
      K3=L3 
      GO TO 24
      END 

