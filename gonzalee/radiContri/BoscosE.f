*module cstes
      MODULE cstes
c Physical constants: hc=hbar*c; e2=e^2/4 pi epsilon_0;hm=hbar^2/2m_N; ci=i
      real*8 :: pi,hc,e2,hm
      complex*16, parameter :: ci=(0,1)
      END MODULE
*module nucleus
      MODULE nucleus
c Data about the nucleus: spin of the fragment, mass and charge of the
c core(c) and fragment(f), reduced c-f mass,...
      integer :: iS,iSc,kZf,kZc,kZT
      real*8 :: Af,Ac,Apr,deuxm
      END MODULE
*module potential
      MODULE potential
c Parameters of the two-body potential Vcf
      integer  :: ntypo,Njpo
      integer, allocatable, dimension (:) :: Jpo,lpo
      real*8, allocatable, dimension (:) :: Vp,rp,ap,Vls,rls,als,rC
      END MODULE

      program Boscos
*     This Fortran code computes the eigenstates of a single-particle
*     Hamiltonian for two particles interacting through a local potential.
*     The programmed potentials are of Woods-Saxon or Gaussian forms.
*     In the former case a spin-orbit coupling term is also possible.
*     The Coulomb term has a point-sphere expression.
*     Both bound and scattering states can be computed.
*     In the latter case, phaseshifts can be provided.
*     The output options are controlled by the 'lec' variable
*     at the end of the input file.

*     new version of 'elamain4.f', in which the second Gaussian term is it's second derivative (r^2 * Gaussian)
*     new version of 'elamain5.f', in which forbidden states can be removed by susy
      USE cstes
      USE nucleus
      USE potential

      implicit none
      integer :: nAc,nAf,lzc,lzf
      integer :: Nel,Neli,NE
      integer, allocatable, dimension (:) :: NRel,Jel,lel,iEm
      integer, allocatable, dimension (:) :: NReli,Jeli,leli
      real*8, allocatable, dimension (:) :: E0,E0i,V,U,S,psm,dps,dpsm,SB
      real*8, allocatable, dimension (:,:) :: ps
      character(len=30) :: fich
      integer :: lec,NNR,iJ,l,ipo,iel,ieli,Nru,jr,iE,ins
      real*8 :: E,Ek,eta,rmu,hu,r
      real*8 :: Emin,hE,delta,rfin
      real*8 :: S2,al,am,z,whit,whit1,eps,vJ,r2,proba,anc
      real*8 :: f(0:20),fp(0:20),g(0:20),gp(0:20)
      real*8, parameter :: EPSW=1d-6

c     Input reading
c Math constants
      pi=acos(-1d0)
c Physical constants: hc=hbar*c; e2=e^2/4 pi epsilon_0;hm=hbar^2/2m_N
      read(*,*)hc,e2,hm
c Mass, charge, and spin of core(c) and the fragment(f)
      read(*,*)Ac,Af,kZc,kZf,iSc,iS
      Apr=Af+Ac
      deuxm=(Af*Ac)/Apr/hm      !=2\mu/\hbar^2
*      deuxm=Af/hm      !=2mN/\hbar^2
c Spin of f and c (up to now Sc is set equal to 0)
      nAc=Ac+0.5d0
      nAf=Af+0.5d0
      lzc=2-mod(nAc,2)
      if(iS==0)then
c     fragment spin is neglected
         lzf=2
      else
         lzf=2-mod(nAf,2)
      endif

      iSc=iSc*lzc
      iS=iS*lzf

c c-f potential
      read(*,*)ntypo,Njpo
c ntypo<10: Woods Saxon
c ntypo<20: Gaussian
c     X1: different potentials for Njpo-1 partial waves + 1 for the others
c     X2: different potentials for l-even and l-odd partial waves

      if(ntypo/=1.and.ntypo/=2.and.ntypo/=11.and.ntypo/=12)then
         write(*,*)'Unprogrammed potential'
         STOP
      endif
      if(mod(ntypo,10)==2) Njpo=2

      allocate(Jpo(Njpo),lpo(Njpo))
      allocate(Vp(Njpo),rp(Njpo),ap(Njpo),rC(Njpo))
      allocate(Vls(Njpo),rls(Njpo),als(Njpo))

      If(mod(ntypo,10)==1)then
         do ipo=1,Njpo-1
            read(*,*)Jpo(ipo),lpo(ipo),Vp(ipo),rp(ipo),ap(ipo),
     a           Vls(ipo),rls(ipo),als(ipo),rC(ipo)
         enddo
         Jpo=Jpo*lzf
         read(*,*)Vp(Njpo),rp(Njpo),ap(Njpo),
     a        Vls(Njpo),rls(Njpo),als(Njpo),rC(Njpo)
      ElseIf(mod(ntypo,10)==2)then
            read(*,*)Vp(1),rp(1),ap(1),Vls(1),rls(1),als(1),rC(1)!1=odd l
            read(*,*)Vp(2),rp(2),ap(2),Vls(2),rls(2),als(2),rC(2)!2=even l
      EndIf

c Bound states to be searched for (etats lie: el)
      read(*,*)Nel
      allocate(NRel(Nel),Jel(Nel),lel(Nel),E0(Nel))
      read(*,*)(NRel(iel),Jel(iel),lel(iel),iel=1,Nel)
      Jel=Jel*lzf

c Forbidden states to be removed by SuSy
*      read(*,*)Neli
      Neli=0
*      allocate(NReli(Neli),Jeli(Neli),leli(Neli),E0i(Neli))
*      read(*,*)(NReli(ieli),Jeli(ieli),leli(ieli),ieli=1,Neli)
*      Jeli=Jeli*lzf

c Numerical parameters: Nru: Nb of point in radial uniform mesh
c rmu: maximum radius in the mesh (in fm), eps: accuracy requested
      read(*,*) Nru,rmu,eps
      allocate(U(Nru),S(Nru),V(Nru))
      hu=rmu/Nru

      write(*,7000)Ac,Af
 7000 format(' Ac=',f8.3,' Af=',f7.3)
      write(*,7001)kZc,kZf
 7001 format(' Zc=',i3,' Zf=',i3)
      write(*,7002)deuxm,1/deuxm
 7002 format(' deuxm=',1pd16.8,' MeV^-1fm^-2 or 1/deuxm=',
     a     2pd16.8,' MeV fm^2')
      
c Energies: NE=Nb of energies (uniform mesh), Emin=first energy, hE=step
c (in MeV)
      read(*,*) NE,Emin,hE

c Printing options
c lec=
c     0: computes all bound states and gives in output their energy and ANC
c     10:prints wf of first state in file fich.dfo
c     100: phaseshifts of the given partial waves
c          NE=1: print wf of first state and first energy in fich.dep
c     110: build a bin wave function on the energy range given
      read(*,*)lec
      read(*,*)fich
      if(lec==10.or.lec==100.and.NE==1.or.lec==110)then
         fich=trim(fich)//'.dfo'
      elseif(lec==100)then
         fich=trim(fich)//'.dep'
      endif

      if(lec==100)then
         allocate(ps(Nel,NE),psm(Nel),dps(Nel),dpsm(Nel),iEm(Nel))
      endif

      DO iel=1,Nel
         NNR=NRel(iel)
         iJ=Jel(iel)
         l=lel(iel)
c Computing potential
         Do jr=1,Nru
            r=jr*hu
            V(jr)=vJ(r,l,iJ)*deuxm
         End Do

c Use a buffer U because some subroutines modify the potential
         U=V

**************************************
c Not used in this version of the code
         Do ieli=1,Neli
            if(l==leli(ieli).AND.iJ==Jeli(ieli))then
               ins=l+1
               call super(Nru,hu,ins,E0i(ieli),NReli(ieli),U,S,eps)
            endif
         EndDo !ieli
**************************************

         IF(lec<100)then !bound states computation
            write(*,7003)l,iJ,NNR
 7003       format(' l=',i3,' J=',i3,' NR=',i3)
c Call 'Num1l'
            E=0
            S2=0
            call num1l(Nru,hu,E,S2,U,S,NNR,eps)
            If(NNR<0)then
               write(*,*)' No solution found for bound state ',iel
            Else
               E0(iel)=E/deuxm
               write(*,7004)E0(iel),E
 7004          format(' E=',f13.8,' MeV; K^2=',f10.5,' fm^-2')
*               write(*,7005)sqrt(r2(Nru,hu,S))
* 7005          format(' <r^2>^{1/2}=',f10.5)
*               write(*,7006)proba(Nru,hu,U,E,S)
* 7006          format(' Proba r>turning point=',f10.5)
*
               Ek=sqrt(-E)
               eta=kZf*kZc*e2*deuxm/Ek/2
               al=-eta
               am=l+0.5d0
               z=2*Ek*(rmu-100*hu)
               call witt(al,am,z,whit,whit1,EPSW)
               anc=S(Nru-100)/whit
               write(*,7007)anc
c calculation with expansion of Whitacker function
*     a           S(Nru)/(exp(-Ek*r)/(2*Ek*r)**eta
*     b           *(1+(l*(l+1)-eta*(eta+1))/(2*Ek*r)
*     c           +(l*(l+1)-eta*(eta+1))*(l*(l+1)-eta*(eta+3)-2)
*     d           /2/(2*Ek*r)**2))
 7007          format(' ANC=',f13.8,' fm^-1/2')

               if(iel==1.and.lec>=10)then
c writes wavefunction of first state, its asymptotic behaviour
c and the potential (in MeV)
                  open(unit=1,status='unknown',file=fich)
                  Do jr=1,Nru
                     r=jr*hu
                     z=2*Ek*r
                     call witt(al,am,z,whit,whit1,EPSW)
                     write(1,'(1f8.3,3d17.8)')r,S(jr),anc*whit,
     a                    U(jr)/deuxm
                  End Do
                  close(unit=1)
               endif
            EndIf
         ELSEIF(lec<110)then       !computation of continuum state(s)

            If(NE==1)then
c give wave function of the continuum state
c at that sole energy in the first partial wave listed,
c its asymptotic behaviour, and the corresponding plane wave. 
               E=Emin
               Ek=sqrt(deuxm*E)
               eta=kZc*kZf*e2*deuxm/Ek/2
               call dephase0(Nru,hu,U,S,eps,delta,l,eta,Ek,rfin)
               open(unit=1,file=fich)
               do jr=1,Nru,10
                  z=Ek*jr*hu
                  call coufra(z,eta,l,l,f,fp,g,gp)
                  write(1,'(1f8.3,3d17.8)')jr*hu,S(jr),
     a                 cos(delta)*f(l)+sin(delta)*g(l),f(l)
               enddo
               close(1)
               write(*,8001)l,iJ,Emin,delta,delta*180/pi
 8001          FORMAT (' l= ',i3,' iJ= ',i3,' E= ',1f8.4,
     a              ' MeV; phaseshift= ',1d18.8,' rad; ',1d18.8,' deg')
            Else
c computation of phase shift as a function of energy, its energy derivative
               Do iE=1,NE
                  E=Emin+(iE-1)*hE
                  if(E==0)then
                     ps(iel,iE)=0
                  else
                     Ek=sqrt(deuxm*E)
                     eta=kZc*kZf*e2*deuxm/Ek/2
                     call dephase0(Nru,hu,V,S,eps,
     a                    ps(iel,iE),l,eta,Ek,rfin)
                     If(iE>1)then
                        do
                           if(abs(ps(iel,iE)-ps(iel,iE-1))<2.5d0) exit
             if(ps(iel,iE)-ps(iel,iE-1)>1d0) ps(iel,iE)=ps(iel,iE)-pi
             if(ps(iel,iE)-ps(iel,iE-1)<-1d0) ps(iel,iE)=ps(iel,iE)+pi
                        enddo
                     EndIf
                  endif
c if one wants to extract scattering lengths:
*               write(*,8002)E,psm1,al(Ek,psm1,eta,l)
* 8002          FORMAT (1f8.4,3d18.8)
               EndDo
            EndIf
         ELSE
c building of the bin wave function (with unit weight) for lec==110
            allocate(SB(Nru))
            SB=0
            Do iE=1,NE
               E=Emin+(iE-1)*hE
               Ek=sqrt(deuxm*E)
               eta=kZc*kZf*e2*deuxm/Ek/2
               call dephase0(Nru,hu,V,S,eps,delta,l,eta,Ek,rfin)
               SB=SB+S
            EndDo
            SB=SB/sqrt(NE*hE)
            open(unit=1,file=fich)
            do jr=1,Nru,10              
               write(1,'(1f8.3,1d17.8)')jr*hu,SB(jr)
            enddo
            close(1)
         ENDIF!lec
      ENDDO !iel

      If(lec==100.and.NE>1)then
c     writes the phaseshift, derivative and contribution to elastic-
c     scattering cross section for each of the requested partial waves
         open(unit=1,status='unknown',file=fich)
         Do iE=1,NE
            E=Emin+(iE-1)*hE
            Ek=sqrt(deuxm*E)
            do iel=1,Nel
               if(iE==1)then
                  dps(iel)=-3*ps(iel,1)+4*ps(iel,2)-ps(iel,3)
                  dpsm(iel)=0
                  psm(iel)=0
                  iEm(iel)=0
               elseif(iE==NE)then
                  dps(iel)=3*ps(iel,NE)-4*ps(iel,NE-1)+ps(iel,NE-2)
               else
                  dps(iel)=ps(iel,iE+1)-ps(iel,iE-1)
                  if(dpsm(iel)<dps(iel))then
                     dpsm(iel)=dps(iel)
                     psm(iel)=ps(iel,iE)
                     iEm(iel)=iE
                  endif
               endif
            enddo
            dps=dps/2/hE
            write(1,'(1f8.3,12(1pd15.6))')E,(ps(iel,iE),dps(iel),
     a           -tan(ps(iel,iE))/Ek**(2*l+1),iel=1,Nel)
         EndDo
         dpsm=dpsm/2/hE
         close(1)

c Print location of eventual resonance
c From the derivative of the phaseshift with energy, provides the location
c and width of a possible resonance.
         do iel=1,Nel
            if(iEm(iel)<NE.and.iEm(iel)>2)then
               write(*,8000)l,iJ
 8000          format(' l=',i3,' J=',i3)
               write(*,8003)Emin+(iEm(iel)-1)*hE,2/dpsm(iel)
 8003          format(' Eres=',f10.6,' Gamma=',f10.6)
            endif
         enddo
         deallocate(ps,psm,dps,dpsm,iEm)
      EndIf

      deallocate(Jpo,lpo,Vp,rp,ap,rC,Vls,rls,als)
      deallocate(NRel,Jel,lel,E0,U,S,V)
 9999 continue

      END

*vJ
      function vJ(r,l,iJ)
c Effective potential between the core and the fragment in the l J channel
c It contains nuclear (Woods-Saxon), Coulomb (point-sphere) and centrifugal
      USE cstes
      USE nucleus
      USE potential
      implicit real*8 (a-h,o-z)

      jls=iJ*(iJ+2)-4*l*(l+1)-iS*(iS+2)

      IF(mod(ntypo,10)==1)then      !Different potentials for various partial waves
         ipo=1
         do
            if(l==lpo(ipo).and.iJ==Jpo(ipo).or.ipo==Njpo)then
c     special channels or all others
               exit
            else
               ipo=ipo+1
            endif
         enddo
      ELSEIF(mod(ntypo,10)==2)then      !Different potentials for odd and even waves
         ipo=2-mod(l,2)
      ELSE
         write(*,*)'Uprogrammed potential'
         STOP
      ENDIF

      IF(ntypo<10)then
         eve=exp(-(r-rp(ipo))/ap(ipo)) !Wodds-Saxon form factor
         vJ=Vp(ipo)*eve/(1+eve)+l*(l+1)/r**2/deuxm
         if(Vls(ipo)/=0)then
            eve=exp(-(r-rls(ipo))/als(ipo))
            vJ=vJ-jls*Vls(ipo)*eve/(1+eve)**2/als(ipo)/r/8
         endif
      ELSEIF(ntypo<20)then
         eve=(r-rp(ipo))/ap(ipo)
         eve=exp(-eve**2/2) !Gaussian form factor
         vJ=Vp(ipo)*eve+l*(l+1)/r**2/deuxm
         if(Vls(ipo)/=0)then !the ls term is used for the derivative
            eve=(r-rls(ipo))/als(ipo)
            eve=(r-rls(ipo))**2*exp(-eve**2/2) !(r-r0)^2*Gaussian form factor
            vJ=vJ+Vls(ipo)*eve
         endif
      ELSE
         write(*,*)'Uprogrammed potential'
         STOP
      ENDIF

      If (kZf/=0) then
c Point-sphere Coulomb potential
         if(r>=rC(ipo))then
            c1=kZf*kZc*e2/r
            vJ=vJ+c1
         else
            c1=kZf*kZc*e2/2/rC(ipo)
            c2=3-(r/rC(ipo))**2
            vJ=vJ+c1*c2
         endif
      Endif
      return
      end
*r2
      function r2(Nr,hr,S)
      implicit real*8(a-h,o-z)
      real*8 ::S(Nr)
      r2=0
      Do ir=1,Nr
         r2=r2+(S(ir)*ir*hr)**2
      End Do
      r2=r2*hr
      return
      end
*proba
      function proba(Nr,hr,U,E,S)
      implicit real*8(a-h,o-z)
      real*8 ::U(Nr),S(Nr)
      jr=1
      Do
         If(U(jr)<E.and.U(jr+1)>E)then
            EXIT
         Else
            jr=jr+1
         Endif
      EndDo
      proba=S(jr)**2/2
      Do ir=jr+1,Nr
         proba=proba+S(ir)**2
      End Do
      proba=proba*hr
      return
      end

*super
      subroutine super (n,h,ins,e0,nrad,u,s,eps)
*     n: number of discretization points
*     h: step
*     ins: singularity modification at the origin
*     e0: approximation of the bound-state energy
*     nrad: number of nodes
*     u: potential multiplied by 2m/hbar^2
*     s: auxiliary vector (wave function, potential modification)
*     eps: accuracy
*     
      implicit real*8 (a-h,o-z)
      dimension u(n),v(n),s(n)

      e=e0
      v=u                       !utilisation d'un vecteur tampon
      s2=0 
*
*     calculation of bound-state wave function
*     s: wave function
*
      call num1l(n,h,e,s2,v,s,nrad,eps)
      write(*,*)'Energie de l''etat lie dans super:',e
*
*     integration of psi^2 by Adams formula
*     (Abramowicz and Stegun: 25.5.5)
*
      s1=s(3) 
      s2=s(2) 
      s3=s(1) 
      k=2*ins+1
      s(1)=log(h*s(1)**2/k)
      ss=2*h*s(2)**2/k
      s(2)=log(ss) 
      ss=ss+h/24*(9*s(3)**2+19*s2**2-5*s3**2) 
      s(3)=log(ss) 
      do 3 j=4,n
         ss=ss+h/24*(9*s(j)**2+19*s1**2-5*s2**2+s3**2) 
         s3=s2 
         s2=s1 
         s1=s(j) 
         s(j)=log(ss) 
 3    continue
*
*     5-point second derivative of the logarithm of the integral
*     (Abramowicz and Stegun: 25.3.24)
*
      s2=s(1)
      s1=s(2)
      s(1)=-k/h**2
      s(2)=s(1)/4
      do 4 j=5,n
         s3=s2
         s2=s1
         s1=s(j-2)
         s(j-2)=(-(s(j)+s3)+16*(s(j-1)+s2)-30*s1)/(12*h**2)
 4    continue
      s(n-1)=0
      s(n)=0
      do 5 j=1,n
         u(j)=u(j)-s(j)*2
 5    continue
      ins=ins+2
      return
      end
      function fcth(eta)
      implicit none
      integer, parameter :: nmax=100
      real*8, parameter :: euler=0.5772156649d0
      integer n
      real*8 fcth,eta,eta2,fn
      eta2=eta**2
      fn=1+eta2
      fcth=1/fn
      do n=2,nmax
         fn=n**3+eta2*n
         fcth=fcth+1/fn
      enddo
      fcth=eta2*fcth
      fcth=fcth-euler-log(eta)
      return
      END

      subroutine fctC(l,eta,C0,Cl)
      USE cstes
      implicit none
      real*8 eta,C0,Cl,dpe,eta2
      integer l,i,l2
      dpe=2*pi*eta
      C0=exp(dpe)-1
      C0=dpe/C0
      if(l==0)then
         Cl=C0
      else
         eta2=eta**2
         Cl=C0*(1+eta2)
         l2=1
         do i=2,l
            l2=l2+2*i-1
            Cl=Cl*(1+eta2/l2)
         enddo
      endif
      END

      function al(Ek,delta,eta,l)
c Computes the scattering "length" of partial wave l from
c phaseshift delta at energy Ek
      implicit none
      integer l
      real*8 al,delta,eta,fcth,fh,C0,Cl,Ek
      If(eta/=0)then
         fh=fcth(eta)
         call fctC(l,eta,C0,Cl)
         al=(1d0/tan(delta)+2*eta*fh/C0)*Cl*Ek**(2*l+1)
         al=-1/al
      Else
         al=-tan(delta)/Ek**(2*l+1)
      EndIf
      return
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
 2000 FORMAT(/27H PAS DE CONVERGENCE POUR L=,I3,5H   K=,F10.6/) 
 2001 FORMAT(/46H LE POTENTIEL N'EST JAMAIS PUREMENT COULOMBIEN)
      END

*NUM1L
      SUBROUTINE NUM1L(N,H,E,S2,U,S,NO,EPS) 
CC  VERSION CORRIGEE LE 21 NOV 72 
C*****INTEGRATION DE L"EQUATION DE SCHROEDINGER PAR LA METHODE DE NUMERO
C*****POUR E NEGATIF  
C*****RECHERCHE DE L"ENERGIE PROPRE PAR LA METHODE DE RAPHSON-NEWTON
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION U(1),S(1) 
*      DIMENSION U(100000),S(100000) 
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

*WITT
      SUBROUTINE WITT(AL,AM,X,W,WP,EPS)
      IMPLICIT REAL*8 (A-H,O-Z) 
      CALL WIT1(AL-1,AM,X,W1,EPS)
      CALL WIT1(AL-2,AM,X,W2,EPS)
      W=(X-2*(AL-1))*W1+(AM+1.5d0-AL)*(AM-1.5d0+AL)*W2
      WP=(AL/X-0.5d0)*W-(AM-AL+0.5d0)*(AM+AL-0.5d0)*W1/X
      RETURN
      END

*WIT1
      SUBROUTINE WIT1(AL,AM,X,W,EPS)
      IMPLICIT REAL*8 (A-H,O-Z) 
      AP=0.5d0+AM-AL
      BP=1+AM+AM
      W1=0
      DO1 L=1,7
      CALL WIT(L,AP,BP,X,W2)
      IF(ABS((W1-W2)/W2).LT.EPS)GOTO2
      W1=W2
    1 CONTINUE
      PRINT1000,AL,AM,X,EPS
    2 W=W2*X**AL*EXP(-X/2)/GAMA(AP)
      RETURN
 1000 FORMAT('PAS DE CONVERGENCE COWITT,AL,AM,X,EPS=',1P,4d13.5)
      END

*WIT
      SUBROUTINE WIT(L,AP,BP,X,W)
      IMPLICIT REAL*8 (A-H,O-Z) 
      dimension KLA(8),AI(136),R(136)
      DATA KLA/ 0, 4,12,24,40,64,96,136/
      DATA (AI(I),I=  1,52 )/
     1  3.2254768961939D-01, 1.7457611011584D+00, 4.5366202969211D+00,
     1  9.3950709123010D+00,
     1  1.7027963230510D-01, 9.0370177679937D-01, 2.2510866298660D+00,
     1  4.2667001702877D+00, 7.0459054023936D+00, 1.0758516010181D+01,
     1  1.5740678641278D+01, 2.2863131736889D+01,
     1  1.1572211735804D-01, 6.1175748451513D-01, 1.5126102697764D+00,
     1  2.8337513377435D+00, 4.5992276394182D+00, 6.8445254531150D+00,
     1  9.6213168424567D+00, 1.3006054993306D+01, 1.7116855187462D+01,
     1  2.2151090379396D+01, 2.8487967250984D+01, 3.7099121044467D+01,
     1  8.7649410478933D-02, 4.6269632891507D-01, 1.1410577748312D+00,
     1  2.1292836450983D+00, 3.4370866338932D+00, 5.0780186145497D+00,
     1  7.0703385350481D+00, 9.4383143363917D+00, 1.2214223368866D+01,
     1  1.5441527368781D+01, 1.9180156856753D+01, 2.3515905693991D+01,
     1  2.8578729742881D+01, 3.4583398702284D+01, 4.1940452647686D+01,
     1  5.1701160339542D+01,
     1  5.9019852181480D-02, 3.1123914619844D-01, 7.6609690554590D-01,
     1  1.4255975908036D+00, 2.2925620586321D+00, 3.3707742642089D+00,
     1  4.6650837034671D+00, 6.1815351187366D+00, 7.9275392471720D+00,
     1  9.9120980150776D+00, 1.2146102711729D+01, 1.4642732289596D+01/
      DATA (AI(I),I= 53,108)/
     1  1.7417992646508D+01, 2.0491460082615D+01, 2.3887329848169D+01,
     1  2.7635937174332D+01, 3.1776041352373D+01, 3.6358405801649D+01,
     1  4.1451720484868D+01, 4.7153106445154D+01, 5.3608574544692D+01,
     1  6.1058531447214D+01, 6.9962240035101D+01, 8.1498279233948D+01,
     1  4.4489365833263D-02, 2.3452610951958D-01, 5.7688462930186D-01,
     1  1.0724487538177D+00, 1.7224087764446D+00, 2.5283367064257D+00,
     1  3.4922132730219D+00, 4.6164567697496D+00, 5.9039585041741D+00,
     1  7.3581267331861D+00, 8.9829409242124D+00, 1.0783018632539D+01,
     1  1.2763697986742D+01, 1.4931139755522D+01, 1.7292454336715D+01,
     1  1.9855860940335D+01, 2.2630889013196D+01, 2.5628636022458D+01,
     1  2.8862101816323D+01, 3.2346629153963D+01, 3.6100494805750D+01,
     1  4.0145719771537D+01, 4.4509207995752D+01, 4.9224394987306D+01,
     1  5.4333721333392D+01, 5.9892509162130D+01, 6.5975377287930D+01,
     1  7.2687628090655D+01, 8.0187446977905D+01, 8.8735340417883D+01,
     1  9.8829542868275D+01, 1.1175139809793D+02,
     1  3.5700394308758D-02, 1.8816228315858D-01, 4.6269428131451D-01,
     1  8.5977296397281D-01, 1.3800108205271D+00, 2.0242091359225D+00,
     1  2.7933693535066D+00, 3.6887026779080D+00, 4.7116411465546D+00,
     1  5.8638508783433D+00, 7.1472479081018D+00, 8.5640170175858D+00/
      DATA (AI(I),I=109,136)/
     1  1.0116634048452D+01, 1.1807892294004D+01, 1.3640933712536D+01,
     1  1.5619285893338D+01, 1.7746905950095D+01, 2.0028232834574D+01,
     1  2.2468249983498D+01, 2.5072560772425D+01, 2.7847480009168D+01,
     1  3.0800145739444D+01, 3.3938657084912D+01, 3.7272245880475D+01,
     1  4.0811492823884D+01, 4.4568603175330D+01, 4.8557763533054D+01,
     1  5.2795611187213D+01, 5.7301863323389D+01, 6.2100179072770D+01,
     1  6.7219370927125D+01, 7.2695158847608D+01, 7.8572802911565D+01,
     1  8.4911231135695D+01, 9.1789874671229D+01, 9.9320808717441D+01,
     1  1.0767244063938D+02, 1.1712230951268D+02, 1.2820184198825D+02,
     1  1.4228004446916D+02/
      DATA (R(I),I=  1,52 )/
     1  8.3273912383786D-01, 2.0481024384543D+00, 3.6311463058215D+00,
     1  6.4871450844076D+00,
     1  4.3772341049290D-01, 1.0338693476656D+00, 1.6697097656587D+00,
     1  2.3769247017586D+00, 3.2085409133480D+00, 4.2685755108251D+00,
     1  5.8180833686718D+00, 8.9062262152922D+00,
     1  2.9720963604447D-01, 6.9646298043062D-01, 1.1077813946158D+00,
     1  1.5384642390428D+00, 1.9983276062742D+00, 2.5007457691008D+00,
     1  3.0653215182824D+00, 3.7232891107827D+00, 4.5298140299818D+00,
     1  5.5972584618353D+00, 7.2129954609262D+00, 1.0543837461910D+01,
     1  2.2503631486429D-01, 5.2583605276235D-01, 8.3196139168709D-01,
     1  1.1460992409637D+00, 1.4717513169668D+00, 1.8131346873813D+00,
     1  2.1755175196946D+00, 2.5657627501650D+00, 2.9932150863713D+00,
     1  3.4712344831021D+00, 4.0200440864446D+00, 4.6725166077329D+00,
     1  5.4874206579862D+00, 6.5853612332889D+00, 8.2763579843641D+00,
     1  1.1824277551659D+01,
     1  1.5149441285947D-01, 3.5325658252990D-01, 5.5678456328812D-01,
     1  7.6268531769732D-01, 9.7187263224651D-01, 1.1853578930378D+00,
     1  1.4042656272844D+00, 1.6298686157570D+00, 1.8636350553320D+00,
     1  2.1072911510814D+00, 2.3629058910419D+00, 2.6330087531638D+00/
      DATA (R(I),I= 53,108)/
     1  2.9207575797277D+00, 3.2301851334923D+00, 3.5665733773686D+00,
     1  3.9370437554551D+00, 4.3515311888634D+00, 4.8244818548978D+00,
     1  5.3780220797891D+00, 6.0484178126197D+00, 6.9008983521807D+00,
     1  8.0699651561471D+00, 9.9027933194843D+00, 1.3820532094793D+01,
     1  1.1418710576812D-01, 2.6606521689761D-01, 4.1879313732490D-01,
     1  5.7253284649977D-01, 7.2764878838098D-01, 8.8453671934022D-01,
     1  1.0436188758921D+00, 1.2053492741523D+00, 1.3702213385218D+00,
     1  1.5387772564687D+00, 1.7116193526864D+00, 1.8894240634494D+00,
     1  2.0729593402465D+00, 2.2631066339970D+00, 2.4608890724883D+00,
     1  2.6675081263971D+00, 2.8843920929220D+00, 3.1132613270395D+00,
     1  3.3562176925959D+00, 3.6158698564843D+00, 3.8955130449484D+00,
     1  4.1993941047115D+00, 4.5331149785341D+00, 4.9042702876116D+00,
     1  5.3235009720236D+00, 5.8063332142336D+00, 6.3766146741598D+00,
     1  7.0735265807070D+00, 7.9676935092955D+00, 9.2050403312780D+00,
     1  1.1163013090768D+01, 1.5390180415262D+01,
     1  9.1625471157176D-02, 2.1342058490497D-01, 3.3571811668027D-01,
     1  4.5854093503343D-01, 5.8206816577905D-01, 7.0649521636717D-01,
     1  8.3202690300345D-01, 9.5887819879439D-01, 1.0872761620305D+00,
     1  1.2174623279777D+00, 1.3496954913567D+00, 1.4842549297768D+00/
      DATA (R(I),I=109,136)/
     1  1.6214441628118D+00, 1.7615953746767D+00, 1.9050746658947D+00,
     1  2.0522883472617D+00, 2.2036905532450D+00, 2.3597925385232D+00,
     1  2.5211741403764D+00, 2.6884980554088D+00, 2.8625278132105D+00,
     1  3.0441506653115D+00, 3.2344070972636D+00, 3.4345293984277D+00,
     1  3.6459928249939D+00, 3.8705845972164D+00, 4.1104986804328D+00,
     1  4.3684687232542D+00, 4.6479589840745D+00, 4.9534461124097D+00,
     1  5.2908484059009D+00, 5.6682046090332D+00, 6.0967964147430D+00,
     1  6.5931088610398D+00, 7.1824959955372D+00, 7.9066663113853D+00,
     1  8.8408924928109D+00, 1.0140899265621D+01, 1.2210021299206D+01,
     1  1.6705520642027D+01/
      save
      NW1=KLA(L)+1
      NW2=KLA(L+1)
      W=0
      DO201 I=NW1,NW2
      XX=AI(I)
  201 W=W+R(I)*XX**(AP-1)*(1+XX/X)**(BP-AP-1)*EXP(-XX)
      RETURN
      END

*GAMMA
      FUNCTION GAMA (XX)
C 
C     GAMMA FUNCTION  
C---                              BY CHEBYSHEV POLYNOMIALS
C     LOGARITHM OF GAMMA FUNCTION 
C---                                             - CDC 6400 - 
C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION CT(20),B(5)
      SAVE  
      DATA B / 12,-360,1260,-1680,1188/
      DATA ALPID/ 0.91893853320467d0/ 
      DATA CK / 0.57721566490153d0/ 
      DATA PI /3.1415926535898d0/ 
C 
      DATA CT/1.8835711955910D0, 0.44153813248410D-2, 
     1        0.56850436815994  D-1,-0.42198353964186  D-2, 
     2        0.13268081812125  D-2,-0.18930245297989  D-3, 
     3        0.36069253274412  D-4,-0.60567619044609  D-5, 
     4        0.10558295463023  D-5,-0.18119673655426  D-6, 
     5        0.3117724964715   D-7,-0.535421963902    D-8, 
     6 0.91932755199D- 9,-0.15779412803D- 9, 0.2707980623 D-10, 
     7-0.464681865  D-11, 0.79733502   D-12,-0.13680782   D-12, 
     8 0.2347319    D-13,-0.402743     D-14/
C 
C *************************** GAMMA (X) *************************** 
C 
      ENTRY GAMMA(XX) 
      X=XX  
      LY   = 1
      IF (X.LE.0) GO TO 99 
  100 IF (X.LE.1.0d-6) GO TO 300
      IF (X.EQ.2) GO TO 10
      Y    = X - 1  
      LG    = 1 
      T     = Y 
      IF (Y) 1,10,2 
    1 IY    = ABS(Y)+1
      T     = Y + IY 
      LG    = 2 
      GO TO 3 
    2 IF (Y.LE.1) GO TO 3 
      IY    = Y 
      T     = Y - IY 
      LG    = 3 
    3 F     = 4*T - 2
      A2    = CT(20)  
      A1    = F * A2      + CT(19)
      DO 4 I=3,18,2 
      J     = 20 - I  
      A2    = F * A1 - A2 + CT(J+1) 
      A1    = F * A2 - A1 + CT(J) 
    4 CONTINUE
      A2    = F * A1 - A2 + CT(2) 
      A2    = F * A2 - A1 + CT(1) 
      A     = (A2-A1)/2     
      IF(LG-2)5,6,8 
    5 GAMA = A
      IF (LY.EQ.2) GO TO 200
      GOTO11
    6 D     = 1 
      DO 7 I=1,IY 
      D     = D * (T - (I-1))
    7 CONTINUE
      GAMA = A / D
      IF (LY.EQ.2) GO TO 200
      GOTO11
    8 D     = 1
      DO 9 I=1,IY 
      D     = D * (T + I)
    9 CONTINUE
      GAMA = A * D
      IF (LY.EQ.1) GOTO11 
  200 GAMA = LOG(GAMA)
      GOTO11
   10 GAMA = 1
      IF (LY.EQ.1) GOTO11 
      GAMA = 0
      GOTO11
C 
C ************************* LOG GAMMA (X) ************************* 
C 
      ENTRY ALGAMA(XX)
      X=XX  
      LY   = 2
      IF (X.LE.0) GO TO 99 
   20 IF (X.LT.15) GO TO 100 
      XL   = LOG(X) 
      GAMA = X*(XL-1)-XL/2+ALPID
      IF (X.GT.2**12) GOTO11
      XL   = X * X
      F    = X
      G    = 1/(X*B(1))
      DO 21 I=2,5 
      F    = XL * F 
      G    = G + 1/(F*B(I))
   21 CONTINUE
      GAMA = GAMA + G 
      GOTO11
  300 ZZ   = (CK * X + 1) * X 
      GAMA = 1/ZZ 
      IF(LY.EQ.2)GOTO200
      GOTO11
   99 X=-X+1 
      IF(X.EQ.INT(X))GOTO12
      IF(LY.EQ.2)GOTO20 
      GOTO100 
   11 IF(XX.GT.0)RETURN
      IF(LY.EQ.1)GAMA=PI/(GAMA*SIN(PI*X)) 
      IF(LY.EQ.2)GAMA=GAMA+LOG(PI/ABS(SIN(PI*X))) 
      RETURN
   12 PRINT1000,INT(XX) 
      GAMA=0 
      RETURN
 1000 FORMAT(8H GAMMA (,I4,15H ) N EXISTE PAS)
      END   
 
