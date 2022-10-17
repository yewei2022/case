      Parameter(nlvl=1,NC=1,NA=2,NR=2000,NF=10) !NLVL is level of height, NC represents aerosol bin, NA = 1 is ammonium sulfate, NA =2 is
!     ice nuclei, NR is hydrometeor's bin, NF is the number of Bessel expasion terms.  
      Common/dm/iunit
      integer iunit
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      Common/maths/pi,cc1,t0,xh,xk
      Common/sk2/dtd,dx,ax,sjrs
      Common/sk/dt,dz
      common/grid/en(nc),e(nr),rl(nr)
      common/distr/f(nr,nc,na,nlvl),fi(nr,nc,na,nlvl),f0(nr,nc,na,nlvl)
      common/switch/icoal,ice,irim,inuc,istop,jrain
      common/compteur/time,xminut,itt
      common/cou/cr(nr,nr),cn(nc,nc),ima(nr,nr),icn(nc,nc)
      common/rim/ckr(nr,nr,na,nlvl)
      common/coa/ckc(nr,nr,na,nlvl),jjmin
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/wat/qv(nlvl),qi(na,nlvl),qc(na,nlvl),xnc(na,nlvl),
     $      xni(na,nlvl),qce(na,nlvl),qie(na,nlvl)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/env/we(nlvl),thetae(nlvl),rhoqve(nlvl),fe(nr,nc,na,nlvl)
     & ,qv0(nlvl),te(nlvl),fie(nr,nc,na,nlvl),qve(nlvl),theta0(nlvl),
     $ wee(nlvl),f0e(nr,nc,1,nlvl)
      common/cyl/uhor(nlvl),bh(nlvl),bv(nlvl),rcyl1,rcyl2,alcar
      common/adv/w(nlvl),wm(nlvl),theta(nlvl),rhoqv(nlvl),zmean(nlvl)
      common/APP/rhos,xs,xnu,xphi,xex
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/supsate/rfde(nlvl),rfolde(nlvl),rfdie(nlvl),rfiolde(nlvl)
      common/oldy/told(nlvl),qvold(nlvl)
      common/oldye/tolde(nlvl),qvolde(nlvl)
      common/xcontact/xn40(nlvl),xnicont(nlvl)                 
      common/aterm/aterm1(nr,nc,na,nlvl),
     $ aterm2(nr,nc,na,nlvl),aterm3(nr,nc,na,nlvl)               
      common/wterm/wterm1(nlvl),wterm2(nlvl),wterm22(nlvl),wterm3(nlvl),
     $wterm4(nlvl),wterm5(nlvl)      
      common/press/uk(nf),uka(nf),ukb(nf),bj1a(nf),bj1b(nf),bjb0(nf)                                  
      common/press2/bja0(nf)
      common/at/at1(nlvl),at2(nlvl),at3(nlvl),at4(nlvl),at5(nlvl) 
      Common/sk3/dxx,axx,zjrs                              
      common/distr3/af(388,178,na,nlvl),afi(388,178,na,nlvl)   
      common/collision/ecc(nlvl,nr,nr),ecs(nlvl,nr,nr),
     $ e_beard(nlvl,nr,nr)
      common/speeds/vt(nlvl,nr),re_snow(nlvl,nr)
      common/multi/xx(nlvl),fdr1(nlvl,nr),fice1(nlvl,nr)
      common/gridab/nab(3),a(nr),b(nr),ab(nr),cpxx(nr)
       common/purb/purt(nlvl),pi0(nlvl),buoy(nlvl),buoye(nlvl),
     $ purt0(nlvl)
      common/graupel/D_graupel(nr),rho_graupel(nr),graupel_bin(nr)
      common/satuiw/ESW,QSW,ESI,QSI
      common/KIND/INA
      common/solubility/ebsilon
      common/fold1/fold(nr,nc,na,nlvl),feold(nr,nc,na,nlvl)
      COMMON/IV_terminal/Mhm 
      common/aconversion/Rauto(nlvl),Rautoa(nlvl),Raccr(nlvl),
     $       Raccra(nlvl)
      common/reference/n_in(nlvl),n_in_organic(nlvl),
     $       H_frac_organic(nlvl)
      integer icont,AGGR,L                                               
      real zmin(NLVL),xn10(na,nlvl),xdbz(nlvl),xmp(nr,na),xrd(nr)
      real xdbzi(nlvl)
      real xmeau(na),xndrop(na),xmglace(na),xncris(na),xn10mic(na)      
      real xmeaue(na),xndrope(na),xmglacee(na)           
      real  fia(nc,nlvl),cph11(nlvl),xnci(nlvl)
      real xcloud(na),xicloud(na),xrain(na),xsnow(na),xnni(na),xnnie(na)       
      real eta(nlvl),lambda(nlvl),xn_in_organic(nlvl),xn_in(nlvl)
      real alpha0,beta0,mass_test1,delta_mass,sh,xnino

      integer i,j,k 
      character *6 nom
      character *20 ifnam
      character *20 resnam
      character *3 myidout
      character *11 fmyidout


c----------------------------------------------------------------------------------
c--------------------------     INITIALISATION       ------------------------------

      ice=0        ! switch on off for ice formation
      INA=0        ! selection for the ice nuclei: 1:for dust,2:black carbon,3.organics (bacteria,pollen,leaf)
      ice_type=1
      ibac=1
      R0=.28705e7
      Rg=8.31e7
      Ra=2.8724e6
      pi=3.14159
      cc1=4./3.*pi
      t0=273.16
      xNa=6.023e23
      xk=1.3804e-16
      xh=6.6252e-27
      xMw=18.
      xl31cp=2834./1.005
      XL21CP=2477.4/1.005
      ebsilon = 0.1     !mass fraction of ammonium sulfate of IN
      xnu=2.            !VANT HOLF FACTOR
      xex=1.
      xphi=1.0
      xs=132.14        !atom mass of amonium safate
      rhos=1.841
      alcar=0.1
       k=1
       p(k)=800.
      qv(k)=1.10497E-02 
      t(k)=285.15

      num=0
      xni=0.              
      xnicont=0.                        
      xn40=0.                            
      nab(1)=10  !minimum grid of splinter
      nab(2)=77  !maximum grid of splinter ice
      nab(3)=nr  !maximum grid of graupel
      rhoaa=1.0  !for bacteria

      sjrs=50.
!2014.09.10      e(1)=2.e-18
      e(1)=4.e-15
      rhog=1.0
      rl(1)=(e(1)/cc1/rhog)**(1./3.)
      dx=alog(2.)/sjrs
      ax=2.**(1./sjrs)
    
      do j=2,nr
        e(j)=e(j-1)*ax
        rl(j)=(e(j)/cc1/rhog)**(1./3.)
        rdr=(e(j)/cc1)**(1./3.)
      IF(rl(j-1) .lt. 0.258 .and. rl(j) .ge. 0.258 ) Mhm = j
      enddo
      IF(rl(nr).lt.0.258) Mhm =nr
      do i=1,nc
         en(i)=e(i)
      enddo
!-------------------------------------------------for ice crystal a axis and b axis
      rhoice=0.9
      do j=1,nr
       cpxx(j)=rl(j)
       a(j)=2.0*rl(j)
       b(j)=2.0*rl(j)
      if(j.gt.nab(1).and.j.lt.nab(2)) then
       vv=e(j)/rhoice/
     $ (1.0-0.6*sin(real(j-nab(1))*pi/real(nab(3)-nab(1))))
       b(j)=(4.0*vv/(pi*(9.0*sin(real(j-nab(1))*pi/real(nab(2)-nab(1)))
     $ +1.0)))**(1./3.)
       a(j)=b(j)*(9.0*sin(real(j-nab(1))*pi/real(nab(2)-nab(1)))+1.0)
       Aab=(a(j)**2-b(j)**2)**(1.0/2.0)
       cpxx(j)=Aab/log((a(j)+Aab)/b(j))
      endif
      if(j.ge.nab(2).and.j.lt.nab(3)) then
       vv=e(j)/rhoice
       a(j)=2.*(vv*3./4.)**(1.0/3.0)
       b(j)=a(j)
       cpxx(j)=a(j)/2.
      endif
       ab(j)=max(a(j),b(j))
      enddo
!------------------------------------------------
          

	
!     2.2 Creation of vertical grids

      dt=2.0                 ! time step is in sec
      dz=100.               ! dz is in m  
      do i=1,NLVL
         zmin(i)=(i-1)*dz
         zmean(i)=zmin(i)+1.0*dz
      enddo

!      3. Initiation of the modelling case
 

      path=0.
      pathw=0.  
	
!      if (istop.eq.1) goto 2426
!       t(1)=285.15
!----------------------------------------------------------------------------------
      call iniap2
      xnci(k)=0
      do j=1,nr
      do i=1,nc
      
       xnci(k)=fi(j,i,1,k)+xnci(k)
           
      end do
      end do
!      print*,'xnci',xnci(k)
!--  --- -----------------------------------------------------------------------------
!      4. Calculation of initial radius of particle for collision
      do j=1,nr
         rdr=(e(j)/cc1)**(1./3.)
         jjmin=j
         if (rdr.gt.1.e-3) goto 9876
      enddo
 9876 continue
      do j=1,nr
         rdr=(e(j)/cc1)**(1./3.)
         jrain=j
         if (rdr.gt.4.e-3) goto 9877
      enddo
 9877 continue


      time=0.
      minu=0
      nom='out'
      write(ifnam,88)nom,minu

 88   format(a3,i3.3,'.dat')      
   
      Open (unit=30,file=ifnam,status='unknown',form='formatted')
      do j=1,nr
      rdr=(e(j)/cc1)**(1./3.) *1.e4
      write(30,*) rdr,f(j,1,1,1)/dx
      enddo
      close(30)
      



        nab(1)=10  !minimum grid of splinter
        nab(2)=77  !maximum grid of splinter ice
        nab(3)=nr  !maximum grid of ice pellets

        rhoice=0.9
      do j=1,nr
      cpxx(j)=rl(j)
      a(j)=2.0*rl(j)
      b(j)=2.0*rl(j)
      if(j.gt.nab(1).and.j.lt.nab(2)) then
      vv=e(j)/rhoice
     $ /(1.0-0.6*sin(real(j-nab(1))*pi/real(nab(3)-nab(1))))
      b(j)=(4.0*vv/(pi*(9.0*sin(real(j-nab(1))*pi/real(nab(2)-nab(1)))
     $+1.0)))**(1./3.)
      a(j)=b(j)*(9.0*sin(real(j-nab(1))*pi/real(nab(2)-nab(1)))+1.0)
      Aab=(a(j)**2-b(j)**2)**(1.0/2.0)
      cpxx(j)=Aab/log((a(j)+Aab)/b(j))
      endif
      if(j.ge.nab(2).and.j.lt.nab(3)) then
       vv=e(j)/rhoice
       a(j)=2.*(vv*3./4.)**(1.0/3.0)
       b(j)=a(j)
       cpxx(j)=a(j)/2.
      endif
      ab(j)=max(a(j),b(j))
      enddo
      
!         theta(1)=285.15
!          rhoqv(k)=0.12e-5
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!--------------------------     INTEGRATION     -----------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
        do itt=1,601         ! time steps for integration
         time=time+dt             ! time in second for simulation
         xminut=xminut+dt/60.     ! time in minutes
!          if(itt .eq.1) then
!         T(k)=theta(k)*(p(k)/1000.)**.286
!          end if
          theta(k)=t(k)*(1000./p(k))**.286
          rho(k)=1000.*p(k)/(R0/10000.)/T(k)
          rhoqv(k)=rho(k)*qv(k)
          call reynoldd(k,ice_type)
!          print*,'r',rho(k)
!         qv(k)=rhoqv(k)/rho(k)
            do l=1,na
            do i=1,nc
               do j=i,nr
                  if (f(j,i,l,k).lt.0.) f(j,i,l,k)=0.
                  if (fe(j,i,l,k).lt.0.) fe(j,i,l,k)=0.
                  if (fi(j,i,l,k).lt.0.) fi(j,i,l,k)=0.
                  if (fie(j,i,l,k).lt.0.) fie(j,i,l,k)=0.
                  if(f0(j,i,l,k).lt.0.) f0(j,i,l,k)=0.
                  if(f0e(j,i,1,k).lt.0.) f0e(j,i,1,k)=0.
               enddo
            enddo
         enddo
       
       
!         print*,'t',t(k),qv(k)
         If(pe(k).eq.0.0) pe(k)=p(k)
            Te(k)=thetae(k)*(pe(k)/1000.)**.286
!            qv(k)=rhoqv(k)/rho(k)
            qve(k)=rhoqve(k)/rho(k)

!----------------------------------------------------------------------------------
!------------------------------------------------
!--------------------------   MICROPHYSICS   -------------------------------------
!----------------------------------------------------------------------------------
       delp=0.0                        !sjm
       dz=100.
       delpe=0.0
!       do k=1,NLVL-3   !nlvl-3   to assure the buffer zone  K-loopK-loopK-loop
         k=1
         cph=0.
	 cphh=0.
         cphi=0.
	 cphnu=0.
         cphnuc=0.
         cphrim=0.
         cphm=0.0
         cphe=0.
         cphie=0.
         xn=0.
         xnn=0.
         rfdx=0.
	 Rauto(k) = 0.0
         Raccr(k) = 0.0
      
         Rautoa(k) = 0.0
         Raccra(k) = 0.0
       
	 call rfrhi(k)
!        write(*,*)t(k),'rfrhi',k,rfdi(k)
!          rfold(k)=rfd(k)
!C         if(k.ne.1) rfdx=rfd(k-1)
!----------------------------------------------------------------------------------
!     1. Microphysical processes
            xnci(k)=0.
          do i=1,nc
          do j=i,nr
            xnci(k)=fi(j,i,1,k)+xnci(k)
!           if (fi(j,i,1,k).ne.0.) PRINT*,'XNCI',j,XNCI(K)
          enddo
          enddo
        
          do l=1,na
            xndrop(l)=0.
            xmeau(l)=0.
            xndrope(l)=0.
            xmeaue(l)=0.            
          enddo
            xndropa=0.
         do l=1,na 
         do j=1,nr
            fist=0.0
            do i=1,nc
              xmeau(l)=f(j,i,l,k)*(e(j)-en(i))+xmeau(l)      !c equal to e(j)-e(i)
              xmeaue(l)=fe(j,i,l,k)*(e(j)-en(i))+xmeaue(l)
              xndrop(l)=f(j,i,l,k)+xndrop(l)
              xndrope(l)=fe(j,i,l,k)+xndrope(l)
            enddo
         enddo
         enddo

  111         format(1x,A6,1x,i2,5(1x,d13.5))     
!----------------------------------------       

        call depdrop(k,nr,cph,ibac)
  
         xL21=2.5E+6-(T(k)-t0)*2.3714E+03
         xL32=4.18E+03*(79.7+0.485*(T(k)-t0)-(2.5e-03)*((T(k)-t0)**2))

         xL31=xL32+xL21
         cp=1005.
         p0=1013.

50      format(2(1x,i4),4(1x,e15.6),1x,f6.1,1x,f8.3)

          if((rfd(k)-1.0).lt.1.e-6) cph=0.0        !just accounting the droplets evaporation because surface covered by liquid water
	  T11=T(k) 
          dd=cph*xl21/cp  
          T(K)=T(K)+ cph*xL21/cp        
         delp=delp-9.8*10000./(R0*t(k)*(1.+0.608*qv(k)))*dz  

          at5(k)=T(k)-T11

         write (82,51) zmean(k),T(k)

 51      format(f7.1,1x,f8.3)
         qv(K)=qv(k)-cph
                       
                                                                            
         call rfrhi(k)

         the1=theta(k)
         
        
 222     format(1x,4e12.4)
         rfold(k)=rfd(k)
         rfiold(k)=rfdi(k)
         told(k)=t(k)
         qvold(k)=qv(k)
c-------------------------------------
         rfolde(k)=rfde(k)
         rfiolde(k)=rfdie(k)
         tolde(k)=te(k)
         qvolde(k)=qve(k)
        

    
      if (xminut.ge.0.99) then        !sjm
!      if (xminut.ge.0.03) then        !sjm
      minu=1+minu
      nom='out'
      write(ifnam,88)nom,minu
      print*,'write dans out dat'     
       
      Open (unit=30,file=ifnam,status='unknown',form='formatted')
      do j=1,nr
      rdr=(e(j)/cc1)**(1./3.) *1.e4
      write(30,*) rdr,f(j,1,1,1)/dx
      enddo
      xminut=0.
      close(30)
      end if 
      
!      end if 
      end do
      
          s=0.
         do j=1,nr
      Open (unit=30,file='ice.txt',status='unknown',form='formatted')
       
     
        rdr=(e(j)/cc1)**(1./3.) *1.e4
         s=s+f(j,1,1,1)
        write(30,*) rdr,f(j,1,1,1)
!         write(30,*) f(j,1,1,1)
         enddo
!         print*,s
 
       close(30)
      
      print*,'end'
      END program


      subroutine depdrop(k,imax,cph,ibac)

      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      Common/maths/pi,cc1,t0,xh,xk
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      Common/sk2/dtd,dx,ax,sjrs
      Common/sk/dt,dz
      common/grid/en(nc),e(nr),rl(nr)
      common/distr/f(nr,nc,na,nlvl),fi(nr,nc,na,nlvl),f0(nr,nc,na,nlvl)
      common/switch/icoal,ice,irim,inuc,istop,jrain
      common/compteur/time,xminut,itt
      common/cou/cr(nr,nr),cn(nc,nc),ima(nr,nr),icn(nc,nc)
      common/rim/ckr(nr,nr,na,nlvl)
      common/coa/ckc(nr,nr,na,nlvl),jjmin
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/wat/qv(nlvl),qi(na,nlvl),qc(na,nlvl),xnc(na,nlvl),
     $      xni(na,nlvl),qce(na,nlvl),qie(na,nlvl)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/adv/w(nlvl),wm(nlvl),theta(nlvl),rhoqv(nlvl),zmean(nlvl)
      common/APP/rhos,xs,xnu,xphi,xex
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/env/we(nlvl),thetae(nlvl),rhoqve(nlvl),fe(nr,nc,na,nlvl)
     & ,qv0(nlvl),te(nlvl),fie(nr,nc,na,nlvl),qve(nlvl),theta0(nlvl),
     $ wee(nlvl),f0e(nr,nc,1,nlvl)
      common/KIND/INA
      real fold(nr,nc,na,nlvl),fist(nr),fd(nr),U(nr+1),qvdep
      integer iloop(nc)
      
      xmold=0.
      xfax=dt/dx
      tdep=t(k)
      qvdep=qv(k)
      open(91,file='m.dat',status='unknown')
      open(92,file='m1.dat',status='unknown')
      if (time.gt.2) then
       rf=(rfd(k)+rfold(k))*.5        
      else
       rf=rfd(k)
      endif
       LL=0
       rf= 1.001   !test                                  
       if(ibac.eq.0) ll=1
       do l=1,na-ll
       do i=1,nc
         do j=i,nr
        
         if(l.eq.1) xmold=f(j,i,l,k)*(e(j)-en(i))+xmold
         fold(j,i,l,k)=f(j,i,l,k)
         enddo
  
!         write(91,*)xmold,f(j,1,1,1)
        enddo
       enddo
!      print*,'e',f   
c-------------------------calculation for the activation size of bacteria
      sv=rf-1.
      if(INA .eq.3) then
       nactive=nc
      if(sv.gt.0.0) then
      xndrop=0.
      a=0.28
      b=0.93
      fa=0.
      xn=0.
        do i=1,nc
        do j=i,nr
               xndrop=f(j,i,1,k)+xndrop
         enddo
         enddo
        xnuo=a*xndrop*(sv*100.)**b
        do i=nc,1,-1
        do j=i,nr
          xn=xn+f(j,i,1,k)
         if(xnuo.lt.xn.and.xnuo.gt.0.0) then
          nactive=i-1
          goto 211
          endif
        enddo
          xx=xn
        enddo
211      continue
c--------------------------------------------------------------
       endif
       endif
 22   format(a8,3i3,3e12.5)
      do l=1,na-LL
        
      do i=1,nc
         iloop(i)=1
         umax=0.
         do j=i,nr-1
            call growth(j,i,l,k,rf)
c            print*,u5(1,1,2)
            ui=abs(u5(j,i,l))
            if (ui.gt.umax) umax=ui
         enddo
           u5(nr,i,l)=0.
         delta=dx/umax
csjm         ilom=IFIX(dt/delta+1.)
           ilom=int(dt/delta+1.)   
         fre=0
         if (ilom.gt.imax) then
            ilom=imax
    
csjm            u2max=FLOAT(ilom)/xfax
            u2max=1.0
            do kk=i,nr-1
         if (abs(u5(kk,i,l)).gt.u2max) u5(kk,i,l)=u2max*sign(1.,
     $     u5(kk,i,l))
      
               wmix=u5(kk+1,i,l)-u5(kk,i,l)
               if (wmix.gt.u2max) then
                u5(kk,i,l)=u2max*sign(.5,u5(kk,i,l))
                u5(kk+1,i,l)=u2max*sign(.5,u5(kk+1,i,l))
                  endif
                if(u5(kk+1,i,l).gt.0..and.u5(kk,i,l).lt.0) then
                 do jk=i,kk
                    u5(jk,i,l)=0.
                 enddo
                 endif
                 if(u5(kk,i,l).gt.0..and.u5(kk+1,i,l).lt.0) then
                 do jk=i,kk
                    u5(jk,i,l)=0.
                 enddo
                 endif
       
                 if(kk.eq.i) u5(kk,i,l)=0.
            enddo
         endif
         iloop(i)=ilom
         if(iloop(i).eq.0)iloop(i)=1
         if(ilom.lt.1) iloop(i)=1
!         print*,'iloop',iloop(i)
      enddo

         nn=1

      do 1111 i=nn,nc

         dtd=dt/iloop(i)
      do 701 it=1,iloop(i)
         U=0.0
       do 111 j=i,nr
111      U(j+1)=u5(j,i,l)*dtd/dx        ! Courant number 
         U(i)=0.
        fist=0.0
         do 222 j=i,nr
222         fist(j)=f(j,i,1,k)       
       call MPDATA1(U,fist,nr,2,1,1)
         do 333 j=i,nr
333       f(j,i,1,k)=fist(j)
 701  continue
1111  continue
      enddo !Na
      xm=0.

      do i=1,nc
         do j=i,nr
     
         xm=f(j,i,1,k)*(e(j)-en(i))+xm
          enddo
         enddo
      cph=xm-xmold
      cph=cph/rho(k)
!      print*,'cph',cph
      qvdep=qvdep-cph
      xL21=2.5E+6-(Tdep-t0)*2.3714E+03
      cp=1005.0
      dete=cph*xL21/cp
      Tdep=Tdep+dete
      P21=6.107*EXP(4028.*(Tdep-273.15)/(234.82*(Tdep-38.33)))
      QES=0.622*P21/(P(k)-0.378*P21)
!      print*,QES
      RF_new=Qvdep*(37.8*QES+62.2)/(QES*(0.378*Qvdep+0.622))
      RF_new=RF_new*.01
      if(rf .lt.1.0 .and. rf_new .ge. 1.0 .or. rf .ge. 1.0 .and.
     $ rf_new .lt. 1.0) then   
          do l=1,na-LL
          do j=1,nr
             do i=1,nc
         
             f(j,i,1,k)=fold(j,i,1,k)
              enddo
          enddo
          enddo
          call depdrop2(k,nr,cph,LL)
      endif
      
 93   format(1x,A30,I3,2F9.3)     
 988  format(1x,i4,3e12.5)
      return
      end




c----------------------------------------------------------------------------------
c	evolution of the spectrum of humid particles by condensation or evaporation
c case of the step of time divides by 10
      subroutine depdrop2(k,imax,cph2,LL)

      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      Common/maths/pi,cc1,t0,xh,xk
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      Common/sk2/dtd,dx,ax,sjrs
      Common/sk/dt,dz
      common/grid/en(nc),e(nr),rl(nr)
      common/distr/f(nr,nc,na,nlvl),fi(nr,nc,na,nlvl),f0(nr,nc,na,nlvl)
      common/switch/icoal,ice,irim,inuc,istop,jrain
      common/compteur/time,xminut,itt
      common/cou/cr(nr,nr),cn(nc,nc),ima(nr,nr),icn(nc,nc)
      common/rim/ckr(nr,nr,na,nlvl)
      common/coa/ckc(nr,nr,na,nlvl),jjmin
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/dep2/u52(nr,nc,na)
      common/wat/qv(nlvl),qi(na,nlvl),qc(na,nlvl),xnc(na,nlvl),
     $        xni(na,nlvl),qce(na,nlvl),qie(na,nlvl)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/adv/w(nlvl),wm(nlvl),theta(nlvl),rhoqv(nlvl),zmean(nlvl)
      common/APP/rhos,xs,xnu,xphi,xex
      common/oldy/told(nlvl),qvold(nlvl)
      common/oldye/tolde(nlvl),qvolde(nlvl)
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/KIND/INA
      real fist(nr),fd(nr),xl21,cp,fold(nr,nc,na,nlvl),U(nr+1)
      integer iloop(nc)
      delt=(t(k)-told(k))/dt
      delqv=(qv(k)-qvold(k))/dt
      dtold2=dt
      dt=dt/10.
      cph2=0.
      qvdep=qvold(k)
      tdep=told(k)
      pdep=p(k)
      do 120 ii=1,10
         qvdep=qvdep+delqv*dt
         tdep=tdep+delt*dt
         P21=6.107*EXP(4028.*(Tdep-273.15)/(234.82*(Tdep-38.33)))
         QES=0.622*P21/(P(k)-0.378*P21)
         RF=Qvdep*(37.8*QES+62.2)/(QES*(0.378*Qvdep+0.622))
         RF=RF*.01
         cph=0.
         xmold=0.
         xfax=dt/dx
         do i=1,NC
           do j=i,nr
             xmold=f(j,i,1,k)*(e(j)-en(i))+xmold
           do l=1,na-ll
             fold(j,i,1,k)=f(j,i,1,k)
           enddo
         enddo
       enddo
c-------------------------calculation for the activation size of bacteria
      sv=rf-1.
!      print*,'sv',sv
      if(INA .eq.3) then
       nactive=nc
      if(sv.gt.0.0) then
      xndrop=0.
      a=0.28
      b=0.93
      fa=0.
       xn=0.
        do i=1,nc
        do j=i,nr
               xndrop=f(j,i,1,k)+xndrop
        enddo
        enddo
        xnuo=a*xndrop*(sv*100.)**b
        do i=nc,1,-1
        do j=i,nr
           xn=xn+f(j,i,2,k)
          if(xnuo.lt.xn.and.xnuo.gt.xx) then
            nactive=i-1
            goto 211
          endif
        enddo
          xx=xn
        enddo
211      continue
       endif
       endif
c--------------------------------------------------------------
         do l=1,na-LL              !!!!!!!!!!!!!!!!!
         do i=1,nc
            iloop(i)=1
            umax=0.
           do j=i,nr-1
             call growth2(j,i,l,k,rf,tdep,pdep)
             
              ui=abs(u52(j,i,l))
             if (ui.gt.umax) umax=ui
              
           enddo
            u52(nr,i,l)=0.
           delta=dx/umax
           ilom=int(dt/delta+1.)   
           fre=0
           if (ilom.gt.imax) then
              ilom=imax
              u2max=FLOAT(ilom/2)/xfax
              do kk=i,nr-1
         if (abs(u52(kk,i,l)).gt.u2max) 
     $      u52(kk,i,l)=u2max*sign(1.,u52(kk,i,l))
          
               wmix=u52(kk+1,i,l)-u52(kk,i,l)
               if (wmix.gt.u2max) then
                u52(kk,i,l)=u2max*sign(.5,u52(kk,i,l))
                u52(kk+1,i,l)=u2max*sign(.5,u52(kk+1,i,l))
               endif 
              if(u52(kk+1,i,l).gt.0..and.u52(kk,i,l).lt.0) then
                 do jk=i,kk
                    u52(jk,i,l)=0.
                 enddo
              endif
              if(u52(kk,i,l).gt.0..and.u52(kk+1,i,l).lt.0) then
                 do jk=i,kk
                    u52(jk,i,l)=0.
                 enddo
              endif

                 if(kk.eq.i) u52(kk,i,l)=0.
              enddo
         endif
           iloop(i)=ilom
         if(iloop(i).eq.0)iloop(i)=1  
         if(ilom.lt.1) iloop(i)=1
        
        enddo

       nn=1
!      endif
      do 1111 i=nn,nc
         dtd=dt/iloop(i)
      do 701 it=1,iloop(i)
         U=0.0
       do 111 j=i,nr
111    U(j+1)=u52(j,i,l)*dtd/dx        ! Courant number
         U(i)=0.
        fist=0.0
         do 222 j=i,nr
222         fist(j)=f(j,i,1,k)
         call MPDATA1(U,fist,nr,2,1,1)
         do 333 j=i,nr
333       f(j,i,1,k)=fist(j)
 701  continue
1111  continue
      enddo !Na
      
       xm=0.
      do j=1,nr
         ii2=min0(j,nc)
         do i=1,ii2
            xm=f(j,i,1,k)*(e(j)-en(i))+xm
         enddo
      enddo
      cph=xm-xmold
      cph=cph/rho(k)
      cph2_old=cph2
      qvdep=qvdep-cph
      xL21=2.5E+6-(Tdep-t0)*2.3714E+03
      cp=1005
      dete=cph*xL21/cp
      Tdep=Tdep+dete
      P21=6.107*EXP(4028.*(Tdep-273.15)/(234.82*(Tdep-38.33)))
      QES=0.622*P21/(P(k)-0.378*P21)
      RF_new=Qvdep*(37.8*QES+62.2)/(QES*(0.378*Qvdep+0.622))
      RF_new=RF_new*.01
      cph2=cph2+cph
       if(rf .le.1.0 .and.rf_new .ge. 1.0 .or. rf .ge. 1.0 .and. 
     $ rf_new .lt. 1.0) then
       cph2=cph2_old
       do i=1,NC
           do j=i,nr
           do l=1,na-ll
             f(j,i,l,k)=fold(j,i,l,k)
           enddo
         enddo
       enddo
   
      exit
      endif  
 120  continue
      dt=dtold2
 988     format(1x,"d2 ",i4,2f8.4,2f9.5,f5.0)
      return
      end

c calculate the speed of growth of the drops dand the general case
      SUBROUTINE growth(j,i,l,k,rf)
      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      Common/maths/pi,cc1,t0,xh,xk
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      common/grid/en(nc),e(nr),rl(nr)
      Common/sk2/dtd,dx,ax,sjrs
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/compteur/time,xminut,itt
      common/APP/rhos,xs,xnu,xphi,xex
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/solubility/ebsilon
      ds=1.77     !density of ammonium sulfate
      e5=e(j)*2.**(.5/sjrs) ! for staggered grid in order to calculate the growth of droplet in the advection agrithom of Smolarkiewicz.(1996)
      wts=en(i)/e5
      call myhre(k,wts,rhosol,sigma)
      xmfpl0=6.6e-6
      rhow=1.
      P21=6.107*EXP(4028.*(T(k)-273.15)/(234.82*(T(k)-38.33)))
      XL21=2.5E+10-(T(k)-273.15)*2.3714E+07
      rdr=(e5/cc1/rhosol)**(1./3.)
      if (l .eq. 1) sol=xnu*xphi*xMw/xs*en(i)/(e5-en(i))
      if (l .eq. 2) sol=xnu*xphi*xMw/xs*ebsilon*en(i)/(e5-en(i))
      curv=2.*sigma*xMw/RG/T(k)/rdr/rhosol
      curva=2.*72.*xMw/RG/T(k)/rdr/rhow

      Ysol=3.3*1.0e-5/T(k)/rdr-sol!curv-sol
      DV1=.211*(T(k)/273.15)**1.94*1013./P(k)
      DVS=(2.*pi*xMw/T(k)/RG)**.5
      ALP=0.04
      xmfpl=xmfpl0*1013.25/293.15*T(k)/P(k)
      delnu=1.3*xmfpl
      DV=DV1/(rdr/(rdr+delnu)+DVS*DV1/ALP/rdr)
      XKA=(5.69+.017*(T(k)-273.15))*4.1868E+02
      xx=(vis(k)/DV)**(1./3.)*reynd(j,k)**(0.5)           
      if(xx.lt.1.4) AF=1.+0.108*xx**2.
      if(xx.ge.1.4) AF=0.78+0.308*xx                    
      AA=RG*T(k)*rhosol/1000./P21/DV/xMw
csjm      BB=XL21/T(k)/XKA*(xMw*XL21/RG/T(k)-1.)*exp(Ysol)
      BB=XL21*rhosol/T(k)/XKA*(xMw*XL21/RG/T(k)-1.)  
!      print*,rdr,rf,Ysol,AF,AA,BB
!      print*,AA+BB            
      if(l.eq.1) then
      u5(j,i,1)=4.*pi*rdr*(rf-exp(Ysol))*AF/(AA+BB)
      u5(j,i,1)=u5(j,i,1)/e5
      endif
      if (l.eq.2) then
      u5(j,i,2)=4.*pi*rdr*(rf-1.)*AF/(AA+BB)
      u5(j,i,2)=u5(j,i,2)/e5
      endif
      RETURN
      END
           


c----------------------------------------------------------------------------------
c recalculer la vitesse de croissance des gouttes si pas de tps divise par 10

      SUBROUTINE growth2(j,i,l,k,rf,tdep,pdep)
      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      Common/maths/pi,cc1,t0,xh,xk
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      common/grid/en(nc),e(nr),rl(nr)
      Common/sk2/dtd,dx,ax,sjrs
      common/APP/rhos,xs,xnu,xphi,xex
      common/dep2/u52(nr,nc,na)
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/solubility/ebsilon
      wts=en(i)/e(j)
      call myhre(k,wts,rhosol,sigma)
      e5=e(j)*2.**(.5/sjrs)
      xmfpl0=6.6e-6
      rhow=1.
      P21=6.107*EXP(4028.*(Tdep-273.15)/(234.82*(Tdep-38.33)))
      XL21=2.5E+10-(Tdep-273.15)*2.3714E+07
      rdr=(e5/cc1/rhosol)**(1./3.)
      if(l .eq.1) sol=xnu*xphi*xMw/xs*en(i)/(e5-en(i))
      if(l .eq.2) sol=xnu*xphi*xMw/xs*ebsilon*en(i)/(e5-en(i))
      curv=2.*sigma*xMw/RG/Tdep/rdr/rhosol
      curva=2.*72.*xMw/RG/Tdep/rdr/rhow

      Ysol=curv-sol
      DV1=.211*(Tdep/273.15)**1.94*1013./Pdep
      DVS=(2.*pi*xMw/Tdep/RG)**.5
      ALP=0.04
      xmfpl=xmfpl0*1013.25/293.15*Tdep/Pdep
      delnu=1.3*xmfpl
      DV=DV1/(rdr/(rdr+delnu)+DVS*DV1/ALP/rdr)
       xx=(vis(k)/DV)**(1./3.)*reynd(j,k)**(0.5)           
      if(xx.lt.1.4) AF=1.+0.108*xx**2.    
      if(xx.ge.1.4) AF=0.78+0.308*xx                    
      XKA=(5.69+.017*(Tdep-273.15))*4.1868E+02
      AA=RG*Tdep*rhosol/1000./P21/DV/xMw
      BB=XL21*rhosol/Tdep/XKA*(xMw*XL21/RG/Tdep-1.)                         
      if (l.eq.1) then
      u52(j,i,1)=4.*pi*rdr*(rf-exp(Ysol))*AF/(AA+BB)
      u52(j,i,1)=u52(j,i,1)/(e5) 
      endif
      if (l.eq.2) then
      u52(j,i,2)=4.*pi*rdr*(rf-1.)/(AA+BB)
      u52(j,i,2)=u52(j,i,2)/(e5)
      endif
      RETURN
      END
! c---------------------------------------------------------------------
c calcul de la densite et la tension de surface d'H2SO4, t et w donnes

      subroutine myhre(k,wts,rhosol,sigma)

      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
C     Calculation sulfate solution density from Myhre et al. (1998).
        w = wts
        C1      = T(k)-273.15
        C2      = C1**2
        C3      = C1**3
        C4      = C1**4
        A0 = 999.8426 + 334.5402e-4*C1 - 569.1304e-5*C2
        A1 = 547.2659 - 530.0445e-2*C1 + 118.7671e-4*C2
     $      + 599.0008e-6*C3
        A2 = 526.295e+1 + 372.0445e-1*C1 + 120.1909e-3*C2
     $      - 414.8594e-5*C3 + 119.7973e-7*C4
        A3 = -621.3958e+2 - 287.7670*C1 - 406.4638e-3*C2
     $      + 111.9488e-4*C3 + 360.7768e-7*C4
        A4 = 409.0293e+3 + 127.0854e+1*C1 + 326.9710e-3*C2
     $      - 137.7435e-4*C3 - 263.3585e-7*C4
        A5 = -159.6989e+4 - 306.2836e+1*C1 + 136.6499e-3*C2
     $      + 637.3031e-5*C3
        A6 = 385.7411e+4 + 408.3717e+1*C1 - 192.7785e-3*C2
        A7 = -580.8064e+4 - 284.4401e+1*C1
        A8 = 530.1976e+4 + 809.1053*C1
        A9 = -268.2616e+4
        A10 = 576.4288e+3
        den = A0 + w*A1 + w**2 * A2 + w**3 * A3 + w**4 * A4
        den = den + w**5 * A5 + w**6 * A6 + w**7 * A7
        den = den + w**8 * A8 + w**9 * A9 + w**10 * A10
        rhosol = den/1000.
        if (rhosol.lt.1.) rhosol=1.
      sula=142.35-0.96525*W*100.-T(k)*(0.22954-0.0033948*W*100.)
      sigma=sula
      return
      end
       subroutine reynoldd(k,ICE_TYPE)

      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      common/grid/en(nc),e(nr),rl(nr)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/graupel/D_graupel(nr),rho_graupel(nr),graupel_bin(nr)
      u0=1.68e-4     !g/cm/s
      t0=273.16
      s=110.5
      visckin=4.301e-02*T(k)**(2.5)/(p(k)*(T(k)+120.))
      vis(k)=u0*(t(k)/t0)**(1.5)*(t0+s)/(t(k)+s)/rho(k)
      call vdrop(k)
      do j=1,nr
!      If(ICE_TYPE .eq. 1) then 
        reynd(j,k)=2.*rl(j)*winf(j,k)/vis(k) 
!      else
!        reynd(j,k)=D_graupel(j)*winf(j,k)/vis(k)
!      endif
      enddo
      return
      end
 
     


!----------------------------------MPDATA
!        codes below are three different MPDATA routines. first goes 1-d, then 2-d on
!cartesian domain and finally 2-d on the sphere. the spherical routine is very
!  close to the cartesian one, except that it explicitly incorporates boundary
!  conditions characteristic for the problems on the sphere. MPDATA2 routine
!  can be used in arbitrary coordinates, as it admitts extra field H that may be
!  considered for the jacobian of coordinate transformation. in truly cartesian
!  coordinates H should be preset to ``1''. these routines were pulled from the
!  codes that deal with nontrivial dynamic applications (i.e., they work). some
! comments are included in MPDATA2 and they apply to all 3 routines. NOTE: both
!  codes MPDATA1 and MPDATA2 require some boundary routines (search for ``call b''
!  call to see example. these may be very simple in the spirit of
!  cyclic, e.g., x(1,j)=x(n-1,j), or zero-gradient, e.g., x(1,j)=x(2,j),
!  conditions. MPDATAS calls xbc which is provided and requires nothing else.
!  customized routines for variety of applications can be easily built from
!  current codes. 


      SUBROUTINE MPDATA1(U,X,M,IORD,ISOR,IFCT)
! THIS SUBROUTINE SOLVES 1-D ADVECTION IN CARTESIAN GEOMETRY ONLY
!SEE MPDATA2 FOR SOME COMMENTS
!      PARAMETER(NM=161+16,N=NM+1)
!      DIMENSION X(M),U(M+1),F(N),V(N),CP(NM),CN(NM)
!      REAL MX(NM),MN(NM)
       DIMENSION X(M),U(M+1),F(M+1),V(M+1),CP(M),CN(M)
      REAL MX(M),MN(M)

      DATA EP /1.E-15/
      DATA IDIV/0/
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1+AMIN1(0.,A)*Y2
      VDYF(X1,X2,A)=(ABS(A)-A**2)*(ABS(X2)-ABS(X1))
     $                           /(ABS(X2)+ABS(X1)+EP)
      VCOR3(X0,X1,X2,X3,A)= -A*(1.-3*ABS(A)+2.*A**2)
     $                       *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))
     $                        /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)/3.
      VCU(A1,A2,A3)=0.25*A2*(A3-A1)
      PP(Y)=AMAX1(0.,Y)
      PN(Y)=AMIN1(0.,Y)
!
       NM=M
      DO 1 I=1,M+1
    1 V(I)=U(I) 
!
      IF(IFCT.EQ.1) THEN
      DO 400 I=2,NM-1
      MX(I)=AMAX1(X(I-1),X(I),X(I+1))
  400 MN(I)=AMIN1(X(I-1),X(I),X(I+1))
      MX(1)=AMAX1(X(1),X(2))
      MX(NM)=AMAX1(X(NM-1),X(NM))
      MN(1)=AMIN1(X(1),X(2))
      MN(NM)=AMIN1(X(NM-1),X(NM))
      ENDIF
!    
!
                         DO 3 K=1,IORD
!
      DO 331 I=2,NM
  331 F(I)=DONOR(X(I-1),X(I),V(I))
      DO 333 I=2,NM-1
  333 X(I)=X(I)-(F(I+1)-F(I))
!      CALL BOUNDARY(X,NM,0)  
!      if (k.eq.1) then 
       x(1)=x(1)-f(2)
!        endif
       x(nm)=x(nm)+f(nm)
!        x(nm)=x(nm)
!        endif 
      IF(K.EQ.IORD) GO TO 6
      DO 501 I=2,NM
  501 F(I)=V(I)
      DO 51 I=2,NM
   51 V(I)=VDYF(X(I-1),X(I),F(I))
!
      IF(IDIV.EQ.1) THEN
      DO 50 I=3,NM-1
      V(I)=V(I)-VCU(F(I-1),F(I),F(I+1))
   50 CONTINUE
      ENDIF
!
      IF(ISOR.EQ.3) THEN
      DO 503 I=3,NM-1
  503 V(I)=V(I)+VCOR3(X(I-2),X(I-1),X(I),X(I+1),F(I))
      ENDIF 
      IF(IFCT.EQ.0) GO TO 3
!                   FLUX LIMITER
!
      DO 401 I=2,NM-1
      MX(I)=AMAX1(X(I-1),X(I),X(I+1),MX(I))
  401 MN(I)=AMIN1(X(I-1),X(I),X(I+1),MN(I))
      MX(1)=AMAX1(X(1),X(2),MX(1))
      MX(NM)=AMAX1(X(NM-1),X(NM),MX(NM))
      MN(1)=AMIN1(X(1),X(2),MN(1))
      MN(NM)=AMIN1(X(NM-1),X(NM),MX(NM))
      DO 399 I=2,NM
  399 F(I)=DONOR(X(I-1),X(I),V(I))
!
      DO 402 I=2,NM-1
      CN(I)=AMIN1(1., (X(I)-MN(I))/(PP(F(I+1))-PN(F(I))+EP))
  402 CP(I)=AMIN1(1., (MX(I)-X(I))/(-PN(F(I+1))+PP(F(I))+EP))
      DO 404 I=3,NM-1
      V(I)=PP(V(I))* 
     $ ( AMIN1(CP(I),CN(I-1))*PP(SIGN(1., X(I-1)))
     $  +AMIN1(CP(I-1),CN(I))*PP(SIGN(1.,-X(I-1))) )
     $    +PN(V(I))*
     $ ( AMIN1(CP(I-1),CN(I))*PP(SIGN(1., X(I  )))
     $  +AMIN1(CP(I),CN(I-1))*PP(SIGN(1.,-X(I  ))) )
  404 CONTINUE
      V(2)=PP(V(2))*(CP(2)*PP(SIGN(1.,X(1)))+CN(2)*PP(SIGN(1.,-X(1))))
     $   +PN(V(2))*(CN(2)*PP(SIGN(1.,X(2)))+CP(2)*PP(SIGN(1.,-X(2))))
      V(NM)=PP(V(NM))*
     $(CN(NM-1)*PP(SIGN(1.,X(NM-1)))+CP(NM-1)*PP(SIGN(1.,-X(NM-1))))
     $    +PN(V(NM))*
     $ (CP(NM-1)*PP(SIGN(1.,X(NM  )))+CN(NM-1)*PP(SIGN(1.,-X(NM  ))))
    3 CONTINUE
    6 CONTINUE
      RETURN
      END



      SUBROUTINE rfrhi(k)





      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      common/wat/qv(nlvl),qi(na,nlvl),qc(na,nlvl),xnc(na,nlvl),
     $      xni(na,nlvl),qce(na,nlvl),qie(na,nlvl)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/supsate/rfde(nlvl),rfolde(nlvl),rfdie(nlvl),rfiolde(nlvl)
      common/env/we(nlvl),thetae(nlvl),rhoqve(nlvl),fe(nr,nc,na,nlvl)
     & ,qv0(nlvl),te(nlvl),fie(nr,nc,na,nlvl),qve(nlvl),theta0(nlvl),
     $ wee(nlvl),f0e(nr,nc,1,nlvl)
      common/satuiw/ESW,QSW,ESI,QSI
!c       
	
      P21 = 6.107*EXP(4028.*(T(k)-273.15)/(234.82*(T(k)-38.33)))
      P31 = 6.1064*10**((9.5*(T(k)-273.15))/(265.5+(T(k)-273.15)))
      QES = 0.622*P21/(P(k)-0.378*P21)
!      print*,'QSS',QES
      RF =QV(k)/QES
      rfd(k) = rf
      rhi    = rf*p21/p31
      rfdi(k)= rhi
     
      ESW = P21
      QSW = QES
      ESI = P31
      QSI = 0.622*P31/(P(k)-0.378*P31)
!----------------------------------
     
!      write(*,*) QES
      return
      end
      
        subroutine iniap2

      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
      common/wat/qv(nlvl),qi(na,nlvl),qc(na,nlvl),xnc(na,nlvl),
     $      xni(na,nlvl),qce(na,nlvl),qie(na,nlvl)
      Common/maths/pi,cc1,t0,xh,xk
      common/APP/rhos,xs,xnu,xphi,xex
      Common/const/Rg,R0,xNa,xMw,xl31cp,xl21cp,ra
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
       Common/sk2/dt,dx,ax,sjrs
      common/adv/w(nlvl),wm(nlvl),theta(nlvl),rhoqv(nlvl),zmean(nlvl)
      common/env/we(nlvl),thetae(nlvl),rhoqve(nlvl),fe(nr,nc,na,nlvl)
     & ,qv0(nlvl),te(nlvl),fie(nr,nc,na,nlvl),qve(nlvl),theta0(nlvl),
     $ wee(nlvl),f0e(nr,nc,1,nlvl)
      common/distr/f(nr,nc,na,nlvl),fi(nr,nc,na,nlvl),f0(nr,nc,na,nlvl)
      common/grid/en(nc),e(nr),rl(nr)
      common/switch/icoal,ice,irim,inuc,istop,jrain
      common/supsat/rfd(nlvl),rfold(nlvl),rfdi(nlvl),rfiold(nlvl)
      common/supsate/rfde(nlvl),rfolde(nlvl),rfdie(nlvl),rfiolde(nlvl)
      common/em/em(nr)
      common/graupel/D_graupel(nr),rho_graupel(nr),graupel_bin(nr)
      common/KIND/INA
      common/solubility/ebsilon
      real N,LAMc,N0
      real fm0(nr),fm(nr),fma(nr),xmeau(nlvl),xmeaue(nlvl)
       yjrs=12.
       sjrs=50.
       nl=500
!2014.9.10      em(1)=1.e-19
      em(1)=4.e-15
      rhog=1.0
      pie=4./3.*pi
      rhoaa=1.0             ! for bacteria          
      ddx=alog(2.)/yjrs/3.0
      aax=2.**(1./yjrs)
      rv=461.5*18.01*1.0e4
       open(2,file='2.dat',status='unknown')
      do j=2,nr
        em(j)=em(j-1)*aax
        rdrr=(em(j)/cc1)**(1./3.)
      enddo
       open (unit=29,file='fi.dat',status='unknown',form='formatted')
        f(:,:,:,:)=0.
        fi(:,:,:,:)=0.
        f0(:,:,:,:)=0.
        fe(:,:,:,:)=0.
	 l=1                                     
         k=1
         delmd=0.
         xnd=0.
         eII=0.
         fm0=0.
         sigma=76.1-.155*(t(k)-273.15)
    
	  

            
          do i=1,nr-1
           rdr=(e(i)/pie)**(1./3.)
           
           rdrk1=(e(i+1)/pie)**(1./3.)
           dddx=alog(2.)/sjrs/3.
          
!            gamma distribution
            N=2.e8/1e6         !cloud droplet number concentration [m-3]
            alpha_c  =  1.    !shape parameter for cloud
            c=pi/6.
!sjm            Qcc=0.00002         !the mixing ratio of cloud
	    Qcc=0.0000016         !the mixing ratio of cloud
	    
            DE= 0.001!0.0006079242      !the density of air
            gam1=0.9999816       !sjm gam1=1
            gam2=23.9993127670914  !sjm gam2=24
            LAMc=(gam2*N*c/(gam1*DE*Qcc))**(1./3.)
	    N0=N*LAMc**(1+alpha_c)/gam1
!            fma(i)=N0*2*rdr*exp(-LAMc*2*rdr)*dddx/1e6
            !fma(i)=N0*4*rdr**2*exp(-LAMc*2*rdr)*dddx
            fma(i)=N*exp(-(log(rdr)-log(5.1e-5))**2/(2*2.16**2))*dddx
     $/2.16/sqrt(2.*pi)

	   enddo

         f(:,1,1,1)=0.   
	 
	 xxx = 0.0        
      
        do j=1,nr
        
              
                f(j,1,1,1)=fma(j)
     
        rdr=(e(j)/pie)**(1./3.) *1.e4
!        print*,rdr
            
        write(29,*) rdr,f(j,1,1,1)

         xxx = xxx + fma(j)
         enddo
        print*, "total cloud droplets=", xxx
	

       
      return
      end
          SUBROUTINE Vdrop(k)
ccc***** NOTE: 5 arrays: Dmcm(Mhm),velftur(Mhm),C_pT(Mhm),velftur_r(Mhm),D_eq(Mhm)     
ccc*******    should be described as arrays in main.for program ***********
ccc=============================================================================     
c            parameter (Mhm=30)                     !number of grid points by radius
      Parameter(nlvl=1,NC=1,NA=2,NR=2000)
c      parameter (Mhm=111)                     !number of grid points by radius 111:maximum diameter for water drops
      COMMON/IV_terminal/Mhm        
      Common/maths/pi,cc1,t0,xh,xk
      common/grid/en(nc),e(nr),rl(nr)
      common/veloc/u5(nr,nc,na),winf(nr,nlvl),
     &   vis(nlvl),winfi(nr,nlvl),ui5(nr),reyn(nr,nlvl),reynd(nr,nlvl)
      common/thermo/p(nlvl),pe(nlvl),t(nlvl),rho(nlvl),p0(nlvl)
        DATA alf/0.524/,bet/3./,gam/0.785/,sig/2./
        DATA del0KC/9.06/,cc0KC/0.29/                           !for drops

         dimension  xx(Mhm),                                       
     # aRe(Mhm), bRe(Mhm), avel(Mhm), bvel(Mhm), velf(Mhm),              
     #  Dmcm(Mhm),Rdcm(Mhm)                                          
     # ,aRetr(Mhm),bRetr(Mhm),D_bRe(Mhm),aksi(Mhm),bveltr(Mhm),
     #  aveltr(Mhm),velftur(Mhm),Psi_crs(Mhm),C_pT(Mhm),
     # vnonspher(Mhm),aveltr_r(Mhm),velftur_r(Mhm),D_eq(Mhm)   
ccc  31/I/2006



ccc================================================================================
cc                                                 
ccc=====================================
ccc===========================================================================================           
ccc                                INPUT PARAMETERS and 1 ARRAY:
ccc     1) alf;  2) bet (area-diameter relations, m=alf*D**bet);       
ccc     3) gam;   4) sig (mass-diameter relations, A=gam*D**sig)
ccc     5) CT=1.6 for turb. corr;   6) PDD=1, diameter in mcm; PDD=1.e3, diameter in mm
ccc     7) PVV=1, Vt in cm/s;     PVV=1.e2, Vt in m/s
ccc     8) pp is pressure in hPa;   9) TT is temperature in K, is needed for CpT correction
ccc    10) del0KC=9.06, drops; =5.83 crystals;   11) cc0KC=0.29, drops; =0.60 crystals;
ccc    12) Mhm is dimensions of arrays;          13) Rdcm(Mhm) - array of droplet radii in cm 
ccc    14) roa is air density (g/cm^3)           
ccc    ****************************************************************************************
ccc                                OUTUT 6 ARRAYS:
ccc    1) velf(Mhm) - Vt without turbul. corr;    2) velftur(Mhm) - Vt with turbul. corr;
ccc    3) C_pT(Mhm) - T-p correction (height);    4) vnonspher(Mhm) - param. nonsphericity
ccc    5) Dmcm(Mhm) - array of diameters, mcm;    6) D_eq(Mhm) - equivalent diameter in cm
ccc******************************************************************************************** 

cc 24 Jan 2006           roa_00=1.2e-3                   !air density at the surface, 0 km, 1000 mb
          roa_00=rho(k)
          roa=rho(k)
          PVV=1.
c            anuvisc=0.11                 !kinem. viscos=0.11 cm^2/s
          TT=T(k)
          PP=P(k)
          TTC=T(k)-273.15
      if (TTC.le.0.) aeta_vis=(1.718+0.0049*TTC-1.2e-5*TTC**2)*1.e-4   !dynam. viscos. air, PK97, p. 417
      if (TTC.gt.0.) aeta_vis=(1.718)*1.e-4
c           anuvisc0=0.15                 !kinem. viscos=0.15 cm^2/s
        anuvisc0=aeta_vis/roa_00            !kinem. viscos
         del24KC=del0KC**2/4.                      !Delta**2/4
          pips3= del0KC**2*cc0KC**0.5/4.           ! 8 July 2004 
c       cc1KC=4./(del0KC**2*cc0KC**05)
        cc1KC=1./pips3                             !coeff. C1
c       write (*,*) ' cc1KC=',cc1KC,' anuvisc0=',anuvisc0,
c     #  ' roa_00=',roa_00
c         pause 66        
c         write (*,*) 'Mhm=',Mhm,' Rdcm(kmh)=',(Rdcm(kmi),kmi=1,Mhm)
c          pause 33   

         do 777 kmh=1,Mhm
c           write (*,*) 'kmh=',kmh
           Rdcm(kmh)=rl(kmh)

          Dmcm(kmh)=2*Rdcm(kmh)*1.e4                           !current diameter Dmcm in mcm
ccc           Dmcm(kmh)=1*Rdcm(kmh)*1.e4                !test 7 Feb 2007; radius, current diameter Dmcm in cm
c        write (*,*) 'Rdcm(kmh)=',Rdcm(kmh),' roa=',roa,' anuvisc0=',
c     # anuvisc0,' Dmcm(kmh)=',Dmcm(kmh),' alf=',alf,' bet=',bet,
c     # ' gam=',gam,' sig=',sig

cc         DD00=0.37        !in cm
cc         DD00=0.47        !in cm
c         DD00=0.57        !in cm
c         DD00=0.17        !in cm
cc         DD22=0.37
cc         DD22=0.47

c        tulip=(1./(1+(dtek(KR)/1.e4)/DDD0))**(+1.0)   !dtek in cm, dtek<8.5 mm
c26 Jan 2006       tulip=(1./(1+((Dmcm(kmh)/1.e4)/DDD0))**2.)**(+0.5)   !dtek in cm, dtek<8.5 mm
c26 Jan 2006       tulip=(1./(1+((Dmcm(kmh)/1.e4)/DDD0))**5.)**(+0.2)   !dtek in cm, dtek<8.5 mm
         div=Dmcm(kmh)/1.e4                  !div=diameter in cm, for non-sphericity 26 Jan 2006
          tulip=1.
cc      tulip=exp(-div/DD22)+(1.-exp(-div/DD22))*(1./(1+div/DD00))**(+1.0)
c          if (tulip.lt.0.55) tulip=0.55

!************  START 7 June 2007, nonsphericity from Testik, Adv. Geophys., 2007 ********
!***              div=diameter in cm
          piv=Dmcm(kmh)/1.e3                                !piv=diameter in mm
!**********    Andsager et al. [1999]   ************************        
c       if (piv.ge.1.1.and.piv.le.4.4) tulip=                     ! Andsager et al. [1999]
c     #        1.012 - 0.144*div -1.03*div**2
        if (piv.ge.1.1) tulip=1.012 - 0.144*div -1.03*div**2         ! Andsager et al. [1999]
!*********************  Beard and Chuang [1987]  **********************************************        
c good!       if (piv.ge.1.1.and.piv.le.4.4) tulip=                         ! Beard and Chuang [1987]
c     # 1.0048 + 0.0057*div - 2.628*div**2 +3.682*div**3 -1.677*div**4
           
c        if (piv.ge.1.1) tulip=                         ! Beard and Chuang [1987]
c     # 1.0048 + 0.0057*div - 2.628*div**2 +3.682*div**3 -1.677*div**4
     
!************  END 7 June 2007, nonsphericity from Testik, Adv. Geophys., 2007 ********
       vnonspher(kmh)=tulip
c        if (Dmcm(kmh).lt.2000.) tulip=1.     !???  26 Feb 2000, corr. for nonspheric.
c        D_eq(kmh)=Dmcm(kmh)*tulip**0.333333/1.e4          !equavalent dimater in cm
        D_eq(kmh)=Dmcm(kmh)*tulip**0.333333                !equavalent dimater in mcm

c       pause 777     
!************** END for drop non-sphericity **************************
        bims=2.*alf*980.
     #     *vnonspher(kmh)                       !for non-sphericity
        boms=(roa*anuvisc0**2*gam)
        bams=(Dmcm(kmh)/1.e4)**(bet+2.-sig) 
c      write (*,*) 'kmh=',kmh,' bims=',bims,' boms=',boms,' bams=',bams,
c     # ' tulip=',tulip 
c         pause 22         
        xx(kmh)=bims/boms*bams                                  !Best param. xx via Dmcm(kmh        
c        xx(kmh)=(2.*alf*980./(roa*anuvisc0**2*gam))*           !Best param. xx via Dmcm(kmh
c     # (Dmcm(kmh)/1.e4)**(bet+2.-sig) 
       
!********************************************************************
          XXtek=xx(kmh)
!***========= calculation of KC02 modified coeffic. ====================       
c  new 19 March 2003, del0KC=5.83 ;  cc0KC=0.60;  del24KC=del0KC**2/4.;  cc1KC=4./(del0KC**2*cc0KC**05)
            braketKC= (1.+cc1KC*XXtek**0.5)**0.5 
        bre(kmh)=cc1KC*XXtek**0.5/2./(braketKC-1)/braketKC                       !coeffic. b1 smooth
        are(kmh)=del24KC*(braketKC-1)**2/XXtek**bre(kmh)
c 24 Jan 2006         Re_x_KC(kmh)= are(kmh)*XXtek**bre(kmh)
c 24 Jan 2006         CD_KC(kmh)=XXtek/(Re_x_KC(kmh))**2                
!**********  8 July 2004, clculation of avel and bvel, eq. (2.11), (2.12) in KC02 ================
c alfagr=0.00145  betagr=1.8  gamagr=0.2285  sigagr=1.88  roa_00=1.2e-3  anuvisc0=0.15
        bvel(kmh)=bre(kmh)*(bet+2.-sig)-1.                        !coef. veloc. avel,  good!
        anum1= 2.*alf*980.                                        ! = 2.842
     #     *vnonspher(kmh)                                                 !for nonspher. 1 Feb 2006
        denom2=roa_00*anuvisc0**2*gam
        fract1=anum1/denom2
          avel(kmh)=are(kmh)*anuvisc0*(fract1)**bre(kmh)              !coef. veloc. bvel,
cccc            !to convert from microns to mm,  16 Sep 2004 
       velf(kmh)=avel(kmh)*(Dmcm(kmh)/1.e4)**bvel(kmh)              !velocity itself velf, no turb. corr
       velf(kmh)=velf(kmh)/PVV                            !scaled velocity in m/s or cm/s, no turb. corr  
c if Dmcm in mm        velf(kmh)=avel(kmh)*(Dmcm(kmh)/1.e1)**bvel(kmh)              !velocity itself velf,        
 777      continue              !end cycle over XX and D     
!**************************************************************
cc
!*******  1 October 2004, turbulent correction         
            akturcr=1.                                 !18 Sep 2004, for crystals
c           akturcr=1.5                                !18 Sep 2004, for crystals
c           akturcr=2.                                 !2 Oct 2004, as in Bohm'92
c         X0cr=2.8e6                     !from Bohm'92 for crystals in turb.,  KC05, after eq. (3.2)
         X0cr=6.7e6                     !from Bohm'92 for drops in turb., KC05, after eq. (3.2), p.4345
c          CT=1.6                       !turb. cor. for CD
           CT=1.6 !SJM
c 24 Jan 2006         TTCC=TT-268.15              !temperature C relative to -5 C
         TTCC=TT-273.15              !temperature C relative to -5 C
         fi_T=1.+2.85e-3*TTCC-6.9e-6*TTCC**2            !correct. from viscosity, PK97, p. 417, ch. 10
c         write (*,*) ' fi_T=',fi_T
c         pause 33
         pp0=1000.                                     !pressure 1000 hPa
         TT0=273.15                                    !T = 0 K
c 24 Jan 2006        TT0=268.15         
!********  START cycle of recalculation with turbul. corr., by KC05 ====================
           do 888 kmh=1,Mhm
             XXtek=xx(kmh)                   ! xx(kmh) calculated in previous cycle
           zz=XXtek/X0cr
           zzkk=zz**akturcr
          Psi_crs(kmh)=(1.+zzkk)/(1.+CT*zzkk)                          !Bohm function for turb. corr.
        D_bRe(kmh)=-akturcr/2.*(CT-1.)*zzkk/(1.+CT*zzkk)/(1.+zzkk)     !turb. cor. for bRe, KC05, eq.(3.5)
        bRetr(kmh)=bRe(kmh) + D_bRe(kmh)                               !bRe with turb. cor., KC05, eq.(3.4)
        aksi(kmh)= ((Psi_crs(kmh))**0.5)/XXtek**D_bRe(kmh)             !correction to aRe, , KC05, eq.(3.7)
        aRetr(kmh)=aRe(kmh)*aksi(kmh)                                  !corrected aRe, , KC05, eq.(3.7)
cc      new:     D_bRe(kmh), aksi(kmh),bveltr(kmh),aveltr(kmh),velftur(kmh)
        bveltr(kmh)=bretr(kmh)*(bet+2.-sig)-1.                        !coef. veloc. bvel,  good!
        anum1= 2.*alf*980.                                             ! = 2.842
     #     *vnonspher(kmh)                                                 !for nonspher. 1 Feb 2006          
        denom2=roa_00*anuvisc0**2*gam
        fract1=anum1/denom2
          aveltr(kmh)=aretr(kmh)*anuvisc0*(fract1)**bretr(kmh)          !coef. veloc. avel,
        aveltr_r(kmh)=aveltr(kmh)*2**bveltr(kmh)                        !Av by radii, 7 Feb 2007
cccc            !to convert from microns to mm,  16 Sep 2004 
          velftur(kmh)=aveltr(kmh)*(Dmcm(kmh)/1.e4)**bveltr(kmh)        !Vt(D) with turb. corr, KC05
          velftur_r(kmh)=velftur(kmh)*2**bveltr(kmh)                    !Vt(r) by radii, 7 Feb 2007
        velftur_r(kmh)=aveltr_r(kmh)*(Dmcm(kmh)/1.e4/2.)**bveltr(kmh)        !Vt(D) with turb. corr, KC05

        pT_rat=((pp0/pp)*(TT/TT0))  
c      C_pT(kmh)=fi_T*pT_rat**(1-bRetr(kmh))                           !p-T correction C_pT
       C_pT(kmh)=(fi_T)**(1-2*bRetr(kmh))*pT_rat**(1-bRetr(kmh))      !p-T correction C_pT
                   
       velftur(kmh)=velftur(kmh)/PVV               !scaled veloc., m/s or cm/s, with turb. corr
cc 21 Nov 2004         *C_pT(kmh)/PVV         
!SJM       Dmcm(kmh)=Dmcm(kmh)/PDD                                !D in mm (PDD=1.e3) or cm (PDD=1.e4)
 888      continue
       
       do i=1,nr
       if(i.le.Mhm) winf(i,k)=velftur_r(i)
       if(i.gt.Mhm) winf(i,k)=velftur_r(Mhm)
       enddo 
ccc====================================================================================================      
       return 
       end
       

