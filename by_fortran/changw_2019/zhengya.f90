!      ��ѹԭʼ����ģʽ
!      1985 10 29 by  Shen tongli
!      m=20 Ϊx����������n=16 Ϊy����������dΪ����࣬rmΪ�Ŵ�ϵ��
!      fΪ��ת������wΪ�������飬cla��clo�ֱ�Ϊ��������γ�Ⱥ;���
!      dtΪʱ�䲽����sΪƽ��ϵ��
!      ua��ub��uc�ֱ�Ϊn-1��n��n+1ʱ����x�������
!      va��vb��vc�ֱ�Ϊn-1��n��n+1ʱ����y�������
!      za��zb��zc�ֱ�Ϊn-1��n��n+1ʱ����λ�Ƹ߶�
!      na���ڿ���12Сʱ��Ԥ����nb���ڼ�¼ʱ����ֲ�����nt2=72�����б�
!      �Ƿ����12Сʱ���Ƿ�����ڵ�ƽ����nt4=6�����ж��Ƿ�����߽�ƽ����
!      nt5�����ж��Ƿ����ʱ��ƽ���� 
!      zo��Ϊ�˼�С���������Ⲩ�Ĳ��٣����Ӳ�ָ�ʽ���ȶ��Զ������λ�Ƹ߶ȡ�
      program shen2
	  parameter(m=20,n=16,d=300000.0,cla=51.0,clo=118.0,dt=600.0)
      dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),&
      uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n),w(m,n)
      zo=2500.0
      s=0.5
      nt2=72
      nt4=6
      nt5=36
      c1=dt/2.0
      c2=dt*2.0
!  Ϊ����Grads��ͼ��������λ�Ƹ߶ȳ������ļ�h.grd��������ʼ����Ԥ������    
      open(10,file='D:\YEWEI\shuzhi\h.grd',form='binary')         

!  ����Ŵ�ϵ���͵�ת����,��д�������ļ���
      call cmf(rm,f,d,cla,m,n)
      open(1,file='D:\YEWEI\shuzhi\rm.dat')
       write(1,101) rm
      close(1)
 101  format(20f10.5)
      open(1,file='D:\YEWEI\shuzhi\f.dat')
       write(1,103) f
	close(1)
 103  format(20e15.5)

!  �����ʼ���ϳ� 
      open(2,file='D:\YEWEI\shuzhi\za.dat',status='old')
  	 read(2,102) za 
      close(2)
 102  format(20f6.0)

!  Ϊ������ͼ,����ʼ��������д������������ļ�h.grd��
 	write(10) ((za(i,j),i=1,m),j=1,n)     !�����ǰ�֮ǰ����dat�е�zaд��grd

!ccccccccccccccccccccccc�����ת���ӳ����˴���Ҫ�޸�cccccccccccccccccccccccccc
!  �����ת���ֵ   
      call cgw(ua,va,za,rm,f,d,m,n)  !�����ת���ֵ
      open(4,file='D:\YEWEI\shuzhi\ua.dat')
      write(4,104) ua  !д���ת���ֵ
      close(4)
      open(5,file='D:\YEWEI\shuzhi\va.dat')
      write(5,104) va  !д���ת���ֵ
      close(5)
 
 104  format(20f10.5)

!   ��ֵ�����ӳ���
      call tbv(ub,vb,zb,ua,va,za,m,n)
      call tbv(uc,vc,zc,ua,va,za,m,n)

!   ��ʼԤ��  
      do 10 na=1,2
       nb=0
!   ŷ��������1Сʱ
      do 20 nn=1,6
       call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,dt,zo,m,n)
       call ti(ua,va,za,ub,vb,zb,ua,va,za,rm,f,d,dt,zo,m,n)
       nb=nb+1
  20  continue

!   �߽�ƽ���ӳ���
      call ssbp(za,w,s,m,n)
      call ssbp(ua,w,s,m,n)
      call ssbp(va,w,s,m,n)

!   ǰ����ְ벽
      call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,c1,zo,m,n)
!   �������ְ벽
      call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
      nb=nb+1
!  ���鴫���ӳ���
      call ta(ub,vb,zb,uc,vc,zc,m,n)
!  ��������һ��,������11Сʱ
      do 30 nn=1,66
       call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,c2,zo,m,n)
       nb=nb+1
!  ��ӡ���ֲ���,na��ѭ����,nbСѭ����
       call pv(na,nb)
!  �ж��Ƿ����12Сʱ
       if(nb.eq.nt2) go to 80
!  �ж��Ƿ����߽�ƽ��
       if(nb/nt4*nt4.eq.nb) go to 40
       go to 50
  40   call ssbp(zc,w,s,m,n)
       call ssbp(uc,w,s,m,n)
       call ssbp(vc,w,s,m,n)

!  �ж��Ƿ���ʱ��ƽ��
  50   if(nb.eq.nt5) go to 60
       if(nb.eq.nt5+1) go to 60
       go to 70
!  ʱ��ƽ���ӳ���
  60   call ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,m,n)

!  ���鴫��,Ϊ��һ�ֻ�����׼��
  70   call ta(ua,va,za,ub,vb,zb,m,n)
       call ta(ub,vb,zb,uc,vc,zc,m,n)
  30  continue

!  �����ڵ�ƽ��
  80  call ssip(zc,w,s,m,n,2)
      call ssip(uc,w,s,m,n,2)
      call ssip(vc,w,s,m,n,2)

!  ��ӡ���ֲ���
      call pv(na,nb)

!  ���鴫��,Ϊ��12Сʱ�Ļ�����׼��
  10  call ta(ua,va,za,uc,vc,zc,m,n)

!  ���Ԥ�����
    open(6,file='D:\YEWEI\shuzhi\zc.dat')
    write(6,102) zc
    close(6)
	write(10) ((zc(i,j),i=1,m),j=1,n)!��Ԥ���ĸ߶ȳ�ֵд��ԭ�ȵ�grd�ļ�
    close(10)    
	open(7,file='D:\YEWEI\shuzhi\uc.dat')
      write(7,104) uc
      close(7)
      open(8,file='D:\YEWEI\shuzhi\vc.dat')
      write(8,104) vc
      close(8)
      open(13,file='D:\YEWEI\shuzhi\uv.grd',form='binary')
      write(13) uc
      write(13) vc
      close(13)
      end

!     �����ͼ�Ŵ�ϵ���Ϳ��ϲ���
!     rkΪԲ׶����,rlqΪ������ͶӰӳ��ƽ���ϳ����������ľ���,aΪ����뾶
!     sitaΪ��׼��γ,psxΪ����������γ,rΪģʽ���ĵ������ľ���
      subroutine cmf(rm,f,d,cla,m,n)
      dimension rm(m,n),f(m,n)
      rk=0.7156
      rlq=11423370.0
      a=6371000.0
      conv=57.29578
      w1=2.0/rk
      sita=30.0/conv
      psx=(90.0-cla)/conv

!  ����ģʽ���ĵ������ľ���r 
      cel0=a*sin(sita)/rk
      cel=(tan(psx/2.0))/(tan(sita/2.0))
      r=cel0*cel**rk

!  ȷ����������ԭ���ڵ�ͼ����ϵ�е�λ��
      xi0=-(m-1)/2.0
      yj0=r/d+(n-1)/2.0

!   ��������������ľ���rl,(xj,yi)Ϊģʽ������ڵ�ͼ����ϵ�е�λ��  
      do 10 i=1,m
      do 10 j=1,n
      xi=xi0+(i-1)
      yj=yj0-(j-1)
      rl=sqrt(xi**2+yj**2)*d

!   ��Ŵ�ϵ��rm�Ϳ��ϲ���f
      w2=(rl/rlq)**w1
      sinl=(1.0-w2)/(1.0+w2)
      rm(i,j)=rk*rl/(a*sqrt(1.0-sinl**2))
  10  f(i,j)=1.4584e-4*sinl
      return
      end

!    �����ת��
!    ��ͬѧ��д��ת���ֵ���ӳ��򣡣���Ӧ�����У�4.134��ʽ��Ҫ��:���������ĵط�һ��Ҫ�������!
     subroutine cgw(ua,va,za,rm,f,d,m,n) 
     dimension ua(m,n),va(m,n),za(m,n),rm(m,n),f(m,n)
     m1=m-1
     n1=n-1
     do i=1,m
       do j=2,n1
          ua(i,j)=(-rm(i,j))*9.8/f(i,j)*((za(i,j+1)-za(i,j))/d) 
          !����������ʽ����γ������ڵ��ֵ
        end do
     end do
     do j=1,n
        do i=2,m1
        va(i,j)=rm(i,j)*9.8/f(i,j)*((za(i+1,j)-za(i,j))/d)
        !����������ʽ���㾭������ڵ��ֵ
        end do
     end do 
     do i=1,m
     ua(i,1)=(-rm(i,1))*9.8/f(i,1)*((za(i,2)-za(i,1))/d) !��ǰ���ʽ���߽�ֵ
     ua(i,n)=(-rm(i,n))*9.8/f(i,n)*((za(i,n)-za(i,n-1))/d)
     end do 
     do j=1,n
     va(1,j)=rm(1,j)*9.8/f(1,j)*((za(2,j)-za(1,j))/d) !��ǰ���ʽ���߽�ֵ
     va(m,j)=rm(m,j)*9.8/f(m,j)*((za(m,j)-za(m-1,j))/d)
     end do
     end subroutine cgw



!     ʱ������
      subroutine ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
      dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),&
      uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n)
       c=0.25/d
       m1=m-1
       n1=n-1
      do 10 i=2,m1
      do 10 j=2,n1
      e=-c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(ub(i+1,j)-ub(i,j))&
      +(ub(i,j)+ub(i-1,j))*(ub(i,j)-ub(i-1,j))&
      +(vb(I,j-1)+vb(i,j))*(ub(i,j)-ub(i,j-1))&
      +(vb(I,j)+vb(i,j+1))*(ub(i,j+1)-ub(i,j))&
      +19.6*(zb(i+1,j)-zb(i-1,j)))+f(i,j)*vb(i,j)
      uc(i,j)=ua(i,j)+e*dt
      g=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(vb(i+1,j)-vb(i,j))&
      +(ub(I,j)+ub(i-1,j))*(vb(i,j)-vb(i-1,j))&
      +(vb(I,j-1)+vb(i,j))*(vb(i,j)-vb(i,j-1))&
      +(vb(I,j)+vb(i,j+1))*(vb(i,j+1)-vb(i,j))&
      +19.6*(zb(i,j+1)-zb(i,j-1)))-f(i,j)*ub(i,j)
  10  vc(i,j)=va(i,j)+g*dt
      do 20 i=2,m1
      do 20 j=2,n1
      h=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(zb(i+1,j)-zb(i,j))&
      +(ub(I,j)+ub(i-1,j))*(zb(i,j)-zb(i-1,j))&
      +(vb(I,j-1)+vb(i,j))*(zb(i,j)-zb(i,j-1))&
      +(vb(I,j)+vb(i,j+1))*(zb(i,j+1)-zb(i,j))&
      +2.0*(zb(i,j)-zo)*(ub(i+1,j)-ub(i-1,j)+vb(i,j+1)-vb(i,j-1)))
  20  zc(i,j)=za(i,j)+h*dt
      return
      end

!     ʱ��ƽ��
      subroutine ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,m,n)
      dimension ua(m,n),va(m,n),za(m,n),&
      ub(m,n),vb(m,n),zb(m,n),&
      uc(m,n),vc(m,n),zc(m,n)
      m1=m-1
      n1=n-1
      do 10 i=2,m1
      do 10 j=2,n1
      ub(i,j)=ub(i,j)+s*(ua(i,j)+uc(i,j)-2.0*ub(i,j))/2.0
      vb(i,j)=vb(i,j)+s*(va(i,j)+vc(i,j)-2.0*vb(i,j))/2.0
  10  zb(i,j)=zb(i,j)+s*(za(i,j)+zc(i,j)-2.0*zb(i,j))/2.0
      return
      end

!   ������5��ƽ��(����ƽ��)
!   ��ͬѧ��д������5��ƽ��(����ƽ��)���ӳ��򣡣���Ӧ�����У�4.126��ʽ
!   ע���˳��������Ƴɿ�����ʽ����֤�ȿ�ѡ������ƽ�����ֿ�ѡ����ƽ��   l=1Ϊִֻ����ƽ����l=2Ϊִ������ƽ��.
   subroutine ssip(a,w,s,m,n,l)
!  call ssip(zc,w,s,m,n,2)�Ǽ���߶ȳ�ʱcall���ӳ���ʱ�ķ�ʽ
!  wΪ�������飬s��ƽ��ϵ����m=20 Ϊx����������n=16 Ϊy��������
   dimension a(m,n),w(m,n)
   m1=m-1
   n1=n-1
   if (l==1) then !���أ�l=1Ϊ��ƽ��
   do  i=2,m1
      do  j=2,n1
          w(i,j)=a(i,j)+s*0.25*(a(i+1,j)+a(i,j+1)+a(i-1,j)+a(i,j-1)-4*a(i,j))
      end do
   end do
   else     !��������ƽ��
   do  i=2,m1
       do  j=2,n1
       w(i,j)=a(i,j)+s*0.25*(a(i+1,j)+a(i,j+1)+a(i-1,j)+a(i,j-1)-4*a(i,j))
       end do
   end do
   do i=2,m1
      do j=2,n1
         a(i,j)=w(i,j)
      end do
   end do
    do  i=2,m1
       do  j=2,n1
       w(i,j)=a(i,j)-s*0.25*(a(i+1,j)+a(i,j+1)+a(i-1,j)+a(i,j-1)-4*a(i,j))
       end do
   end do
   end if
   do i=2,m1
      do j=2,n1
         a(i,j)=w(i,j)
      end do
   end do
   end subroutine ssip






!   �߽�ŵ�ƽ��
      subroutine ssbp(a,w,s,m,n)
      dimension a(m,n),w(m,n)
      m1=m-1
      m3=m-3
      n1=n-1
      n2=n-2
      n3=n-3
      do 10 i=2,m1
      do 10 j=2,n1,n3
  10  w(i,j)=a(i,j)+0.5*s*(1.0-s)*&
      (a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))+0.25*s*s*&
      (a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
      do 20 i=2,m1,m3
      do 20 j=3,n2
  20  w(i,j)=a(i,j)+0.5*s*(1.0-s)*&
     (a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))+0.25*s*s*&
     (a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
      do 30 i=2,m1
      do 30 j=2,n1,n3
  30  a(i,j)=w(i,j)
      do 40 i=2,m1,m3
      do 40 j=3,n2
  40  a(i,j)=w(i,j)
      return
      end


!     ���鴫��
      subroutine ta(ua,va,za,ub,vb,zb,m,n)
      dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n)
      do 10 i=1,m
      do 10 j=1,n
      ua(i,j)=ub(i,j)
      va(i,j)=vb(i,j)
  10  za(i,j)=zb(i,j)
      return
      end

!     ���̶��߽�ֵ
      subroutine tbv(ua,va,za,ub,vb,zb,m,n)
      dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n)
      m1=m-1
      n1=n-1
      do 10 i=1,m
      do 10 j=1,n,n1
      ua(i,j)=ub(i,j)
      va(i,j)=vb(i,j)
  10  za(i,j)=zb(i,j)
      do 20 i=1,m,m1
      do 20 j=1,n
      ua(i,j)=ub(i,j)
      va(i,j)=vb(i,j)
  20  za(i,j)=zb(i,j)
      return
      end

!     ��ӡ���ֲ���
      subroutine pv(na,nb)
      write(*,100)na,nb
 100  format(//////5x,3hna=,i3,5x,3hnb=,i2/)
      return
      end 
