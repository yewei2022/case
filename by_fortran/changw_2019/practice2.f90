 program practice2
 implicit none
 external stand,simplecor,autocorr,maxum
 integer i,j,m0,n0,amm,wmm
 parameter(m0=20)
 parameter(n0=10)
 integer k(n0)
 real amean,wmean,simpcor,amax,wmax
 real anu(m0),win(m0),axd(m0),wxd(m0),asta(m0),wsta(m0)
 real aautocor(n0),wautocor(n0)
 open(1,file='D:\YEWEI\practice2\annual_data.dat')
 read(1,"(20f5.2)") anu
 close(1)
 open(2,file='D:\YEWEI\practice2\winter_data.dat')
 read(2,"(20f5.2)") win
 close(2)
 call stand(m0,anu,amean,axd,asta) !把年平均气温化为标准化变量
 call stand(m0,win,wmean,wxd,wsta) !把冬季平均气温化为标准化变量
 call simplecor (m0,asta,wsta,simpcor)!计算简单相关系数
 print *,"年平均气温和冬季平均气温的相关系数为",simpcor 
 call autocorr(amean,m0,n0,axd,aautocor)!计算年平均气温的自相关系数
 call autocorr(wmean,m0,n0,wxd,wautocor)!计算冬季平均气温的自相关系数
 call maxum(aautocor,n0,amm,amax)!找相关系数最大
 call maxum(wautocor,n0,wmm,wmax)!找相关系数最大
 do i=1,n0
 k(i)=i
 end do
 print *,"年平均气温的自相关系数:"
 write(*,100),k
 write(*,200),aautocor
 print *,"年平均气温自相关系数最大为",amax,"此时的滞后长度为",amm
 print *,"冬季平均气温的自相关系数:"
 write(*,100),k
 write(*,200),wautocor
 print *,"冬季平均气温自相关系数最大为",wmax,"此时的滞后长度为",wmm
 100 format (10i5)
 200 format (10f4.2)
 pause
end program
subroutine stand(m0,x,mean,xd,sta) !子例行程序，计算均方差，距平，并化为标准化变量
implicit none
integer i,m0
real x(m0),xd(m0),sta(m0),ave,mean,sum
sum=0
do i=1,m0
sum=sum+x(i)
end do
ave=sum/m0
mean=0
do i=1,m0
xd(i)=x(i)-ave
mean=xd(i)**2+mean
end do
mean=sqrt(mean/m0)
do i=1,m0
sta(i)=xd(i)/mean
end do
end subroutine stand
subroutine simplecor (m0,x,y,simpcor)!子例行程序，计算相关系数
implicit none
integer i,m0
real x(m0),y(m0),simpcor
simpcor=0
do i=1,m0
simpcor=x(i)*y(i)+simpcor
end do
simpcor=simpcor/m0
end subroutine simplecor
subroutine autocorr(mean,m0,n0,x,autocor)!子例行程序，计算自相关系数
implicit none 
integer i,j,n,m0,n0
real mean,x(m0),y(m0),auto,autocor(n0)
do j=1,n0
n=m0-j
auto=0
  do i=1,n
  y(i)=x(i+j)
  auto=auto+x(i)*y(i)
  end do
auto=auto/(n*mean**2)
autocor(j)=auto
end do
end subroutine autocorr
subroutine maxum(x,n0,mm,max)!子例行程序，找相关系数最大
implicit none
integer i,n0,mm
real x(n0),max
max=0.0
do i=1,n0
  if(abs(x(i))>=max) then
  max=abs(x(i))
  mm=i
  end if
end do
end subroutine maxum