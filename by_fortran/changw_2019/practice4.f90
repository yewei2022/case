program practice4
implicit none
external mean,covar
integer i,m0
real Fa,F
parameter (m0=20,Fa=4.41)
integer year(m0)
real xi(m0),yi(m0),a,b,r
open (1,file="D:\YEWEI\practice4\pracitce4_data.dat")
do i=1,20
read (1,*) year(i),yi(i),xi(i)
end do
close(1)
call covar(m0,xi,yi,a,b,r)!调用子程序计算系数
F=r**2/((1-r**2)/(m0-2))
print *,"y=",a,"+",b,"*x"
if (F>Fa)then
print *,'F=',F,'>','Fa=4.41','拒绝原假设，回归方程显著'
else
print *,'F=',F,'<','Fa=4.41','接受原假设，回归方程不显著'
end if
pause
end program
subroutine covar(m0,x,y,a,b,r) !子例行程序,计算回归系数、相关系数
implicit none
integer i,m0
real x(m0),y(m0),xd(m0),yd(m0),xsta(m0),ysta(m0)
real a,b
real xmean,ymean,xave,yave,r
call mean(m0,x,xave,xd,xsta,xmean)
call mean(m0,y,yave,yd,ysta,ymean)
b=0
do i=1,m0
   b=b+xd(i)*yd(i)
end do
b=b/(m0*xmean**2)
a=yave-b*xave
r=0
do i=1,m0
r=xsta(i)*ysta(i)+r
end do
r=r/m0
end subroutine covar
subroutine mean(m0,x,ave,xd,xsta,xmean)!子例行程序，计算均值，距平，均方差,并化为标准化变量
implicit none
integer i,m0
real x(m0),xd(m0),xsta(m0),xmean,ave,sum
sum=0
 do i=1,m0
     sum=sum+x(i)
 end do
 ave=sum/m0
 do i=1,m0
   xd(i)=x(i)-ave
   xmean=xd(i)**2+xmean
 end do
 xmean=sqrt(xmean/m0)
 do i=1,m0
    xsta(i)=xd(i)/xmean
 end do
 end subroutine mean

