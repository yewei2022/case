program moving_ave_axd !计算时间序列的滑动平均
integer i,j,k,n
parameter(k=11,n=85)
integer yr(n)
real x(n),x1(n-k+1),ny(n),axd(n)
open(2,file="D:\YEWEI\practice7\MA.DAT",form='formatted')
read(2,*)(x(i),i=1,n)
close(2)
do i=1,n
   yr(i)=0
   yr(i)=yr(i)+i+1922
end do
do j=1,n-k+1
  x1(j)=0
  do i=1,k
     x1(j)=x1(j)+x(i+j-1)
  end do
  x1(j)=x1(j)/k
end do
do i=1,n-k+1
   ny(i+5)=x1(i)
end do
call accuxd(n,x,axd)!算累积距平
open(3,file='D:\YEWEI\practice7\ma.grd',form='binary')
do i=1,n
   write(3) x(i)
   write(3) ny(i)
   write(3) axd(i)
end do
close(3)
open(4,file='D:\YEWEI\practice7\ma.txt')
write(4,*)  yr
write(4,*)  x
write(4,*)  ny
write(4,*)  axd
close(4)
end 
subroutine accuxd(m0,x,axd)!子例行程序，计算累积距平
implicit none
integer i,j,m0
real x(m0),xd(m0),axd(m0),ave,sum
sum=0
 do i=1,m0
     sum=sum+x(i)
 end do
 ave=sum/m0
 do i=1,m0
   xd(i)=x(i)-ave
 end do
 do j=1,m0
   axd(j)=0
   do i=1,j
      axd(j)=axd(j)+xd(i)
   end do
 end do
 end subroutine accuxd
