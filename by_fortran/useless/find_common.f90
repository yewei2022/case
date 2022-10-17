!挑选出12个时次的公共站点，并写入txt文件
program find_common
implicit none
integer nt
parameter(nt=12)
integer::sta_num(nt)
integer ::i,j,k,n1,n2,com
integer,allocatable::sta1(:),sta2(:),sta_com1(:),sta_com(:)
!使用动态数组之前记得先分配
!读取各个时次的站点数量
open(1,file='D:\ncl_related\data\precipitation\6h\sta_num.txt') 
do i=1,nt 
  read(1,*) sta_num(i)
end do
close(1)
!先挑选第1、2个时次的相同站点数
open(2,file='D:\ncl_related\data\precipitation\6h\rain.txt')
n1=sta_num(1) 
n2=sta_num(2) 
!根据站点数n1、n2分配动态数组
allocate(sta1(n1))  
allocate(sta2(n2))
do k=1,n1
  read(2,*) sta1(k) !第一个时次站号,这里读不进来
end do
do k=1,n2
  read(2,*) sta2(k) 
end do
call pick(n1,n2,sta1,sta2,com,sta_com)
deallocate(sta1)
deallocate(sta2)
do i=3,nt
  n1=com !前两个文件的公共站点数
  n2=sta_num(i) !下一个时次站点数
  allocate(sta_com1(n1))  
  allocate(sta2(n2))
  do k=1,n1
    sta_com1(k)=sta_com(k) !前两个文件的公共站号
  end do
  deallocate(sta_com) !释放前两个文件的公共站号数组
  do k=1,n2    
    read(2,*) sta2(k) !下一个时次的站号
  end do
  call pick(n1,n2,sta_com1,sta2,com,sta_com)
  deallocate(sta_com1)
  deallocate(sta2) 
end do
close(2)
open(3,file="D:\ncl_related\data\precipitation\6h\station_com.txt")
do i=1,com
 write(3,'(i4,i8)') i,sta_com(i)     
end do
close(3)
pause
end program find_common


!挑选出两个时次的公共站点数量，公共站号
subroutine pick(n1,n2,sta1,sta2,com,sta_com)
integer ::n1,n2,stat,com
integer ::sta1(:),sta2(:)
integer,allocatable::sta_com(:)
!统计公共站点数
com=0
do i=1,n1
    stat=0
    do j=1,n2
        if (sta1(i).eq.sta2(j)) then 
            stat=1
            exit
        endif
    enddo
    if (stat.eq.1)  then 
        com = com+1  
    end if
enddo
allocate(sta_com(com)) 
!存入公共站号
com=0
do i=1,n1
    stat=0
    do j=1,n2
        if (sta1(i).eq.sta2(j)) then 
            stat=1
            exit
        endif
    enddo
    if (stat.eq.1)  then 
        com = com+1
        sta_com(com)=sta1(i)   
    end if
enddo
end subroutine pick

