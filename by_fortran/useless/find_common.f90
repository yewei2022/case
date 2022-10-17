!��ѡ��12��ʱ�εĹ���վ�㣬��д��txt�ļ�
program find_common
implicit none
integer nt
parameter(nt=12)
integer::sta_num(nt)
integer ::i,j,k,n1,n2,com
integer,allocatable::sta1(:),sta2(:),sta_com1(:),sta_com(:)
!ʹ�ö�̬����֮ǰ�ǵ��ȷ���
!��ȡ����ʱ�ε�վ������
open(1,file='D:\ncl_related\data\precipitation\6h\sta_num.txt') 
do i=1,nt 
  read(1,*) sta_num(i)
end do
close(1)
!����ѡ��1��2��ʱ�ε���ͬվ����
open(2,file='D:\ncl_related\data\precipitation\6h\rain.txt')
n1=sta_num(1) 
n2=sta_num(2) 
!����վ����n1��n2���䶯̬����
allocate(sta1(n1))  
allocate(sta2(n2))
do k=1,n1
  read(2,*) sta1(k) !��һ��ʱ��վ��,�����������
end do
do k=1,n2
  read(2,*) sta2(k) 
end do
call pick(n1,n2,sta1,sta2,com,sta_com)
deallocate(sta1)
deallocate(sta2)
do i=3,nt
  n1=com !ǰ�����ļ��Ĺ���վ����
  n2=sta_num(i) !��һ��ʱ��վ����
  allocate(sta_com1(n1))  
  allocate(sta2(n2))
  do k=1,n1
    sta_com1(k)=sta_com(k) !ǰ�����ļ��Ĺ���վ��
  end do
  deallocate(sta_com) !�ͷ�ǰ�����ļ��Ĺ���վ������
  do k=1,n2    
    read(2,*) sta2(k) !��һ��ʱ�ε�վ��
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


!��ѡ������ʱ�εĹ���վ������������վ��
subroutine pick(n1,n2,sta1,sta2,com,sta_com)
integer ::n1,n2,stat,com
integer ::sta1(:),sta2(:)
integer,allocatable::sta_com(:)
!ͳ�ƹ���վ����
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
!���빫��վ��
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

