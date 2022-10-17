!步骤1 读12个时次文件的站点数，站号，经度，纬度，降水，并全部写入rain.txt
!将各个时次站点数写入sta_num.txt，每次重新运行都需删掉之前生成的文件
program rain_read
  parameter(nt=12) !12个时次数据
  integer n,stat !n文本行数，-1=站点数
  real,allocatable::lon(:),lat(:),rd(:)    !动态数组分别用来储存站点的经纬度和降水量
  integer::sta_num(nt) !
  real:: temp(5) !读数据时存取不需要的数据
  character*8,allocatable::sta(:)      !动态数组，储存站号，站号必须是8个字符，5个的话会stnmap出不出来   
  character*70 filename(nt)        !用于储存各时次的文件名
  open(1,file='D:\ncl_related\data\precipitation\6h\filename_6h.txt')
    do i=1,nt
      read(1,*) filename(i)
      print*,'Filename:',filename(i)
    end do
  close(1) 
  !pause           
  open(2,file='D:\ncl_related\data\precipitation\6h\rain.txt')  !储存12个时次的降水
  open(5,file='D:\ncl_related\data\precipitation\6h\sta_num.txt')  !储存12个时次的站点数
  !读取各个时次站点数
  do k=1,nt
    !读第1个时次,将该时次的站点数赋值于n
    n=0
    open(3,file=filename(k))
    do while (.TRUE.)
      read(3,*,iostat=stat)
      if (stat .NE.0) exit
        n=n+1      
    end do
    !将站点数写入文件
    sta_num(k)=n-1 
    !print*,sta_num(k)
    !pause
    write(5,'(I5)') sta_num(k)  
    close(3)
  end do
  !pause
  do k=1,nt
    n=0 
    n=sta_num(k) 
    !根据站点数n分配动态数组
    allocate(lat(n))  
    allocate(lon(n))
    allocate(rd(n))
    allocate(sta(n)) 
    !获知n后，重新读取数据
    open(4,file=filename(k))   
    read(4,*)  !读第一行
    do i=1,n        !只读降水，别的要素暂略
        read(4,*) (temp(j),j=1,4),sta(i),lat(i),lon(i),temp(5),rd(i)
        write(2,*)  sta(i),lat(i),lon(i),rd(i)  
        !一个时次的降水场输入完毕 ，'(I5,1x,f5.2,1x,f6.2,1x,f3.1)',不知道为啥这样写不进去
    end do
    !释放动态数组 
    deallocate(lat)
    deallocate(lon)
    deallocate(rd)
    deallocate(sta)
    !关掉该时次文件
    close(4)
  enddo
close(2)
close(5)
end program rain_read