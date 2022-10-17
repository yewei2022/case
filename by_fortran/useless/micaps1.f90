!步骤1 读四个时次文件的站号，经度，纬度，降水，并全部写入rain.txt
!换时次文件只需改filename90.txt
program micaps_surface_data_read
    parameter(nt=4) !4个时次数据
    integer::n,num,n1,n2,n3,n4
    real,allocatable::lon(:),lat(:),rd(:)    !动态数组分别用来储存站点的降水
    real:: temp(9)
    character*8,allocatable::sta(:)      !动态数组，储存站号，站号必须是8个字符，5个的话会stnmap出不出来   
    character*44 filename(nt)        !用于储存4个时次的文件名
    open(1,file='D:\ncl_related\data\micaps_rain\filename90.txt')
        do i=1,nt
               read(1,*) filename(i)
               print*,'Filename:',filename(i)
        end do
    close(1) 
    !pause             
    open(2,file='D:\ncl_related\data\micaps_rain\rain.txt')  !储存4个时次的降水
      do k=1,nt
      !读第1个时次,将该时次的站点数赋值于n
        open(3,file=filename(k))
        read(3,*)
        read(3,*) n1,n2,n3,n4,n    
      !分配动态数组
        allocate(lat(n))  
        allocate(lon(n))
        allocate(rd(n))
        allocate(sta(n))
        close(3)
        !获知n后，重新读取数据
        open(4,file=filename(k))   
        read(4,*)
        read(4,*)
          do i=1,n        !只读降水，别的要素暂略
             read(4,*) sta(i),lon(i),lat(i),(temp(j),j=1,9),rd(i)
             write(2,*)  sta(i),lat(i),lon(i),rd(i)  !一个时次的降水场输入完毕 
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
      end program micaps_surface_data_read