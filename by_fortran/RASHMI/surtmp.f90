!读出9个时次文件的55690（错那），56227（波密），56434（察隅）的地面温度
program micaps_surface_tmp_read
    parameter(nt=9) !9个时次数据
    integer::n,num,n1,n2,n3,n4
    integer,allocatable::T(:)    !动态数组分别用来储存站点的降水
    real:: useless(18)
    character*5,allocatable::sta(:)      !动态数组，储存站号
    character*36 filename(nt)        !用于储存9个时次的文件名
    open(1,file='D:\ncl_related\data\tmp\filename.txt')
        do i=1,nt
               read(1,*) filename(i)
               print*,'Filename:',filename(i)
        end do
    close(1) 
    !pause             
    open(2,file='D:\ncl_related\data\tmp\surtmp56434.txt')  !储存9个时次的地面温度
      do k=1,nt
        !读第1个时次,将该时次的站点数赋值于n
        open(3,file=filename(k))
        read(3,*)
        read(3,*) n1,n2,n3,n4,n    !n储存站点数
        !分配动态数组
        allocate(T(n))
        allocate(sta(n))
        close(3)
        !获知n后，重新读取数据
        open(4,file=filename(k))   
        read(4,*)
        read(4,*)
          do i=1,n        !只读温度(第20个)，别的要素暂略,
             read(4,*) sta(i),(useless(j),j=1,18),T(i)
             if (sta(i).eq."56434") then
             write(2,*)  sta(i),T(i)  !一个时次的降水场输入完毕  
             endif
          end do
        !释放动态数组 
        deallocate(T)
        deallocate(sta)
        !关掉该时次文件
        close(4)
      enddo
    close(2)
      end program micaps_surface_tmp_read