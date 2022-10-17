program AWX2GRD
!转码awx资料，用于ncl绘图

    implicit none    
    integer,parameter :: N=1201

    character*1 :: ch(N,N+2)
    integer :: tbb(N,N),i,j

    open(10,file='D:\ncl_related\data\TBB\FY2D_DATA\FY2D_TBB_IR1_OTG_20081028_0630.AWX',&
    form='unformatted',access='direct',recl=1201*1023)
    open(20,file='D:\ncl_related\data\TBB\FY2D_TBB_IR1_OTG_20081028_0630.grd',form='unformatted')
       
    read(10,rec=1) ch

    tbb = ichar(ch(:,3:N+2)) !ichar是取字符串的ASCII码
    tbb = tbb+100

    write(20)((tbb(j,i),j=1,N),i=N,1,-1)
    PAUSE

end program AWX2GRD