!����1 ���ĸ�ʱ���ļ���վ�ţ����ȣ�γ�ȣ���ˮ����ȫ��д��rain.txt
!��ʱ���ļ�ֻ���filename90.txt
program micaps_surface_data_read
    parameter(nt=4) !4��ʱ������
    integer::n,num,n1,n2,n3,n4
    real,allocatable::lon(:),lat(:),rd(:)    !��̬����ֱ���������վ��Ľ�ˮ
    real:: temp(9)
    character*8,allocatable::sta(:)      !��̬���飬����վ�ţ�վ�ű�����8���ַ���5���Ļ���stnmap��������   
    character*44 filename(nt)        !���ڴ���4��ʱ�ε��ļ���
    open(1,file='D:\ncl_related\data\micaps_rain\filename90.txt')
        do i=1,nt
               read(1,*) filename(i)
               print*,'Filename:',filename(i)
        end do
    close(1) 
    !pause             
    open(2,file='D:\ncl_related\data\micaps_rain\rain.txt')  !����4��ʱ�εĽ�ˮ
      do k=1,nt
      !����1��ʱ��,����ʱ�ε�վ������ֵ��n
        open(3,file=filename(k))
        read(3,*)
        read(3,*) n1,n2,n3,n4,n    
      !���䶯̬����
        allocate(lat(n))  
        allocate(lon(n))
        allocate(rd(n))
        allocate(sta(n))
        close(3)
        !��֪n�����¶�ȡ����
        open(4,file=filename(k))   
        read(4,*)
        read(4,*)
          do i=1,n        !ֻ����ˮ�����Ҫ������
             read(4,*) sta(i),lon(i),lat(i),(temp(j),j=1,9),rd(i)
             write(2,*)  sta(i),lat(i),lon(i),rd(i)  !һ��ʱ�εĽ�ˮ��������� 
          end do
       !�ͷŶ�̬���� 
        deallocate(lat)
        deallocate(lon)
        deallocate(rd)
        deallocate(sta)
        !�ص���ʱ���ļ�
        close(4)
      enddo
    close(2)
      end program micaps_surface_data_read