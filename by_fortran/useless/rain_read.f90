!����1 ��12��ʱ���ļ���վ������վ�ţ����ȣ�γ�ȣ���ˮ����ȫ��д��rain.txt
!������ʱ��վ����д��sta_num.txt��ÿ���������ж���ɾ��֮ǰ���ɵ��ļ�
program rain_read
  parameter(nt=12) !12��ʱ������
  integer n,stat !n�ı�������-1=վ����
  real,allocatable::lon(:),lat(:),rd(:)    !��̬����ֱ���������վ��ľ�γ�Ⱥͽ�ˮ��
  integer::sta_num(nt) !
  real:: temp(5) !������ʱ��ȡ����Ҫ������
  character*8,allocatable::sta(:)      !��̬���飬����վ�ţ�վ�ű�����8���ַ���5���Ļ���stnmap��������   
  character*70 filename(nt)        !���ڴ����ʱ�ε��ļ���
  open(1,file='D:\ncl_related\data\precipitation\6h\filename_6h.txt')
    do i=1,nt
      read(1,*) filename(i)
      print*,'Filename:',filename(i)
    end do
  close(1) 
  !pause           
  open(2,file='D:\ncl_related\data\precipitation\6h\rain.txt')  !����12��ʱ�εĽ�ˮ
  open(5,file='D:\ncl_related\data\precipitation\6h\sta_num.txt')  !����12��ʱ�ε�վ����
  !��ȡ����ʱ��վ����
  do k=1,nt
    !����1��ʱ��,����ʱ�ε�վ������ֵ��n
    n=0
    open(3,file=filename(k))
    do while (.TRUE.)
      read(3,*,iostat=stat)
      if (stat .NE.0) exit
        n=n+1      
    end do
    !��վ����д���ļ�
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
    !����վ����n���䶯̬����
    allocate(lat(n))  
    allocate(lon(n))
    allocate(rd(n))
    allocate(sta(n)) 
    !��֪n�����¶�ȡ����
    open(4,file=filename(k))   
    read(4,*)  !����һ��
    do i=1,n        !ֻ����ˮ�����Ҫ������
        read(4,*) (temp(j),j=1,4),sta(i),lat(i),lon(i),temp(5),rd(i)
        write(2,*)  sta(i),lat(i),lon(i),rd(i)  
        !һ��ʱ�εĽ�ˮ��������� ��'(I5,1x,f5.2,1x,f6.2,1x,f3.1)',��֪��Ϊɶ����д����ȥ
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
close(5)
end program rain_read