!����9��ʱ���ļ���55690�����ǣ���56227�����ܣ���56434�����磩�ĵ����¶�
program micaps_surface_tmp_read
    parameter(nt=9) !9��ʱ������
    integer::n,num,n1,n2,n3,n4
    integer,allocatable::T(:)    !��̬����ֱ���������վ��Ľ�ˮ
    real:: useless(18)
    character*5,allocatable::sta(:)      !��̬���飬����վ��
    character*36 filename(nt)        !���ڴ���9��ʱ�ε��ļ���
    open(1,file='D:\ncl_related\data\tmp\filename.txt')
        do i=1,nt
               read(1,*) filename(i)
               print*,'Filename:',filename(i)
        end do
    close(1) 
    !pause             
    open(2,file='D:\ncl_related\data\tmp\surtmp56434.txt')  !����9��ʱ�εĵ����¶�
      do k=1,nt
        !����1��ʱ��,����ʱ�ε�վ������ֵ��n
        open(3,file=filename(k))
        read(3,*)
        read(3,*) n1,n2,n3,n4,n    !n����վ����
        !���䶯̬����
        allocate(T(n))
        allocate(sta(n))
        close(3)
        !��֪n�����¶�ȡ����
        open(4,file=filename(k))   
        read(4,*)
        read(4,*)
          do i=1,n        !ֻ���¶�(��20��)�����Ҫ������,
             read(4,*) sta(i),(useless(j),j=1,18),T(i)
             if (sta(i).eq."56434") then
             write(2,*)  sta(i),T(i)  !һ��ʱ�εĽ�ˮ���������  
             endif
          end do
        !�ͷŶ�̬���� 
        deallocate(T)
        deallocate(sta)
        !�ص���ʱ���ļ�
        close(4)
      enddo
    close(2)
      end program micaps_surface_tmp_read