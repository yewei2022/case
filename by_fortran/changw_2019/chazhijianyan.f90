!****************************************************************************
!
!  PROGRAM: T639������ϲ�ֵ��վ������
!
!****************************************************************************

   module pra_mod
   implicit none
       integer,parameter :: Nlon=90                               ! T639 E0-180  2
       integer,parameter :: Nlat=45                               ! T639 N0-90   2 
       integer,parameter :: Nsta_max=5000                         ! ���5000վ                      
       real,parameter :: yuzhi=0.1                               ! ���޽�ˮ����ֵ  mm
       character(len = *), parameter :: dir_T639='E:\��ֵʵϰ\��ֵ����\t639-05\'      ! input t639 data 
       character(len = *), parameter :: file1_t639_name='16112102.027'       
	   character(len = *), parameter :: file2_t639_name='16112202.027'       
       character(len = *), parameter :: dir_station='E:\��ֵʵϰ\��ֵ����\'             ! input station data
       character(len = *), parameter :: filename_station='STATIONS.DAT'      
       character(len = *), parameter :: dir_obs='E:\��ֵʵϰ\��ֵ����\sur-05-24\'       ! input original_data 
   	   character(len = *), parameter :: file1_obs_name='16112205.000'        ! input observation_data 
	   character(len = *), parameter :: file2_obs_name='16112305.000'      
	   character(len = *), parameter :: dir_results='E:\��ֵʵϰ\ʵϰ4\'      

	end module pra_mod

	module data_mod
	   use pra_mod
	   implicit none
	   real lon_t639(0:Nlon)
	   real lat_t639(0:Nlat)
	   real Rain_t639(0:Nlon,0:Nlat) !�����㽵ˮ����
	   integer Nsta_t639,Nsta_obs,Nsta_common
	   integer sta_number_t639(Nsta_max) !����վ�����
	   real sta_lat_t639(Nsta_max),sta_lon_t639(Nsta_max),Rain_t639_sta(Nsta_max)  !վγ�ȡ�վ���ȡ�վ�㽵ˮ����
       character(len=10) sta_name_t639(Nsta_max) !վ����
	   real Rain_obs(Nsta_max)
	   integer sta_number_obs(Nsta_max)
	   integer sta_number_common(Nsta_max)
	   real Rain_t639_common(Nsta_max),Rain_obs_common(Nsta_max)

    end module data_mod

	 
	program data_chazhi   ! main program
	   use pra_mod
	   use data_mod
       implicit none
       character(len=50)  cha
	   character(len=10)  cha10
	   integer station,stat,templat,templon,temp1,temp2
	   real lon,lat,gh,rain,ts
	   real  rf,lf,uf,df,lonP,latP
	   integer rn,ln,un,dn,i,j

       !��ֵ����lon_t639(0,2,4,6...180)
	   do i=0,Nlon
	      lon_t639(i)=i*2
		  !write(*,*) "lon",i,lon_t639(i)
	   enddo
       !��ֵγ��lat_t639(0,2,4,6...90)
	   do j=0,Nlat
	      lat_t639(j)=j*2
		  !write(*,*) "lat",j,lat_t639(j)
	   enddo
       !��ȡfile2_t639_name��ÿ2��γ��1�����ݣ��Ľ�ˮ��Ԥ��ֵ��Rain_t639
       write(*,*)'step1'
       Rain_t639(:,:)=0.
       open(22,file=dir_T639//file2_t639_name,form="formatted")
       read(22,*) cha
	   !write(*,*) cha
	   read(22,*) cha
	   !write(*,*) cha
	   do while(1)
          read(22,'(i6,f8.3,f8.3,f8.1,2x,f8.1)',iostat=stat)   station,lon,lat,gh,rain
		  if (stat.lt.0)  exit  ! stat.lt.0  �ļ�����
          i=lon/2
		  j=lat/2
		  !write(*,*) i,j,rain
          Rain_t639(i,j)=rain
	   enddo
	   close(22)
	   i=90
	   j=45
       write(*,'(f8.1,f8.1,4x,f8.1)')   i*2.0,j*2.0,Rain_t639(i,j)
       !��ȡfilename_station��sta_number_t639��վ�ţ�sta_lat_t639��վγ�ȣ�sta_lon_t639��վ���ȣ�sta_name_t639��վ����
       Nsta_t639=0 !���վ�����
       open(22,file=dir_station//filename_station,form="formatted")
	   do while(1)
          read(22,*,iostat=stat)   station,templat,templon,gh,temp1,temp2,cha10
		  if (stat.lt.0)  exit  ! stat.lt.0  �ļ�����
		  Nsta_t639=Nsta_t639+1
          !write(*,*) Nsta_t639,station,templat,templon,cha10
          sta_number_t639(Nsta_t639)=station !��վ��
          sta_lat_t639(Nsta_t639)=templat/100. !д��վ��γ��
          sta_lon_t639(Nsta_t639)=templon/100. !д��վ�㾭��
		  sta_name_t639(Nsta_t639)=cha10	 
	   enddo
	   close(22)
	   write(*,'(a20,i6)')   "stations number",Nsta_t639 !���վ����

       ! test integrate_point
       !lonP=115.7
	   !latP=36.2
       !call integrate_point(lon_t639,Nlon,lat_t639,Nlat,lonP,latP, &
	   !                         rf,lf,uf,df,rn,ln,un,dn)

       !��ֵRain_t639����㣩��Rain_t639_sta��վ�㣩��д��Rain_t639_station2.txt
	   open(22,file=dir_results//"Rain_t639_station2.txt",form="formatted")
	   do i=1,Nsta_t639
          lonP=sta_lon_t639(i) !��ֵվ�����ݾ������꣬����ֵ֪������ʼֵ 
	      latP=sta_lat_t639(i)
          call integrate_point(lon_t639,Nlon,lat_t639,Nlat,lonP,latP, &
	                            rf,lf,uf,df,rn,ln,un,dn)
		  Rain_t639_sta(i)= Rain_t639(rn,un)*rf*uf & !˫���Բ�ֵ��ʽ
		              +Rain_t639(rn,dn)*rf*df &
					  +Rain_t639(ln,dn)*lf*df &
                      +Rain_t639(ln,un)*lf*uf 
          if (((rn+ln).eq.0).or.((un+dn).eq.0)) then 
		     write(*,'(a10,i5,2x,2f5.1)')  "bianjie:",i,sta_lon_t639(i),sta_lat_t639(i)   ! �����߽�
             sta_number_t639(i)=0
		  endif
		  write(22,'(i5,2x,2f5.1,2x,4f6.1,x,4f7.2,2x,i6,f7.2)')  i,sta_lon_t639(i),sta_lat_t639(i), &
		      lon_t639(ln),lon_t639(rn),lat_t639(un),lat_t639(dn), &
		      Rain_t639(rn,un),Rain_t639(rn,dn),Rain_t639(ln,un),Rain_t639(ln,dn),sta_number_t639(i),Rain_t639_sta(i)

	   enddo
	   close(22)


	   !��ȡfile2_obs_name�Ľ�ˮ���۲�ֵ��sta_number_obs��վ�ţ�Rain_obs����ˮ����
       write(*,*)'step2'
       Rain_obs(:)=0
	   open(22,file=dir_obs//file2_obs_name,form="formatted")
	   do i=1,14
          read(22,*) cha
	      !write(*,*) cha
       enddo
	   Nsta_obs=0
	   do while(1)
          read(22,'(i6,f8.3,f8.3,f6.0,f6.1)',iostat=stat)   station,lon,lat,gh,rain
		  if (stat.lt.0)  exit  ! stat.lt.0  �ļ�����
		  Nsta_obs=Nsta_obs+1
          sta_number_obs(Nsta_obs)=station
          Rain_obs(Nsta_obs)=rain
		  !write(*,*)  Nsta_obs,sta_number_obs(Nsta_obs), Rain_obs(Nsta_obs)
	   enddo
	   close(22)
       write(*,'(a20,i6)')   "stations_obs number",Nsta_obs

       !����Ԥ��վ��sta_number_t639(Rain_t639_station2.txt)��˳��
       !��Ԥ��͹۲⽵ˮ��д��Rain_t639_obs_common2.txt��Rain_t639_common,Rain_obs_common
       open(22,file=dir_results//"Rain_t639_obs_common2.txt",form="formatted")
       Nsta_common=0
       do i=1,Nsta_t639
          stat=0
		  do j=1,Nsta_obs
             if (sta_number_t639(i).eq.sta_number_obs(j)) then ! ɸѡ��ͬ��վ�������
			    stat=1
				exit
             endif
		  enddo
		  if (stat.eq.1)  then !Ԥ��ֵ�͹۲�ֵ����
			 Nsta_common = Nsta_common+1  !����ͬվ�����
             write(22,'(i5,2i6,2x,2f7.1)') Nsta_common,sta_number_t639(i),sta_number_obs(j),Rain_t639_sta(i),Rain_obs(j)
			 Rain_t639_common(Nsta_common)=Rain_t639_sta(i) !��ͬһ��վ������ݷֱ�ŵ�t639Ԥ��͹۲����ݱ����д���������
			 Rain_obs_common(Nsta_common)=Rain_obs(j)	 
		  endif
	   enddo
	   write(*,'(a20,i6)')   "Nsta_common number",Nsta_common


       call QETS(Nsta_common,Rain_t639_common(1:Nsta_common),Rain_obs_common(1:Nsta_common),ts) 
       write(*,*) "ts(had rain):",ts
       close(22)
       pause
	endprogram data_chazhi

	!�н�ˮ��TS���ֺ���
	!����ֵ��ts���ֽ��
	subroutine QETS(N,forc,obse,ts)  
	  use pra_mod,only:yuzhi
      implicit none 
	  integer :: N                ! �������ֵ�վ����Ŀ
      real :: forc(N),obse(N)     ! Ԥ���͹۲�Ľ�ˮ����
      real :: ts                  ! ts�÷�
	  integer :: A,B,C            ! ���У��ձ���©��
	  integer :: i
	  A=0
	  B=0
	  C=0
	  write(*,*)'yuzhi',yuzhi
      do i=1,N
         if(forc(i)>=yuzhi.and.obse(i)>=yuzhi) A=A+1 !����
         if(forc(i)>=yuzhi.and.obse(i)<yuzhi)  B=B+1 !�ձ�
         if(forc(i)<yuzhi.and.obse(i)>=yuzhi)  C=C+1 !©��
      enddo
      ts=A*1.0/(A+B+C)
      write(*,*) A,B,C
    end subroutine QETS

    !˫���Բ�ֵ׼������
    !����ֵ��rn,ln,un,dn�Ҿ�������lon_t639����������γ����γ��lat_t639��������rf,lf,uf,df˫���Բ�ֵ���ĸ���γ��ϵ��
	subroutine integrate_point(lonMap,Ncols,latMap,Nrows,lonP,latP, &
	                            rf,lf,uf,df,rn,ln,un,dn)
	   implicit none
	   integer,intent(in) :: Ncols,Nrows
	   real,intent(in) :: lonMap(0:Ncols),latMap(0:Nrows)
	   real,intent(in) :: lonP,latP
	   real,intent(out) :: rf,lf,uf,df
	   integer,intent(out) :: rn,ln,un,dn
	   integer i,j

	   if ((lonP.le.lonMap(0)).or.(lonP.ge.lonMap(Ncols))) then   ! �����߽� rn=ln=0
          rn=0
		  ln=0
		  rf=0.
		  lf=0.
		  write(*,*) " (lonP.le.lonMap(0)).or.(lonP.ge.lonMap(Ncols))",lonP
		  return
	   endif
	   if ((latP.le.latMap(0)).or.(latP.ge.latMap(Nrows))) then   ! �����߽�  un=dn=0
          un=0
		  dn=0
		  uf=0.
		  df=0.
		  write(*,*) "(latP.le.latMap(0)).or.(latP.ge.latMap(Nrows))",latP
		  return
	   endif

       rn=1
       ln=1
       rf=1.
       lf=0.
       do i=1,Ncols
          if (lonP.le.lonMap(i)) exit
       enddo
       rn=i
       ln=i-1
       rf=(lonP-lonMap(i-1))/(lonMap(i)-lonMap(i-1))
       lf=(lonMap(i)-lonP)/(lonMap(i)-lonMap(i-1))             
	   !write(*,'(2i5,6f8.3)') ln,rn,lf,rf,lonP,lonMap(ln),lonMap(rn),lf+rf
       un=1
       dn=1
       uf=1.
       df=0.
       do j=1,Nrows
          if (latP.le.latMap(j)) exit
       enddo
       un=j
       dn=j-1
       uf=(latP-latMap(j-1))/(latMap(j)-latMap(j-1))
       df=(latMap(j)-latP)/(latMap(j)-latMap(j-1))  
       !write(*,'(2i5,6f8.3)') dn,un,df,uf,latP,latMap(dn),latMap(un),df+uf	        
    end subroutine integrate_point






    
