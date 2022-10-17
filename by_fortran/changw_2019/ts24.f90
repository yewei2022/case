program main
parameter(nsta=360,nz=7)      ! 用于评分的站点数目
real :: ts(nz),obse(nsta),forc(nsta)      ! 七种降水情况ts得分,观测和预报的降水序列 
integer :: station1(nsta),station2(nsta)
      
!读取数据文件
open (1,file='E:\数值实习\实习4\Rain_t639_obs_common.txt')
    do i=1,nsta
        read(1,*)nn,station1(i),station2(i),obse(i),forc(i)
    enddo
close(1)      

   print *,obse(nsta),forc(nsta)     !检验读取数据是否正确
      
!计算 ts 评分

 yuzhi=0.1
 call qets2(nsta,obse,forc,yuzhi,ts(1))          !无雨
 pause
 
 yuzhi1=0.1
 yuzhi2=10.0
 call qets1(nsta,obse,forc,yuzhi1,yuzhi2,ts(2))  !小雨
 pause

 yuzhi1=10.1
 yuzhi2=25.0
 call qets1(nsta,obse,forc,yuzhi1,yuzhi2,ts(3))  !中雨
  pause
 yuzhi1=25.1
 yuzhi2=50.0
 call qets1(nsta,obse,forc,yuzhi1,yuzhi2,ts(4))  !大雨 
  pause
 yuzhi1=50.1
 yuzhi2=100.0
 call qets1(nsta,obse,forc,yuzhi1,yuzhi2,ts(5))  !暴雨
  pause
 yuzhi1=100.1
 yuzhi2=200.0
 call qets1(nsta,obse,forc,yuzhi1,yuzhi2,ts(6)) !大暴雨
  pause
 yuzhi=200.0
 call qets3(nsta,obse,forc,yuzhi,ts(7))       !特大暴雨
  pause
!输出数据到文件
open(2,file='E:\数值实习\实习4\ts2105_2205_24.txt')
      do i=1,7
        write(2,"(f9.7)") ts(i)     
      enddo
close(2)
end program main
    
subroutine qets1(n,obse,forc,yuzhi1,yuzhi2,s) 
integer :: n                ! 用于评分的站点数目
real :: forc(n),obse(n)     ! 预报和观测的降水序列
real :: s                   ! ts得分
real :: yuzhi1,yuzhi2       ! 小、中、大、暴、大暴雨降水情况下：yuzhi1为最小阈值，yuzhi2为最大阈值      
integer :: A,B,C            ! 命中，空报，漏报
integer :: i
      
     A=0
     B=0
     C=0
     do i=1,n
         
     if(obse(i).ge.yuzhi1.and.obse(i).le.yuzhi2) then
        if(forc(i).ge.yuzhi1.and.forc(i).le.yuzhi2) A=A+1
        if(forc(i).lt.yuzhi1.or.forc(i).gt.yuzhi2) B=B+1  
     else
     if(forc(i).ge.yuzhi1.and.forc(i).le.yuzhi2) C=C+1    
     endif
     
     enddo
     
      write(*,*) A,B,C    
      s=A*1.0/(A+B+C)   
      write(*,*) s         
    end subroutine qets1
    
subroutine qets2(n,obse,forc,yuzhi,s) 
integer :: n                ! 用于评分的站点数目
real :: forc(n),obse(n)     ! 预报和观测的降水序列
real :: s                   ! ts得分
real :: yuzhi               ! 微量或者无雨：yuzhi      
integer :: A,B,C            ! 命中，空报，漏报
integer :: i
     A=0
     B=0
     C=0
     do i=1,n
         
     if(obse(i).lt.yuzhi) then
        if(forc(i).lt.yuzhi) A=A+1
        if(forc(i).ge.yuzhi) B=B+1  
     else
     if(forc(i).lt.yuzhi) C=C+1    
     endif
     
     enddo
     
      write(*,*) A,B,C    
      s=A*1.0/(A+B+C)   
      write(*,*) s         
    end subroutine qets2
      
subroutine qets3(n,obse,forc,yuzhi,s) 
integer :: n                ! 用于评分的站点数目
real :: forc(n),obse(n)     ! 预报和观测的降水序列
real :: s                   ! ts得分
real :: yuzhi               ! 特大暴雨：yuzhi      
integer :: A,B,C            ! 命中，空报，漏报
integer :: i
     A=0
     B=0
     C=0
     do i=1,n
         
     if(obse(i).ge.yuzhi) then
        if(forc(i).ge.yuzhi) A=A+1
        if(forc(i).lt.yuzhi) B=B+1  
     else
     if(forc(i).ge.yuzhi) C=C+1    
     endif
     
     enddo
     
      write(*,*) A,B,C    
      s=A*1.0/(A+B+C)   
      write(*,*) s         
    end subroutine qets3

