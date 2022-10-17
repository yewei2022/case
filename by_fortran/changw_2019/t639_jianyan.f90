!  PROGRAM: EC��T639���ϼ���
program t639_jianyan   ! main program
    implicit none
    integer nxy,xy !������������
    parameter(nxy=181*91,xy=41*71)
    real irf(nxy),iro(nxy),irfr(xy),iror(xy) !Ԥ��ֵ�͹۲�ֵ
    integer   i,j,k
    real  acc,RMTC !��ƽ���ϵ���;��������
    open(1,file='E:\��ֵʵϰ\����\t639-500\16112408.000',status='old',action='read')
    open(2,file='E:\��ֵʵϰ\����\t639-500\16112108.072',status='old',action='read')
    !��ȡ�۲�ֵ��Ĭ�϶����ĸ�ʽ����Ƕ��ѭ���ľ�γ�ȸ�ʽ����
    read(1,*)
    read(1,*)
    read(1,*) iro

    !��ȡԤ��ֵ��irf(nxy)
    read(2,*)
    read(2,*)
    read(2,*) irf

    !ɸѡ20-60��N;70-140��E�������Ԥ��ֵirfr�͹۲�ֵiror
    k=1
    do i=21,61
        do j=250,320
            irfr(k)=irf(i*181+j)
            iror(k)=iro(i*181+j)
            k=k+1
        enddo
    enddo
 
    call  SACC(xy,irfr,iror,ACC)
    call  RmTEST(xy,irfr,iror,RMTC)

open(3,file='E:\��ֵʵϰ\ʵϰ3\t639_72_acc_rmtc.txt',status='old',access='append')
    !write(3,100) ACC,RMTC
    print 100,ACC,RMTC
    100 format(2(f7.5,1X))
    close(1)
    close(2)
    close(3)
    pause
end program t639_jianyan


!�����ƽ���ϵ��
!����ֵ����ƽ���ϵ��ACC
SUBROUTINE SACC(MM,irf,iro,ACC)
    implicit none
    integer  MM,I
    real ACC,irf(MM),iro(MM),HAM,FAM,SFH,SF,SH,DFA,DHA
    HAM=0.
    FAM=0.
    SFH=0.
     SF=0.
    SH=0.
    DFA=0.
    DHA=0.
    DO I=1,MM
            HAM=HAM+irf(I)
            FAM=FAM+iro(I)
    enddo
    !�����ƽ���ϵ��acc             
        HAM=HAM/real(MM) !Ԥ��ֵ��ֵ
        FAM=FAM/real(MM) !�۲�ֵ��ֵ
        DO I=1,MM
            DFA=iro(I)-FAM
            DHA=irf(I)-HAM
            SFH=SFH+DFA*DHA
            SF=SF+DFA*DFA
            SH=SH+DHA*DHA
        enddo
            ACC=SFH/(SQRT(SF*SH))
    
end subroutine SACC

!������������
!����ֵ�����������RMTC
subroutine  RmTEST(MM,irf,iro,RMTC)
    implicit none
    integer MM,i
    real irf(MM),iro(MM),RMTC,sm
    sm=0.
    do  i=1,MM
        sm=sm+(irf(i)-iro(i))**2    
    enddo
        RMTC=SQRT(sm/real(MM))
end subroutine RmTEST

