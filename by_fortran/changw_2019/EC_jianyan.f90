!  PROGRAM: EC��T639���ϼ���
program data_jianyan   ! main program
    implicit none
    integer nxy,xy !������������
    real  quece !ȱʧֵ
    parameter(nxy=145*37,quece=999999,xy=17*29)
    real irf(nxy),iro(nxy),irfr(xy),iror(xy) !Ԥ��ֵ�͹۲�ֵ
    integer   i,j,k,aa,dd !����վλ����
    real  bb,cc,ee !������վλ����
    real  acc,RMTC !��ƽ���ϵ���;��������
    open(1,file='E:\��ֵʵϰ\����\ec500\16112408.000',status='old',action='read')
    open(2,file='E:\��ֵʵϰ\����\ec500\16112108.072',status='old',action='read')
    read(1,*)
    read(1,*)
    !��ȡ�۲�ֵ(16112108.000)��iro(nxy)
    do  j=1,37
      do  i=1,145
        read(1,*)aa,bb,cc,dd,iro((j-1)*145+i) 
      enddo
      !read(1,*)
    enddo
    !��ȡԤ��ֵ��irf(nxy)
    read(2,*)
    read(2,*)
    do  j=1,37
      do  i=1,145
        read(2,*)aa,bb,cc,dd,irf((j-1)*145+i)
      enddo
      !read(2,*)
    enddo
    !ɸѡ20-60��N;70-140��E�������Ԥ��ֵirfr�͹۲�ֵiror
    k=1
    do i=9,25
        do j=100,128
            irfr(k)=irf(i*145+j)
            iror(k)=iro(i*145+j)
            k=k+1
        enddo
    enddo
 
    call  SACC(xy,irfr,iror,ACC,quece)
    call  RmTEST(xy,irfr,iror,RMTC,quece)

open(3,file='E:\��ֵʵϰ\ʵϰ3\EC72_acc_rmtc.txt',status='old',access='append')
    write(3,100) ACC,RMTC
    100 format(2(f7.5,1X))

end program data_jianyan


!�����ƽ���ϵ��
!����ֵ����ƽ���ϵ��ACC
SUBROUTINE SACC(MM,irf,iro,ACC,quece)
    implicit none
    integer  MM,LT(MM),nn,I1,I
    real quece,ACC,irf(MM),iro(MM),HAM,FAM,SFH,SF,SH,DFA,DHA
    nn=0
    HAM=0.
    FAM=0.
    SFH=0.
     SF=0.
    SH=0.
    I1=0
    DFA=0.
    DHA=0.
    !���Ԥ��ֵ�͹۲�ֵ����Ч���������ͷŵ�HAM��FAM
    DO I=1,MM
        if(abs(irf(i)).ne.quece.and.abs(iro(i)).ne.quece) then
            I1=I1+1
            LT(I1)=I !��¼ȱ��ֵ���
            nn=nn+1
            HAM=HAM+irf(I)
            FAM=FAM+iro(I)
        endif
    enddo
    !WRITE(*,*) 'nn=',nn, HAM,  FAM
    !�����ƽ���ϵ��acc
    if(nn.eq.0) then    
        ACC=-999.0
    elseif(nn.ne.0) then              
        HAM=HAM/real(nn) !Ԥ��ֵ��ֵ
        FAM=FAM/real(nn) !�۲�ֵ��ֵ
        DO I=1,nn
            I1=LT(I) !ֻ����Чֵ
            DFA=iro(I1)-FAM
            DHA=irf(I1)-HAM
            SFH=SFH+DFA*DHA
            SF=SF+DFA*DFA
            SH=SH+DHA*DHA
        enddo

        if(SF*SH.gt.0) then
            ACC=SFH/(SQRT(SF*SH))
        elseif(SF*SH.le.0) then
            ACC=-999.0
        end if

    endif         
    !write(*,*) 'nn=',nn,'acc=',acc
    
end subroutine SACC

!������������
!����ֵ�����������RMTC
subroutine  RmTEST(MM,irf,iro,RMTC,quece)
    implicit none
    integer MM,NN,i
    real irf(MM),iro(MM),RMTC,quece,sm
    NN=0
    sm=0.
    do  i=1,MM
        if(abs(irf(i)).eq.quece.or.abs(iro(i)).eq.quece) cycle
        NN=NN+1
        sm=sm+(irf(i)-iro(i))**2    
    enddo
    if(NN.eq.0) then
        RMTC=-999.0
    else
        RMTC=SQRT(sm/real(NN))
    endif

end subroutine RmTEST

