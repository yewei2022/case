 PROGRAM EOFPW
 implicit none 
 integer it,i,j,ii,iiii,kv,iik,jj,kk,m,n,mnh,ks,kvt,irec,kkkk
 integer nx,ny
 PARAMETER(m=516,n=216,mnh=216,ks=-1,kv=10,kvt=2)
 !kv为10代表只输出10个特征值
 !kvt为2代表只输出两个模态的eof分析值
 real ff
 PARAMETER(ff=-999.0,nx=18,ny=12)
 DIMENSION F(n,m),A(mnh,mnh),S(mnh,mnh),ER(mnh,4),DF(n),V(mnh),AVF(n),evf(n,kvt),tCF(m,kvt),dat(nx,ny),nf(n)
 real F,A,S,ER,DF,V,AVF,evf,tCF,dat
 integer nf
 open(11,file="D:\YEWEI\practice8\sstpx.grd",form='unformatted',access='direct',recl=nx*ny)
 do 132 it=1,m
 read(11,rec=it)((dat(i,j),i=1,nx),j=1,ny)
 do 132 jj=1,ny
    do 132 ii=1,nx
        kkkk=nx*(jj-1)+ii
        f(kkkk,it)=dat(ii,jj)
        132  continue
 close(11)
      CALL Test1(n,m,ff,f,nf) !输入数据
	write(*,*)'ok2'
      CALL TRANSF(N,M,F,nf,AVF,DF,KS)
      	write(*,*)'ok3'
      CALL FORMA(N,M,MNH,F,A)
		write(*,*)'ok4'
      CALL JCB(MNH,A,S,0.00001)
		write(*,*)'ok5'
      CALL ARRANG(KV,MNH,A,ER,S)
		write(*,*)'ok6'
      CALL TCOEFF(KVT,KV,N,M,MNH,S,F,V,evf,tcf,ER)
		write(*,*)'ok7'
      call test3(N,ff,nf,evf,kvt)
		write(*,*)'ok8'
      open(21,file='D:\YEWEI\practice8\evf.grd',form='unformatted',access='direct',recl=nx*ny)
       irec=0
      do 668 kk=1,kvt
      irec=irec+1
668   write(21,rec=irec)((evf(nx*(j-1)+i,kk),i=1,nx),j=1,ny) 
      close(21)
      open(22,file='D:\YEWEI\practice8\tcf.grd',form='unformatted',access='direct',recl=kvt)
      irec=0
      do 345 it=1,m
      irec=irec+1
345   write(22,rec=irec) tcf(it,1)!写入特征向量的时间系数序列
      close(22)
106   format(10f8.4)
      open(23,file='D:\YEWEI\practice8\dats.dat')
      write(23,106)(er(iiii,3),iiii=1,kv) !写入每个特征向量上投影的总和，即方差贡献率
      close(23)
      END


      SUBROUTINE Test1(n1,m,ff,f,nf)
      DIMENSION F(N1,M)
      DIMENSION nF(N1)
      do i=1,n1
      nf(i)=0.0
      enddo
      do i=1,n1
      do k=1,m
      if(f(i,k).eq.ff)then
      f(i,k)=0.0
      nf(i)=nf(i)+1
      endif
      enddo
      enddo
      return
      end



      SUBROUTINE TRANSF(n1,m,f,nf,avf,df,ks)
! THIS SUBROUTINE PROVIDES INITIAL F BY KS and kv.
      DIMENSION F(N1,M),AVF(N1)
      DIMENSION DF(N1)
      DIMENSION nF(N1)
      if(ks.eq.-1)then
      goto 30
      endif
      do i=1,n1
      avf(i)=0.0
      enddo
      if(ks.eq.0)then
      goto 5
      endif
      do i=1,n1
      df(i)=0.0
      enddo
5      continue
      DO 141 I=1,N1
      if(nf(i).ne.0) goto 141
      do 12 j=1,m
  12  AVF(I)=AVF(I)+F(I,J)
      AVF(I)=AVF(I)/M
      DO 14 J=1,M
 14    F(I,J)=F(I,J)-AVF(I)
 141  CONTINUE
      IF(KS.EQ.0) THEN
      RETURN
      ELSE
      DO 241 I=1,N1
      if(nf(i).ne.0) goto 241
      DO 22 J=1,M
  22  DF(I)=DF(I)+F(I,J)*F(I,J)
      DF(I)=SQRT(DF(I)/M)
      DO 24 J=1,M
  24  F(I,J)=F(I,J)/DF(I)
 241  CONTINUE
      ENDIF
  30  CONTINUE
      RETURN
      END


      SUBROUTINE FORMA(N,M,MNH,F,A)
!!THIS SUBROUTINE FORMS A BY F
      DIMENSION F(N,M),A(MNH,MNH)
      IF(M-N) 40,50,50
  40  DO 44 I=1,MNH
      DO 44 J=I,MNH
      A(I,J)=0.0
      DO 42 IS=1,N
  42  A(I,J)=A(I,J)+F(IS,I)*F(IS,J)
      A(J,I)=A(I,J)
  44  CONTINUE
      RETURN
  50  DO 54 I=1,MNH
      DO 54 J=I,MNH
      A(I,J)=0.0
      DO 52 JS=1,M
  52  A(I,J)=A(I,J)+F(I,JS)*F(J,JS)
      A(J,I)=A(I,J)
  54  CONTINUE
      RETURN
      END


      SUBROUTINE JCB(N,A,S,EPS)
!!THIS SUBROUTINE COMPUTS EIGENVALUES AND EIGENVECTORS OF A
      DIMENSION A(N,N),S(N,N)
      DO 30 I=1,N
      DO 30 J=1,I
      IF(I-J) 20,10,20
  10  S(I,J)=1.
      GO TO 30
  20  S(I,J)=0.
      S(J,I)=0.
  30  CONTINUE
      G=0.
      DO 40 I=2,N
      I1=I-1
      DO 40 J=1,I1
  40  G=G+2.*A(I,J)*A(I,J)
      S1=SQRT(G)
      S2=EPS/FLOAT(N)*S1
      S3=S1
      L=0
  50  S3=S3/FLOAT(N)
  60  DO 130 IQ=2,N
      IQ1=IQ-1
      DO 130 IP=1,IQ1
      IF(ABS(A(IP,IQ)).LT.S3) GOTO 130
      L=1
      V1=A(IP,IP)
      V2=A(IP,IQ)
      V3=A(IQ,IQ)
      U=0.5*(V1-V3)
      IF(U.EQ.0.0) G=1.
      IF(ABS(U).GE.1E-10) G=-SIGN(1.,U)*V2/SQRT(V2*V2+U*U)
      ST=G/SQRT(2.*(1.+SQRT(1.-G*G)))
      CT=SQRT(1.-ST*ST)
      DO 110 I=1,N
      G=A(I,IP)*CT-A(I,IQ)*ST
      A(I,IQ)=A(I,IP)*ST+A(I,IQ)*CT
      A(I,IP)=G
      G=S(I,IP)*CT-S(I,IQ)*ST
      S(I,IQ)=S(I,IP)*ST+S(I,IQ)*CT
  110 S(I,IP)=G
      DO 120 I=1,N
      A(IP,I)=A(I,IP)
  120 A(IQ,I)=A(I,IQ)
      G=2.*V2*ST*CT
      A(IP,IP)=V1*CT*CT+V3*ST*ST-G
      A(IQ,IQ)=V1*ST*ST+V3*CT*CT+G
      A(IP,IQ)=(V1-V3)*ST*CT+V2*(CT*CT-ST*ST)
      A(IQ,IP)=A(IP,IQ)
  130 CONTINUE
      IF(L-1) 150,140,150
  140 L=0
      GO TO 60
  150 IF(S3.GT.S2) GOTO 50
      RETURN
      END


      SUBROUTINE ARRANG(KV,MNH,A,ER,S)
      DIMENSION A(MNH,MNH),ER(mnh,4),S(MNH,MNH)
      TR=0.0
      DO 200 I=1,MNH
      TR=TR+A(I,I)
  200 ER(I,1)=A(I,I)
      MNH1=MNH-1
      DO 210 K1=MNH1,1,-1
      DO 210 K2=K1,MNH1
      IF(ER(K2,1).LT.ER(K2+1,1)) THEN
      C=ER(K2+1,1)
      ER(K2+1,1)=ER(K2,1)
      ER(K2,1)=C
      DO 205 I=1,MNH
      C=S(I,K2+1)
      S(I,K2+1)=S(I,K2)
      S(I,K2)=C
  205 CONTINUE
      ENDIF
  210 CONTINUE
      ER(1,2)=ER(1,1)
      DO 220 I=2,KV
      ER(I,2)=ER(I-1,2)+ER(I,1)
  220 CONTINUE
      DO 230 I=1,KV
      ER(I,3)=ER(I,1)/TR
      ER(I,4)=ER(I,2)/TR
  230 CONTINUE
      WRITE(6,250) TR
  250 FORMAT(/5X,'TOTAL SQUARE ERROR=',F20.5)
      RETURN
      END


      SUBROUTINE TCOEFF(KVT,KV,N,M,MNH,S,F,V,evf,tcf,ER)
!    !THIS SUBROUTINE PROVIDES STANDARD EIGENVECTORS (M.GE.N,SAVED IN S;
!M.LT.N,SAVED IN F) AND ITS TIME COEFFICENTS SERIES (M.GE.N,
!SAVED IN F; M.LT.N,SAVED IN S)
      DIMENSION S(MNH,MNH),F(N,M),V(MNH),ER(mnh,4),evf(n,kvt),tcf(m,kvt)
      DO 360 J=1,KVT
      C=0.
      DO 350 I=1,MNH
  350 C=C+S(I,J)*S(I,J)
      C=SQRT(C)
      DO 160 I=1,MNH
      s(I,J)=S(I,J)/C
  160 evf(I,J)=S(I,J)/C
  360 CONTINUE
      IF(N.LE.M) THEN
      DO 390 J=1,M
      DO 370 I=1,N
      V(I)=F(I,J)
      F(I,J)=0.
  370 CONTINUE
      do 371 is=1,kvt
	tcf(j,is)=0.
371	continue
      DO 380 IS=1,KVT
      DO 380 I=1,N
      f(is,j)=F(IS,J)+V(I)*S(I,IS)
  380 tcf(j,is)=tcf(J,is)+V(I)*S(I,IS)
  390 CONTINUE

      ELSE
      DO 410 I=1,N
      DO 400 J=1,M
      V(J)=F(I,J)
      F(I,J)=0.
  400 CONTINUE
      DO 410 JS=1,KVT
      DO 410 J=1,M
      f(I,JS)=F(I,JS)+V(J)*S(J,JS)
  410 CONTINUE
      DO 430 JS=1,KVT
      DO 420 J=1,M
      tcf(J,JS)=S(J,JS)*SQRT(ER(JS,1))
  420 CONTINUE
      DO 430 I=1,N
      evf(I,JS)=F(I,JS)/SQRT(ER(JS,1))
  430 CONTINUE
      t=0.0
      do 3650 i=1,m
3650   t=t+tcf(i,1)*tcf(i,2)
      write(*,*)t
	t=0.0
      do 3651 i=1,n
3651   t=t+evf(i,1)*evf(i,2)
      write(*,*)t
      ENDIF
      RETURN
      END


      SUBROUTINE test3(N1,ff,nf,evf,kvt)
!this subroutine sent undefine value ff to evf
      dimension nf(n1),evf(n1,kvt)
      do i=1,n1
      if(nf(i).ne.0)then
      do k=1,kvt
      evf(i,k)=ff
      enddo
      endif
      enddo
      return
      end

      SUBROUTINE OUTER(KV,ER,mnh)
!THIS SUBROUTINE PRINTS ARRAY ER
!ER(KV,1) FOR  SEQUENCE OF EIGENVALUE FROM BIG TO SMALL
!ER(KV,2) FOR  EIGENVALUE FROM BIG TO SMALL
!ER(KV,3) FOR  SMALL LO=(LAMDA/TOTAL VARIANCE)
!ER(KV,4) FOR  BIG LO=SUM OF SMALL LO)
      DIMENSION ER(mnh,4)
      WRITE(16,510)
  510 FORMAT(/30X,'EIGENVALUE AND ANALYSIS ERROR',/)
      WRITE(16,520)
  520 FORMAT(10X,1HH,8X,5HLAMDA,10X,6HSLAMDA,11X,2HPH,12X,3HSPH)
      WRITE(16,530) (IS,(ER(IS,J),J=1,4),IS=1,KV)
  530 FORMAT(1X,I10,4F15.5)
      WRITE(16,540)
  540 FORMAT(//)
      RETURN
      END


      SUBROUTINE OUTVT(KVT,N,M,MNH,S,F,EGVT,ECOF)
! THIS SUBROUTINE PRINTS STANDARD EIGENVECTORS
!AND ITS TIME-COEFFICENT SERIES
      DIMENSION F(N,M),S(MNH,MNH),EGVT(N,KVT),ECOF(M,KVT)
      WRITE(16,560)
  560 FORMAT(30X,'STANDARD EIGENVECTORS',/)
      WRITE(16,570) (IS,IS=1,KVT)
  570 FORMAT(3X,80I7)
      DO 550 I=1,N
      IF(M.GE.N) THEN
      WRITE(16,580) I,(S(I,JS),JS=1,KVT)
  580 FORMAT(1X,I3,80F7.3,/)
      DO 11 JS=1,KVT
      EGVT(I,JS)=S(I,JS)
   11 CONTINUE
      ELSE
      WRITE(16,590) I,(F(I,JS),JS=1,KVT)
  590 FORMAT(1X,I3,80F7.3)
      DO 12 JS=1,KVT
      EGVT(I,JS)=F(I,JS)
   12 CONTINUE
      ENDIF
  550 CONTINUE
      WRITE(16,720)
  720 FORMAT(//)
      WRITE(16,610)
  610 FORMAT(30X,'TIME-COEFFICENT SERIES OF S. E.'/)
      WRITE(16,620) (IS,IS=1,KVT)
  620 FORMAT(3X,80I7)
      DO 600 J=1,M
      IF(M.GE.N) THEN
      WRITE(16,630) J,(F(IS,J),IS=1,KVT)
  630 FORMAT(1X,I3,80F7.1)
      DO 13 IS=1,KVT
      ECOF(J,IS)=F(IS,J)
   13 CONTINUE
      ELSE
      WRITE(16,640) J,(S(J,IS),IS=1,KVT)
  640 FORMAT(1X,I3,80F7.1)
      DO 14 IS=1,KVT
      ECOF(J,IS)=S(J,IS)
   14 CONTINUE
      ENDIF
  600 CONTINUE
      RETURN
      END



