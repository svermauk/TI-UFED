PROGRAM mtd_probability
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddiff1, gridmin2, gridmax2, griddiff2
REAL*8 :: gridmin3, gridmax3, griddiff3, gridmin4, gridmax4, griddiff4
REAL*8 :: kt0, kt, ktb, bias_fact, alpha_cv, deltaT, den
REAL*8 :: diff_s2, ds2, ss, hh, dum, dum1, num, gamma_, s1, s2, s3, s4
REAL*8 :: phi_z, num_phi_z, data1, prob1
REAL*8, ALLOCATABLE :: prob(:,:,:), av_dhdl(:,:,:)
INTEGER, ALLOCATABLE :: norm_dhdl(:,:,:)
REAL*8, ALLOCATABLE :: cv1(:), cv2(:), vbias(:), dhdl(:)
REAL*8, ALLOCATABLE :: cv3(:), cv4(:), ct(:)
INTEGER :: mtd_steps, md_steps, i, t_min, t_max
INTEGER :: mtd_max, w_cv, w_hill, i_mtd, i_md
INTEGER :: i_s1, i_s2, i_s3, i_s4, nbin1, nbin2, nbin3, nbin4
INTEGER :: index1, index2, index3, index4

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

OPEN(1,FILE='input',STATUS='old')
!OPEN(2,FILE='PROB.dat',STATUS='replace',form='unformatted')
OPEN(3,FILE='Pu.dat',STATUS='unknown')
OPEN(11,FILE='COLVAR',STATUS='old')
OPEN(12,FILE='HILLS',STATUS='old')
OPEN(13,FILE='ti001.out',STATUS='unknown')
OPEN(21,FILE='data_ct.dat',STATUS='old')

CALL get_steps(11,md_steps)
CALL get_steps(12,mtd_steps)

print *, 'md_steps=', md_steps, 'mtd_steps=', mtd_steps

read(1,*) kt0, kt, bias_fact
read(1,*) t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'

read(1,*) gridmin1, gridmax1, griddiff1
read(1,*) gridmin2, gridmax2, griddiff2
read(1,*) gridmin3, gridmax3, griddiff3
read(1,*) w_cv, w_hill

! deltaT=(bias_fact - 1.d0)*kt
! deltaT = (bias_fact - 1.d0)*kt0
! alpha_cv = (kt + deltaT)/deltaT
WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max

kt = kb*kt
kt0 = kb*kt0
! gamma_ = 1/alpha
gamma_ = (bias_fact - 1.0)/bias_fact
write(*,*) 'gamma_=', gamma_

ALLOCATE(cv1(md_steps),cv2(md_steps),cv3(md_steps),dhdl(md_steps))
ALLOCATE(vbias(md_steps),ct(mtd_steps))

DO i_md=1,t_max
   READ(11,*) dum, cv1(i_md),cv2(i_md),cv3(i_md),vbias(i_md)
!   READ(11,*) dum, cv1(i_md),cv2(i_md),cv3(i_md)
   READ(13,*)dhdl(i_md)

   IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
   IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
!   IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
!   IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0
END DO

nbin1 = NINT((gridmax1-gridmin1)/griddiff1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddiff2)+1
nbin3 = NINT((gridmax3-gridmin3)/griddiff3)+1
write(*,*) nbin1, nbin2, nbin3

!write(*,*) "reading vbias.dat file"

!do i_md=1,t_max
!  read(22,*) dum, vbias(i_md)
!end do

write(*,*) "reading ct.dat file"
do i_mtd=1,mtd_steps
  read(21,*) dum, ct(i_mtd)
end do

ALLOCATE(prob(nbin1,nbin2,nbin3))
ALLOCATE(av_dhdl(nbin1,nbin2,nbin3))
ALLOCATE(norm_dhdl(nbin1,nbin2,nbin3))

WRITE(*,*) 'calculating  probability'

den=0.d0
prob=0.d0
av_dhdl=0.d0
dum1=0.d0
DO i_md=1,md_steps
 IF((i_md.GT.t_min).AND.(i_md.LE.t_max))THEN
    index1 = nint((cv1(i_md)-gridmin1)/griddiff1) +1
    index2 = nint((cv2(i_md)-gridmin2)/griddiff2) +1
    index3 = nint((cv3(i_md)-gridmin3)/griddiff3) +1
    if(index1.eq.nbin1)index1=1 ! WARNING! this is to avoid -pi and +pi bins to be counted seperately
    if(index2.eq.nbin2)index2=1 ! WARNING! this is to avoid -pi and +pi bins to be counted seperately
    if(index3.eq.nbin3)index3=1 ! WARNING! this is to avoid -pi and +pi bins to be counted seperately
    IF(index1.gt.0.and.index2.gt.0.and.index3.gt.0.and.index1.le.nbin1.and.index2.le.nbin2.and.index3.le.nbin3) then
       if((gridmin1.le.cv1(i_md)).and.(gridmin2.le.cv2(i_md)).and.(gridmin3.le.cv3(i_md)).and.   &
          (gridmax1.ge.cv1(i_md)).and.(gridmax2.ge.cv2(i_md)).and.(gridmax3.ge.cv3(i_md))) then
!         i_mtd=(i_md*w_cv/w_hill)+1
          i_mtd=((i_md-1)*w_cv/w_hill)!+1
          dum=(vbias(i_md)*kj_to_kcal) - ct(i_mtd)                      ! Use this when vbias is calculated from PLUMED
!          dum=vbias(i_md) - ct(i_md)                                    ! Use this when vbias and ct is calculated from PLUMED
!          dum=dum*kj_to_kcal                                            ! Use this when vbias and ct is calculated from PLUMED
          prob(index1,index2,index3) = prob(index1,index2,index3) +  DEXP(dum/kt0)
          av_dhdl(index1,index2,index3)=av_dhdl(index1,index2,index3)+dhdl(i_md)
          norm_dhdl(index1,index2,index3)=norm_dhdl(index1,index2,index3)+1
!          print*,dum,kt0,dum/kt0,DEXP(dum/kt0) 
       end if
    END IF
 END IF
END DO

DO index1=1,nbin1
   DO index2=1,nbin2
     DO index3=1,nbin3
       den=den+prob(index1,index2,index3)
     end do
   end do
end do

!  den=den*griddiff1*griddiff2


DO i_s1=1,nbin1
   s1=DFLOAT(i_s1-1)*griddiff1+gridmin1
   DO i_s2=1,nbin2
      s2=DFLOAT(i_s2-1)*griddiff2+gridmin2
      DO i_s3=1,nbin3
         s3=DFLOAT(i_s3-1)*griddiff3+gridmin3
!        prob(i_s1,i_s2)=prob(i_s1,i_s2)/den
         prob1=prob(i_s1,i_s2,i_s3)
         dum1=dfloat(norm_dhdl(i_s1,i_s2,i_s3)) !dum1 has the number of visits on a given bin
         data1=0.d0
         if (prob1 > 0.0) data1 = av_dhdl(i_s1,i_s2,i_s3)/dum1 ! Local average <du/dl> in bin
         if(prob1+1.eq.prob1) prob1=1.0d-16  !remove infinity
         phi_z = -(kt)*dlog(prob1/(den*griddiff1*griddiff2*griddiff3))
         if (phi_z > 0.0) num_phi_z=dexp(-kt0*phi_z)                     ! Numerator of A_\lambda(Z) in our paper
          
!        WRITE(2)prob(i_s1,i_s2)
         WRITE(3,'(5E16.8)')s1,s2,s3,num_phi_z,data1
      END DO
      WRITE(3,*)
   END DO
   WRITE(3,*)
END DO
      WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

close(1)
close(2)
close(3)
close(11)
close(12)
close(13)
close(21)

DEALLOCATE(cv1, cv2, cv3, vbias, ct)
DEALLOCATE(prob)
DEALLOCATE(av_dhdl)

END PROGRAM

SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER :: iunit, nsteps
INTEGER :: ios
nsteps=0
REWIND(iunit)
  Read_Loop: DO
  READ(iunit,*,IOSTAT=ios)
    IF(ios.ne.0)EXIT Read_Loop
    nsteps=nsteps+1
  END DO Read_Loop
REWIND(iunit)
END SUBROUTINE
