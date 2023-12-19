PROGRAM temp
IMPLICIT NONE
INTEGER :: Nz, counter, n,Tbase_air, Tamp_air,stat
INTEGER :: i,j, info, nstep, Pswr_n, N_winter, T_e
REAL :: dt, dz, t1, t2,alpha,k,z_R,Delta_R, Delta_z,z_e
REAL :: rho, c, Deltat_SS, sigma, eps_s, a,L,somme,hc
REAL:: Tol,power,alpha1,c1,rho1,k1,rk1,rk2
REAL::c2,rho2,k2,alpha2,L1,L2,rdz2,rkdz,rkdrdz,rk1drdz,rk2drdz
REAL:: Ly,dy,rdy2,minalpha,mindiscr,ribbony,ribbonz
REAL:: rribbonyribbonz
INTEGER:: indexl1,indexl2,Ny,m,nribbon, indexribbony1, indexribbony2
INTEGER:: indice_R,Middleribbon,nlinescase,indexribbony1z,indexribbony2z
REAL(KIND(1.0E0)),DIMENSION(:,:),ALLOCATABLE:: w_R,Told,Tnew,Tplot
REAL,DIMENSION(3):: listalpha
REAL,DIMENSION(2):: listdiscr
REAL,DIMENSION(:),ALLOCATABLE:: f_s,T_air,Maxsurf,Tdepth, fD
INTEGER :: countdays,nstepbegin,ndays,indexfak1,indexfak2,indexfak3
CHARACTER(LEN=20)::Filename,Filenamebis
REAL:: avgrho1,avgrho2,avgc1,avgc2,avgk1,avgk2,avgalpha1,avgalpha2
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: c_eff,k_eff,theta,d_T_theta
REAL::latentheat,theta_0_1,theta_0_2,avgd_T_theta1,avgd_T_theta2
REAL(KIND(1.0D0))::c1unfrozen, c1frozen, k1unfrozen, k1frozen,beta
REAL(KIND(1.0D0))::c2unfrozen, c2frozen, k2unfrozen, k2frozen
REAL::exponentlatent,temp_ref, numberindexribbon
INTEGER:: indexNyleft,indexNyright
REAL::depthsnowinit, absorb_snow, c_snow, k_snow, rho_snow
INTEGER::indexsnowinit, indexstartsnow, allocationdepth,parite
INTEGER,DIMENSION(:),ALLOCATABLE:: depthsnow, snowhelp
REAL::alpha_snow,emit_snow,theta_snow
INTEGER,DIMENSION(:),ALLOCATABLE:: indexsnowtot,indexmeltup,indexmeltdown
INTEGER,DIMENSION(:),ALLOCATABLE:: indexsnowdead,indexmeltuptot,indexmeltdowntot
REAL::c_snow_u,c_snow_f,k_snow_u,k_snow_f,c_avg_snow,k_avg_snow,avgd_T_theta_snow
REAL::c_air,k_air,rho_air,temp_melted,T_l,para_gauss
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: rho_eff
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: Taux_value, Taux_time
INTEGER,DIMENSION(:),ALLOCATABLE:: indexmeltedold
INTEGER::indexmeltedold_right, indexsnowtot_right, shift_right,shift_right_cumul
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: Taux_value_right, Taux_time_right
INTEGER::cumul_lateral,m_loop_1,m_loop_2,m_loop_3
INTEGER,DIMENSION(:),ALLOCATABLE::cumul_melted,cumul_vertical
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE::Taux_new,Taux_old
REAL(KIND(1.0D0))::avgrho0
REAL,DIMENSION(:),ALLOCATABLE::heatflux_surf,V_air
REAL::power_acc


!READ PARA SUBROUTINE

NAMELIST /para/dt,Nz,nstep,rho,c,k,Deltat_SS,Pswr_n,Tbase_air,Tamp_air, &
               L,z_R,Delta_R,Delta_z,N_winter,sigma,eps_s,a,T_e,Tol,&
               power,c1,rho1,k1,c2,rho2,k2,z_R,L1,L2,z_e,Ly,Ny,nribbon,ribbony,&
               ribbonz,latentheat,c2unfrozen,c2frozen,k2unfrozen,k2frozen,&
               c1unfrozen,c1frozen,k1unfrozen,k1frozen,theta_0_1,theta_0_2,&
               beta,temp_ref,depthsnowinit,absorb_snow,c_snow,k_snow,rho_snow,&
               emit_snow,theta_snow,c_snow_u,c_snow_f,k_snow_u,k_snow_f,&
               c_air,k_air,rho_air,temp_melted,T_l
               

OPEN(20,FILE='para.dat')

READ(20,NML=para)
!WRITE(*,NML=para)

CLOSE(20)

!END READ PARA SUBROUTINE

!INIT SUBROUTINE

nlinescase = 0

OPEN(36,FILE='Windspeed.dat',STATUS='old')

DO
   READ(36,*,END=22)
   nlinescase = nlinescase + 1
END DO
22 CLOSE(36)

ALLOCATE(V_air(nlinescase))

OPEN(36,FILE='Windspeed.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) V_air(i)
END DO

CLOSE(36)


PRINT*,'...'

nlinescase = 0

OPEN(36,FILE='Airtemp.dat',STATUS='old')

DO
   READ(36,*,END=23)
   nlinescase = nlinescase + 1
END DO
23 CLOSE(36)

ALLOCATE(T_air(nlinescase))

OPEN(36,FILE='Airtemp.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) T_air(i)
END DO

CLOSE(36)


PRINT*,'...'

nlinescase = 0

OPEN(36,FILE='Longwave.dat',STATUS='old')

DO
   READ(36,*,END=24)
   nlinescase = nlinescase + 1
END DO
24 CLOSE(36)

ALLOCATE(fD(nlinescase))

OPEN(36,FILE='Longwave.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) fD(i)
END DO

CLOSE(36)


PRINT*,'...'


alpha = k/(rho*c)
alpha1 = k1/(rho1*c1)
alpha2 = k2/(rho2*c2)

avgrho1 = 2.0*rho*rho1/(rho+rho1)
avgrho2 = 2.0*rho1*rho2/(rho1+rho2)

dz = L/(Nz-1.0)
dy = Ly/(Ny-1.0)

indice_R = INT(z_R/dz)

allocationdepth = Nz

ALLOCATE(Told(allocationdepth,Ny))
ALLOCATE(Tnew(allocationdepth,Ny))
ALLOCATE(c_eff(allocationdepth,Ny))
ALLOCATE(k_eff(allocationdepth,Ny))
ALLOCATE(rho_eff(allocationdepth,Ny))
ALLOCATE(d_T_theta(allocationdepth,Ny))
ALLOCATE(theta(allocationdepth,Ny))
ALLOCATE(Tplot(allocationdepth,Ny))
ALLOCATE(Taux_value(allocationdepth,Ny))
ALLOCATE(Taux_time(allocationdepth,Ny))
ALLOCATE(Taux_value_right(allocationdepth,Ny))
ALLOCATE(Taux_time_right(allocationdepth,Ny))
ALLOCATE(Taux_new(allocationdepth,Ny))
ALLOCATE(Taux_old(allocationdepth,Ny))
ALLOCATE(heatflux_surf(Ny))

DO i = 1,allocationdepth
DO j = 1,Ny

 Told(i,j) = 4.3*(i-1.0)*dz + 1.7

END DO
END DO

Tnew = Told

ALLOCATE(w_R(allocationdepth,Ny))

DO i = 1,allocationdepth
DO j = 1,Ny

w_R(i,j) = 0.0

END DO
END DO

indice_R = (INT(z_R/dz)) + 1
indexribbony1z = indice_R - INT(ribbonz*0.5/dz) + 1
indexribbony2z = indice_R + INT(ribbonz/dz)


indexribbony1 = INT(0.0/dy) + 1
indexribbony2 = INT(ribbony*0.5/dy) + 1

DO i = indexribbony1,indexribbony2
   DO m = indexribbony1z,indexribbony2z
     w_R(m,i) = 1.0
   END DO
END DO


indexl1 = INT(L1/dz) + 1
indexl2 = INT((L1+L2)/dz) + 1

DO i = 1,allocationdepth
DO j = 1,Ny
   numberindexribbon = numberindexribbon + w_R(i,j)
END DO
END DO

indexNyleft = INT(0.0/dy) + 1
indexNyright = INT(Delta_R*0.5/dy) + 1

!indexNyleft = INT((Ly*0.5 - Delta_R)/dy) + 0
!indexNyright = INT((Ly*0.5 + Delta_R)/dy) + 1

para_gauss = 0.035

DO i = 1,indexl1-1
DO j = indexNyleft,indexNyright

 k_eff(i,j) = k !+ EXP(-((j - indexNyleft)*dy-Delta_R*0.5)**2.0/(2*para_gauss*para_gauss))*0.5

END DO
END DO

IF (1.EQ.0) THEN

OPEN(78, IOSTAT=stat, FILE='Temperature2D.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(78, STATUS='delete')
END IF

OPEN(26,FILE='Temperature2D.dat',POSITION='append')

DO j = 1,indexl1
DO m = indexNyleft,indexNyright
WRITE(26,*) (j-1)*dz,(m-1)*dy,w_R(j,m)
END DO
WRITE(26,'(A)')
END DO

PRINT*,'HEY'

CLOSE(26)

STOP

END IF

!END INIT SUBROUTINE

! MEASURING CPU TIME
CALL CPU_TIME(t1)

! DEFINING INVERSE OF COMMONLY USED VARIABLES TO OPTIMIZE SERIAL PERFORMANCE
rdz2 = 1.0/(dz**2.0)

rdy2 = 1.0/(dy**2.0)

rkdz = 1.0/(k*dz)

rkdrdz = 1.0/(k*Delta_R*dz)

rk1 = 1.0/(k1)

rk2 = 1.0/(k2)

avgrho0 = 2.0*rho*rho_snow/(rho_snow+rho)
avgrho1 = 2.0*rho*rho1/(rho+rho1)
avgrho2 = 2.0*rho2*rho1/(rho2+rho1)

avgc1 = 2.0*c*c1/(c+c1)
avgc2 = 2.0*c2*c1/(c2+c1)

avgk1 = 2.0*k*k1/(k+k1)
avgk2 = 2.0*k2*k1/(k2+k1)

DO m = indexNyleft,indexNyright
DO j = indexl1 + 1, indexl2-1 ! indexl1 + 1, indexl2-1
   k_eff(j,m) = k1unfrozen
   c_eff(j,m) = c1unfrozen
END DO
END DO

DO m = indexNyleft,indexNyright
DO j = indexl2 + 1, Nz-1 
   k_eff(j,m) = k2unfrozen
   c_eff(j,m) = c2unfrozen
END DO
END DO

! END OF DEFINING INVERSE

! SUBROUTINE INIT INITIALIZE THE MATRIX Tnew and Told + BOUNDARY CONDITIONS

!nlinescase = 11*10000

avgalpha1 = 2.0*alpha*alpha1/(alpha+alpha1)
avgalpha2 = 2.0*alpha2*alpha1/(alpha2+alpha1)

rribbonyribbonz = REAL(nribbon)/REAL(dy*dz*numberindexribbon)

! ALLOCATE MEMORY FOR Maxsurf AND Tdepth
! Maxsurf is used to track the surface temperature and plot the min and max values at the end
! Tdepth is used to track the temperature at one time step along the depth

listalpha(1) = alpha
listalpha(2) = alpha1
listalpha(3) = alpha2
minalpha = MAXVAL(listalpha)

countdays = 0

!nlinescase = 3600*3

IF(0.EQ.0) THEN
DO i = 1,1500
   !IF((MOD(i,1)).EQ.0)THEN
      IF(countdays.LT.10) THEN
         WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE
         WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      END IF
      countdays = countdays + 1
   !END IF
END DO
END IF

countdays = 0

IF(0.EQ.0) THEN
DO i = 1,1500
   !IF((MOD(i,1)).EQ.0)THEN
      IF(countdays.LT.10) THEN
         WRITE(FileNamebis,'(A,I1.1,A)') 'Hwidth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileNamebis, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileNamebis,'(A,I2.2,A)') 'Hwidth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileNamebis, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileNamebis,'(A,I3.3,A)') 'Hwidth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileNamebis, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE
         WRITE(FileNamebis,'(A,I4.4,A)') 'Hwidth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileNamebis, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      END IF
      countdays = countdays + 1
   !END IF
END DO
END IF

countdays = 0

!STOP

OPEN(77, IOSTAT=stat, FILE='Heatflux.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(77, STATUS='delete')
END IF

OPEN(16,FILE='Heatflux.dat',POSITION='append')

! OPEN PARALLELIZATION REGION
!!$OMP Parallel DEFAULT(__auto)
!!$OMP DO

!power = 10.0

k = 0.5 - 0.5

Delta_R = 0.10 - 0.02

z_R = 0.03 - 0.01

PRINT*,'-----------------------------------'
PRINT*,'------START OF CALCULATIONS--------'
PRINT*,'-----------------------------------'

DO m_loop_1 = 1,7 ! LOOP FOR THERMAL CONDUCTIVITY

k = k + 0.5

Delta_R = 0.10 - 0.02

z_R = 0.03 - 0.01

DO m_loop_2 = 1,6 ! LOOP FOR RIBBONS SPACING

Delta_R = Delta_R + 0.02

z_R = 0.03 - 0.01

DO m_loop_3 = 1,8 ! LOOP FOR RIBBONS EMBEDMENT DEPTH

z_R = z_R + 0.01

IF(countdays.LT.10) THEN
   WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.100) THEN
   WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.1000) THEN
   WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.10000) THEN
   WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.100000) THEN
   WRITE(FileName,'(A,I5.5,A)') 'Tdepth.',countdays,'.dat'
END IF

avgk1 = 2.0*k*k1/(k+k1)

DO i = 1,allocationdepth
DO j = 1,Ny

 Told(i,j) = 4.3*(i-1.0)*dz + 1.7

END DO
END DO

Tnew = Told

DO i = 1,allocationdepth
DO j = 1,Ny

w_R(i,j) = 0.0

END DO
END DO

indice_R = (INT(z_R/dz)) + 1
indexribbony1z = indice_R + 1
indexribbony2z = indice_R + INT(ribbonz/dz) + 1
IF(m_loop_3.EQ.1) THEN
   indexribbony1z = indice_R + 1
ELSE
   indexribbony1z = indice_R + 2
END IF

indexribbony1 = INT(0.0/dy) + 1
indexribbony2 = INT(ribbony*0.5/dy) + 1

DO i = indexribbony1,indexribbony2
   DO m = indexribbony1z,indexribbony2z
     w_R(m,i) = 1.0
   END DO
END DO

indexNyleft = INT(0.0/dy) + 1
indexNyright = INT(Delta_R*0.5/dy) + 1

!OPEN(78, IOSTAT=stat, FILE='Temperature2D.dat', STATUS='old')
!IF (stat.EQ.0) THEN
!CLOSE(78, STATUS='delete')
!END IF

!OPEN(26,FILE='Temperature2D.dat',POSITION='append')

!DO j = 1,indexl1
!DO m = indexNyleft,indexNyright
!WRITE(26,*) (j-1)*dz,(m-1)*dy,w_R(j,m)
!END DO
!WRITE(26,'(A)')
!END DO

!PRINT*,"CHECK FILE"
!READ(*,*)

power_acc = 0.0

! LOOP OVER THE NUMBER OF STEPS
DO i = 1,nlinescase

   hc = 1.78*V_air(i) + 2.21 ! The "2.0" is for air speed 1.78*V_air(i) + 2.21

!!$OMP DO

IF(0.EQ.0) THEN
IF(((MOD(i-1,600)).EQ.0)) THEN
   OPEN(26,FILE=Filename,POSITION='append')
   WRITE(26,*) (i-1)*dt,Tnew(1,indexNyleft),Tnew(1,indexNyright),&
   heatflux_surf(indexNyleft),heatflux_surf(indexNyright),&
   k, Delta_R, z_R, power_acc
   CLOSE(26)
END IF
END IF

IF(Tnew(1,indexNyright).GT.5.0) THEN
   IF(power.GT.5.0) THEN
      power = 0.0
   ELSE
      power = 0.0
   END IF
ELSE
   IF(power.GT.5.0) THEN
      power = 40.0
   ELSE
      power = 40.0
   END IF
   power_acc = power_acc + 1.0
END IF


DO m = indexNyleft+1,indexNyright-1 

   k_eff(1,m) = k
   !k_eff(1,m)  = k + EXP(-((m - indexNyleft)*dy-Delta_R*0.5)**2.0/(2*para_gauss*para_gauss))*0.5
   c_eff(1,m) = c
   d_T_theta(1,m) = 0.0
   
   !Tnew(1,m) = Told(1,m) + ( 2.0*&
   !k_eff(1,m)*dt*rdz2 * (Told(2,m)&
   !- Told(1,m) +  dz/k_eff(1,m) * &
   !(a*f_s(i) + hc*(T_air(i) - Told(1,m)) + eps_s * &
   !(fD(i) - sigma*(Told(1,m) + 273.15)**4.0)  )  ) &
   !+ k_eff(1,m)*dt*rdy2 * (Told(1,m+1)&
   !- 2.0*Told(1,m) + Told(1,m-1) ) ) /&
   !( c_eff(1,m) * rho + latentheat * &
   !d_T_theta(1,m))
   

   Tnew(1,m) = Told(1,m) + ( 2.0*&
   k_eff(1,m)*dt*rdz2 * (Told(2,m)&
   - Told(1,m) +  dz/k_eff(1,m) * &
   (a*0.0 + hc*(T_air(i) - Told(1,m)) + eps_s * &
   (fD(i) - sigma*(Told(1,m) + 273.15)**4.0)  )  ) &
   + k_eff(1,m)*dt*rdy2 * (Told(1,m+1)&
   - 2.0*Told(1,m) + Told(1,m-1) ) ) /&
   ( c_eff(1,m) * rho + latentheat * &
   d_T_theta(1,m))
   
   !Tnew(1,m) = +3.0

   heatflux_surf(m) = -k_eff(1,m)*(Tnew(1,m)-Tnew(2,m))/dz
   heatflux_surf(indexNyleft) = heatflux_surf(indexNyleft+1)
   heatflux_surf(indexNyright) = heatflux_surf(indexNyright-1)

   DO j = 2, indexl1 - 1

      c_eff(j,m) = c
      
      k_eff(j,m) = k
	  !k_eff(j,m) = k + EXP(-((m - indexNyleft)*dy-Delta_R*0.5)**2.0/(2*para_gauss*para_gauss))*0.5

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) + dt*power*w_R(j,m)*rribbonyribbonz )/( c_eff(j,m) * rho )
	  

   END DO
!   !$OMP END DO

!   !$OMP END DO
! CALCULATIONS FOR FIRST INTERFACE

   avgc1 = 2.0*(c_eff(indexl1-1,m)*&
   c_eff(indexl1+1,m))/(c_eff( &
   indexl1 - 1,m)+c_eff(indexl1 + 1,m))

   avgk1 = 2.0*(k_eff(indexl1-1,m)*&
   k_eff(indexl1+1,m))/(k_eff( &
   indexl1 - 1,m)+k_eff(indexl1 + 1,m))

   IF((d_T_theta(indexl1-1,m).EQ.0.0).AND.&
   (d_T_theta(indexl1+1,m).EQ.0.0)) THEN
      avgd_T_theta1 = 0.0
   ELSE
      avgd_T_theta1 = 2.0*(d_T_theta(indexl1-1,m)*&
      d_T_theta(indexl1+1,m))/(d_T_theta(&
      indexl1-1,m)+d_T_theta(indexl1+1,m))
   END IF

   Tnew(indexl1,m) = Told(indexl1,m) &
   + ( dt*rdz2 * ( k_eff(indexl1+1,m)*Told(&
   indexl1+1,m)-(k_eff(indexl1+1,m) + k_eff(&
   indexl1-1,m))*Told(indexl1,m) + k_eff(&
   indexl1-1,m) *Told(indexl1-1,m) ) + avgk1 * dt * rdy2 * &
   ( Told(indexl1,m+1) + Told(indexl1,m-1) &
   -2.0*Told(indexl1,m) ) ) / (avgc1*avgrho1 + latentheat*avgd_T_theta1)
   
   !IF((MOD(i-1,1000).EQ.0)) THEN
   !PRINT*,avgc1,avgk1,avgd_T_theta1,Tnew(indexl1,indexNyleft + 3), i
   !READ(*,*)
   !END IF

! DO LOOP FOR SECOND LAYER

   DO j = indexl1 + 1,indexl2 - 1

      IF(Told(j,m).GE.temp_ref) THEN
         theta(j,m) = theta_0_1
         d_T_theta(j,m) = 0.0
      ELSE
         theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         theta(j,m) = theta_0_1*(1.0 - theta(j,m))
         d_T_theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         d_T_theta(j,m) = theta_0_1*beta*d_T_theta(j,m)/&
                        (temp_ref - Told(j,m))
      END IF
	  
	  IF(d_T_theta(j,m).GT.10.0**(30.0)) THEN
	     d_T_theta(j,m) = 10.0**(30.0)
	  END IF

      exponentlatent = theta(j,m)/theta_0_1

      c_eff(j,m) = c1frozen*(1.0-exponentlatent) + &
      c1unfrozen*exponentlatent
      
      k_eff(j,m) = k1frozen**(1.0-exponentlatent) * &
      k1unfrozen**exponentlatent

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) )/( c_eff(j,m) * rho1 + latentheat*d_T_theta(j,m) )

   END DO

! CALCULATIONS FOR SECOND INTERFACE

   avgc2 = 2.0*(c_eff(indexl2-1,m)*&
   c_eff(indexl2+1,m))/(c_eff( &
   indexl2-1,m)+c_eff(indexl2+1,m))

   avgk2 = 2.0*(k_eff(indexl2-1,m)*&
   k_eff(indexl2+1,m))/(k_eff( &
   indexl2-1,m)+k_eff(indexl2+1,m))

   IF((d_T_theta(indexl2-1,m).EQ.0.0).AND.&
   (d_T_theta(indexl2+1,m).EQ.0.0)) THEN
      avgd_T_theta2 = 0.0
   ELSE
      avgd_T_theta2 = 2.0*(d_T_theta(indexl2-1,m)*&
      d_T_theta(indexl2+1,m))/(d_T_theta(&
      indexl2-1,m)+d_T_theta(indexl2+1,m))
   END IF

   Tnew(indexl2,m) = Told(indexl2,m) &
   + ( dt*rdz2 * ( k_eff(indexl2+1,m)*Told( &
   indexl2+1,m)-  (k_eff(indexl2+1,m) + k_eff( &
   indexl2-1,m))*Told(indexl2,m) + k_eff( &
   indexl2-1,m) * Told(indexl2-1,m) ) + avgk2 * dt * rdy2 * &
   ( Told(indexl2,m+1) + Told(indexl2,m-1) &
   -2.0*Told(indexl2,m) ) ) / (avgc2*avgrho2 + latentheat*avgd_T_theta2)

! DO LOOP FOR THIRD LAYER

   DO j = indexl2+1, Nz-1

      IF(Told(j,m).GE.temp_ref) THEN
         theta(j,m) = theta_0_2
         d_T_theta(j,m) = 0.0
      ELSE
         theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         theta(j,m) = theta_0_2*(1.0 - theta(j,m))
         d_T_theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         d_T_theta(j,m) = theta_0_2*beta*d_T_theta(j,m)/&
                        (temp_ref - Told(j,m))
      END IF
	  
	  IF(d_T_theta(j,m).GT.10.0**(30.0)) THEN
	     d_T_theta(j,m) = 10.0**(30.0)
	  END IF

      exponentlatent = theta(j,m)/theta_0_2

      c_eff(j,m) = c2frozen*(1.0-exponentlatent) + &
      c2unfrozen*exponentlatent
      
      k_eff(j,m) = k2frozen**(1.0-exponentlatent) * &
      k2unfrozen**exponentlatent

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) )/( c_eff(j,m) * rho2 + latentheat*d_T_theta(j,m) )

   END DO

! UPDATE ARRAY 
!   !$OMP WORKSHARE
   
!   !$OMP END WORKSHARE

   Tnew(Nz,m) =&
   Tnew(Nz-1,m)

END DO ! END FOR m (Ny)

DO j = 1,Nz
   Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)
   Told(j,indexNyleft) = Told(j,indexNyleft+1)
END DO

DO j = 1,Nz
   Tnew(j,indexNyright) = Tnew(j,indexNyright-1)
   Told(j,indexNyright) = Told(j,indexNyright-1)
END DO

Told = Tnew

!IF((MOD(i-1,60).EQ.0)) THEN

   !WRITE(16,*) (i-1)*dt,Tnew(indexribbony1z-5,indexNyleft),Tnew(indexribbony1z-4,indexNyleft),&
   !Tnew(indexribbony1z-3,indexNyleft),Tnew(indexribbony1z-2,indexNyleft),&
   !Tnew(indexribbony1z-1,indexNyleft),Tnew(indexribbony1z,indexNyleft),&
   !Tnew(indexribbony1z+1,indexNyleft),Tnew(indexribbony1z+2,indexNyleft)
   
!   WRITE(16,*) (i-1)*dt,T_air(i),V_air(i),fD(i),Tnew(1,indexNyleft+1)
!   PRINT*,i-1

!END IF

   IF(1.EQ.0) THEN
   IF((MOD(i-1,600).EQ.0)) THEN 
      IF(countdays.LT.10) THEN
        WRITE(FileName,'(A,I1.1,A)') 'Twidth.',countdays,'.dat'
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileName,'(A,I2.2,A)') 'Twidth.',countdays,'.dat'
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileName,'(A,I3.3,A)') 'Twidth.',countdays,'.dat'
      ELSE 
         WRITE(FileName,'(A,I4.4,A)') 'Twidth.',countdays,'.dat'
      END IF
      OPEN(36,FILE=Filename,POSITION='append')
      !DO j = 1,(INT(0.1/dz)) + 1
         DO m = indexNyleft,indexNyright
            WRITE(36,*) (m-indexNyleft)*dy,Tnew(1,m)
         END DO
         !WRITE(36,'(A)')
      !END DO
      !IF(countdays.LT.10) THEN
      !  WRITE(FileNamebis,'(A,I1.1,A)') 'Hwidth.',countdays,'.dat'
      !ELSE IF(countdays.LT.100) THEN
      !   WRITE(FileNamebis,'(A,I2.2,A)') 'Hwidth.',countdays,'.dat'
      !ELSE IF(countdays.LT.1000) THEN
      !   WRITE(FileNamebis,'(A,I3.3,A)') 'Hwidth.',countdays,'.dat'
      !ELSE 
      !   WRITE(FileNamebis,'(A,I4.4,A)') 'Hwidth.',countdays,'.dat'
      !END IF
      !OPEN(36,FILE=FileNamebis,POSITION='append')
      !DO j = 1,(INT(0.1/dz)) + 1
      !   DO m = indexNyleft,indexNyright
      !      WRITE(36,*) (m-indexNyleft)*dy,heatflux_surf(m)
      !   END DO
         !WRITE(36,'(A)')
      !END DO
     countdays = countdays + 1
	 !PRINT*,i-1
   END IF
   END IF

END DO ! END FOR nlinescase

countdays = countdays + 1

PRINT*,'-----------------------------------'
PRINT*,'--------SIMULATION NUMBER',countdays,'--------'
PRINT*,'-----------------------------------'

END DO ! END FOR m_loop_1

END DO ! END FOR m_loop_2

END DO ! END FOR m_loop_3

PRINT*,'-----------------------------------'
PRINT*,'--------END OF CALCULATIONS--------'
PRINT*,'-----------------------------------'

!DO m = indexNyleft,indexNyright
!   WRITE(16,*) (m-indexNyleft)*dy,heatflux_surf(m)
!END DO

CLOSE(16)

STOP

!(j-1)*dz,(m-indexNyleft)*dy,Tnew(j,m)
!(m-indexNyleft)*dy,heatflux_surf(m)

!!$OMP END DO
!!$OMP END Parallel

OPEN(28,FILE='Temperature2D.dat',POSITION='append')

DO j = 1,indexsnowinit + indice_R + 10
   DO m = indexNyleft,indexNyright
      WRITE(28,*) (j-1)*dz,(m-1)*dy,Tnew(j,m)
   END DO
   WRITE(28,'(A)')
END DO

CLOSE(28)


! END MAIN PROGRAM
END PROGRAM temp
