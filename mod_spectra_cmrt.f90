Module mod_spectra_cmrt
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J


!! Output
complex*16 dip_mom_corr

!! 2d spectra
real*8 t1,t2,t3,dt_t1
real*8 t1_max,t3_max
integer nst_t1,nst_t3,nst_avg
complex*16,allocatable :: mu_t3(:,:,:)

!! System specific
integer nquant,nsite,ndou_ex
real*8 gamma,eta,temperature
complex*16,allocatable :: dip_moment_pos(:,:),dip_moment_neg(:,:),Hamil_e(:,:)
real*8,allocatable :: Hamil_diab(:,:),lambda(:)
complex*16,allocatable::dou_ex(:,:),Hamil_s(:,:),rho_cmrt(:,:,:)


!! Evolution
integer nsteps
real*8 dt,tim_tot,curr_time(3,6)

!! Parallelization
integer iparallel,iflow,iwait,ifolder,iproc

!! Misc
real*8 tim_ind

!! CMRT
real*8,allocatable:: H_diab(:,:),E_exc(:),lambda_diab(:,:)
real*8,allocatable::lambda_exc(:,:,:,:),c_tr(:,:)
real*8,allocatable::RR(:,:,:,:)
complex*16,allocatable::gg(:,:,:,:),gdot(:,:,:,:),g2dot(:,:,:,:)
complex*16,allocatable::g_diab(:,:),gdot_diab(:,:),g2dot_diab(:,:)
complex*16,allocatable::sigma(:,:),sigma_diab(:,:)
real*8,allocatable:: q_mat(:,:,:)

!!diagonalization
integer nold
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)


contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  integer i,j
  character st_ch

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar

  open(10,file="spectra_cmrt.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) nsite
  read(10,*) t2
  read(10,*) t1_max
  read(10,*) t3_max
  read(10,*) nst_avg
  read(10,*) dt
  read(10,*) dt_t1
  read(10,*) tim_tot
  allocate(Hamil_s(nsite,nsite))
  read(10,*)
  read(10,*)
  do i=1,nsite
    read(10,*) Hamil_s(i,:)
  enddo
  read(10,*)
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  !---------------------------------------------------------- 
!!ndou_ex=no. of double excitations!!
!!nquant=no. of total states in the dynamics!!

ndou_ex=nsite*(nsite-1)/2
nquant=1+nsite+ndou_ex

  nsteps=nint(tim_tot/dt)
  nst_t1=nint(t1_max/dt_t1)
  nst_t3=nint(t3_max/dt)


  i=nst_t1/nst_avg
  j=nst_t3/nst_avg

  allocate(dip_moment_pos(nquant,nquant),dip_moment_neg(nquant,nquant),Hamil_e(nquant,nquant))
  allocate(mu_t3(nquant,nquant,6))
  allocate(rho_cmrt(nquant,nquant,6))
  allocate(Hamil_diab(nquant,nquant))

  allocate(H_diab(nquant,nquant),E_exc(nquant),lambda_diab(nquant,nquant),lambda(nquant))
  allocate(lambda_exc(nquant,nquant,nquant,nquant),c_tr(nquant,nquant))
  allocate(gg(nquant,nquant,nquant,nquant),gdot(nquant,nquant,nquant,nquant),g2dot(nquant,nquant,nquant,nquant))
  allocate(g_diab(nquant,nquant),gdot_diab(nquant,nquant),g2dot_diab(nquant,nquant))
  allocate(sigma(nquant,nquant),RR(nquant,nquant,3,6),sigma_diab(nquant,nquant))
  allocate(q_mat(1+nsite,nquant,nquant))
  allocate(dou_ex(ndou_ex,ndou_ex))
!  allocate(Hamil_s(nsite,nsite))

  RR=0d0
  curr_time=0d0

!!q_mat are the coupling matrix of system to bath!!

q_mat=0d0

q_mat(2,2,2)=1d0
q_mat(2,5,5)=1d0
q_mat(2,6,6)=1d0
 
q_mat(3,3,3)=1d0
q_mat(3,5,5)=1d0
q_mat(3,7,7)=1d0

q_mat(4,4,4)=1d0
q_mat(4,6,6)=1d0
q_mat(4,7,7)=1d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer clck_counts_beg, clck_counts_end, clck_rate
  real*8 t1_cpu,t2_cpu
  integer i,j,nst_t2
  complex*16 rho_save(nquant,nquant,6)

  call system_clock ( clck_counts_beg, clck_rate )
  call cpu_time(t1_cpu)

  open(101,file="response.out")
  open(103,file="se.out")
  open(104,file="gsb.out")
  open(105,file="esa.out")
  call setup_parameters
  call init_cond
 call strike_lightning(0)  !! t=0
  rho_save=rho_cmrt
  t1=0.d0
  !! parallelizing over t1
  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    t1_max=t1_max/real(iparallel)
    nst_t1=nint(t1_max/dt_t1)
    if(ifolder>1) then
      do i=1,ifolder-1
        call evolve_cmrt(rho_cmrt,nst_t1,dt_t1,1)
        t1=t1+t1_max
      enddo
    endif
  endif

  rho_save=rho_cmrt

  do i=1,nst_t1/nst_avg  !! evolve for t1
   
      call strike_lightning(1)  !! second laser hit
    nst_t2=nint(t2/dt)
curr_time(2,:)=0d0
RR(:,:,2,:)=0d0
    call evolve_cmrt(rho_cmrt,nst_t2,dt_t1,2)  !! evolve for t2
    
    call strike_lightning(2)  !! third laser hit
curr_time(3,:)=0d0
RR(:,:,3,:)=0d0
   
    do j=1,nst_t3 !! evolve for t3
      t3=(j-1)*t3_max/real(nst_t3-1)
      mu_t3=rho_cmrt
      if(mod(i,nst_avg)==1.and.mod(j,nst_avg)==1) call compute_response(i,j) !! measurement
     call evolve_cmrt(rho_cmrt,1,dt,3)
    enddo

    if(mod(i,nst_avg)==1)write(101,*)
    if(mod(i,nst_avg)==1)write(103,*)
    if(mod(i,nst_avg)==1)write(104,*)
    if(mod(i,nst_avg)==1)write(105,*)

    rho_cmrt=rho_save

    call evolve_cmrt(rho_cmrt,1,dt_t1,1) !! evolve for t1

    t1=t1+dt_t1

    rho_save=rho_cmrt

  enddo
  close(101)
  close(103)
  close(104)
  close(105)

  call system_clock ( clck_counts_end, clck_rate )
  call cpu_time(t2_cpu)
  write(6,*) "total time=",(clck_counts_end - clck_counts_beg) / real(clck_rate),t2_cpu-t1_cpu
  write(6,*) "tim_index=",tim_ind



end subroutine main
!---------------------------------------------------------- 

subroutine setup_parameters
  implicit none
  integer i,j,k
  real*8 tmp,omg,cc
  real*8:: en(nquant),eve(nquant,nquant),dipole_e(nsite)

  lambda=0.d0
  
  gamma=1.d0/50.d-15!*au2s
  lambda(2:1+nsite)=35.d0*wave_to_J/pi
  !lambda(2+nsite:nquant)=2d0*35.d0*wave_to_J/pi
  temperature=77d0

call dipole_matrix_upper
call hamiltonian_full 

!do i=1,nquant
!do j=1,nquant
!  write(1,*) i,j,(Hamil_diab(i,j))
!  write(1,*) i,j,real(dip_moment_neg(i,j))
!enddo
!enddo
!stop

Hamil_diab=Hamil_diab*wave_to_J

 call diag(Hamil_diab,nquant,en,eve,nquant)
  c_tr = eve
  E_exc = en

  lambda_diab=0.d0
  do i=1,nquant
    lambda_diab(i,i)=lambda(i)*pi
  enddo


write(6,*) "Parameters Set ..."

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i

  rho_cmrt=0.d0 
  rho_cmrt(1,1,:)=1.d0

  curr_time=0.d0
  dip_mom_corr=0.d0

  tim_ind=0.d0

write(6,*) "Intitial Conditions set ... "

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine strike_lightning(iflag)
  implicit none
  integer,intent(in) :: iflag
  integer i

 !! dipole moment on the right --> 0
 !! dipole moment on the left  --> 1

  if(iflag==0) then   !! t=0
    call mult(rho_cmrt(:,:,1),dip_moment_neg,0)
    call mult(rho_cmrt(:,:,2),dip_moment_neg,0)
    call mult(rho_cmrt(:,:,3),dip_moment_neg,0)
    call mult(rho_cmrt(:,:,4),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,5),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,6),dip_moment_pos,1)
  endif

  if(iflag==1) then   !! t=t1
    call mult(rho_cmrt(:,:,1),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,2),dip_moment_pos,0)
    call mult(rho_cmrt(:,:,3),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,4),dip_moment_neg,0)
    call mult(rho_cmrt(:,:,5),dip_moment_neg,1)
    call mult(rho_cmrt(:,:,6),dip_moment_neg,0)
  endif

  if(iflag==2) then   !! t=t2
    call mult(rho_cmrt(:,:,1),dip_moment_pos,0)
    call mult(rho_cmrt(:,:,2),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,3),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,4),dip_moment_pos,0)
    call mult(rho_cmrt(:,:,5),dip_moment_pos,1)
    call mult(rho_cmrt(:,:,6),dip_moment_pos,1)
  endif

end subroutine strike_lightning
!-----------------------------------------------------------------  

subroutine mult(rho,mu,iflag)
  implicit none
  complex*16,intent(inout)::rho(nquant,nquant)
  complex*16,intent(in)::mu(nquant,nquant)
  integer,intent(in)::iflag
  integer i

  if(iflag==0) then !!  dipole moment on the right
      rho(:,:)=matmul(rho(:,:),mu)
  endif

  if(iflag==1) then !!  dipole moment on the left
      rho(:,:)=matmul(mu,rho(:,:))
  endif

end subroutine mult
!-----------------------------------------------------------------  

subroutine compute_response(i,j)
  implicit none
  integer,intent(in) :: i,j
  integer k,l
  complex*16 mat(nquant,nquant),response(6)

!! 6 pathways!!
  do l=1,6
    mat=matmul(dip_moment_neg,mu_t3(:,:,l))

    response(l)=0.d0
    do k=1,nquant
      response(l)=response(l)+mat(k,k)
    enddo
  enddo

write(101,'(T10, E12.4, E12.4, F12.7, F12.7, F12.7, F12.7)') t1/au2s,t3/au2s,response(1)+response(2)-response(3),response(4)+response(5)-response(6)
write(103,'(T10, E12.4, E12.4, F12.7, F12.7, F12.7, F12.7)') t1/au2s,t3/au2s,response(1),response(4)
write(104,'(T10, E12.4, E12.4, F12.7, F12.7, F12.7, F12.7)') t1/au2s,t3/au2s,response(2),response(5)
write(105,'(T10, E12.4, E12.4, F12.7, F12.7, F12.7, F12.7)') t1/au2s,t3/au2s,response(3),response(6)

end subroutine compute_response
!-----------------------------------------------------------------  

function commute(mat1,mat2) result(mat3)
  complex*16,intent(in) :: mat1(:,:),mat2(:,:)
  complex*16,allocatable :: mat3(:,:)
  integer i1

  i1=size(mat1,1)
  allocate(mat3(i1,i1))

  mat3=matmul(mat1,mat2)-matmul(mat2,mat1)

end function commute
!-----------------------------------------------------------------  

integer function factorial(n)
  implicit none
  integer,intent(in) :: n
  integer i

  if(n==0) then
    factorial=1
  else
    factorial=1
    do i=1,n
      factorial=factorial*i
    enddo
  endif

end function factorial
!-----------------------------------------------------------------  
subroutine evolve_cmrt(rho_diab,nsteps,dt,iflag)
  implicit none
  complex*16,intent(inout)::rho_diab(nquant,nquant,6)
  complex*16 :: rho(nquant,nquant,6)
  integer,intent(in)::nsteps
  real*8,intent(in)::dt
  complex*16,dimension(nquant,nquant) :: vec,vec_diab
  integer i,j
  integer,intent(in)::iflag

!! 6 pathways !!
  do i=1,6
  do j=1,nsteps
  vec_diab(:,:)=rho_diab(:,:,i)

  vec=matmul(transpose(c_tr),matmul(vec_diab,c_tr)) 
  rho(:,:,i)=run_cmrt(vec,dt,iflag,i)
  rho_diab(:,:,i)=matmul(c_tr,matmul(rho(:,:,i),transpose(c_tr)))
  enddo
  enddo
end subroutine evolve_cmrt
!-----------------------------------------------------------------  

function run_cmrt(initialDensityMatrix, dt,iflag,i) result(A)
  complex*16, intent(inout) :: initialDensityMatrix(nquant,nquant)
  real*8, intent(in) :: dt
  integer ::  j
  complex*16 sigma_dot(nquant,nquant,3,6),A(nquant,nquant)
  integer,intent(in)::iflag,i

!write(1,*) iflag,curr_time(1,1),curr_time(3,1)

sigma=initialDensityMatrix

  call calculate_gg(iflag,i)
  call calculate_sigmadot(sigma_dot,iflag,i)

  sigma=sigma+sigma_dot(:,:,iflag,i)*dt


  A=sigma
  curr_time(iflag,i)=curr_time(iflag,i)+dt
!if (iflag==1) then 
!curr_time(iflag,i)=curr_time(iflag,i)+curr_time(3,i)
!endif
end function run_cmrt
!-----------------------------------------------------------------  

subroutine calculate_gg(iflag,w)
  implicit none
  integer i,j,k,l,n,m
  integer,intent(in)::iflag,w

  call calculate_g_diab(iflag,w)


  do i=1,nquant
    do j=1,nquant
      do k=1,nquant
        do l=1,nquant
          gg(i,j,k,l)=0.d0
          gdot(i,j,k,l)=0.d0
          g2dot(i,j,k,l)=0.d0
          lambda_Exc(i,j,k,l)=0.d0
          do n=1,1+nsite
              lambda_exc(i,j,k,l)=lambda_exc(i,j,k,l)+dot_product(c_tr(:,i),matmul(q_mat(n,:,:),c_tr(:,j)))*dot_product(c_tr(:,k),matmul(q_mat(n,:,:),c_tr(:,l)))*lambda_diab(n,n)
            do m=1,1+nsite
              gg(i,j,k,l)=gg(i,j,k,l)+dot_product(c_tr(:,i),matmul(q_mat(n,:,:),c_tr(:,j)))*dot_product(c_tr(:,k),matmul(q_mat(m,:,:),c_tr(:,l)))*g_diab(n,m)
              gdot(i,j,k,l)=gdot(i,j,k,l)+dot_product(c_tr(:,i),matmul(q_mat(n,:,:),c_tr(:,j)))*dot_product(c_tr(:,k),matmul(q_mat(m,:,:),c_tr(:,l)))*gdot_diab(n,m)
              g2dot(i,j,k,l)=g2dot(i,j,k,l)+dot_product(c_tr(:,i),matmul(q_mat(n,:,:),c_tr(:,j)))*dot_product(c_tr(:,k),matmul(q_mat(m,:,:),c_tr(:,l)))*g2dot_diab(n,m)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine calculate_gg
!-----------------------------------------------------------------  

subroutine calculate_g_diab(iflag,w)
  implicit none
  integer i,j,k
  integer,parameter::nw=10000
  real*8 dw,w_max
  real*8 omg(nw),wt(nw),spec(nquant,nquant,nw)
  integer,intent(in)::iflag,w

  w_max=100*gamma
  dw=w_max/real(nw)
  do i=1,nw
    omg(i)=i*dw
    do j=1,nquant
      do k=1,nquant
        spec(j,k,i)=spectral(j,k,omg(i))
!write(20,*) omg(i)/(1*pi*clight),spec(1,1,i)/wave_to_J
      enddo
    enddo
  enddo
!  stop
  wt=omg*curr_time(iflag,w)

  do i=1,nquant
    do j=1,nquant
      g_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(1-cos(wt))+iota*(sin(wt)-wt)))
      g_diab(i,j)=g_diab(i,j)/(hbar)
      gdot_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(sin(wt)*omg)+iota*(cos(wt)*omg-omg)))
      gdot_diab(i,j)=gdot_diab(i,j)/(hbar)
      g2dot_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(cos(wt)*omg*omg)-iota*(sin(wt)*omg*omg)))
      g2dot_diab(i,j)=g2dot_diab(i,j)/(hbar)
    enddo
  enddo



end subroutine calculate_g_diab
!-----------------------------------------------------------------  

real*8 function spectral(i,j,w)
  implicit none
  integer,intent(in)::i,j
  real*8,intent(in)::w

  if(i==j) then
    spectral=2*lambda(i) * gamma*w/(w**2+gamma**2)
  else
    !spectral=2*lambda * gamma*w/(w**2+gamma**2)
    spectral=0.d0
  endif


end function spectral
!-----------------------------------------------------------------  

subroutine calculate_sigmadot(sigma_dot,iflag,w)
  implicit none
  complex*16,intent(out):: sigma_dot(nquant,nquant,3,6)
  complex*16 RR_pd(nquant,nquant)
  integer i,j,f
  integer,intent(in):: iflag,w
 
  call calculate_RR(iflag,w)
  call calculate_RR_pd(RR_pd)

  do i=1,nquant
    do j=1,nquant
      sigma_dot(i,j,iflag,w)=-iota*(E_exc(i)-E_exc(j))/hbar * sigma(i,j)
      if(i==j) then
        do f=1,nquant
          sigma_dot(i,j,iflag,w)=sigma_dot(i,j,iflag,w)+RR(i,f,iflag,w)*sigma(f,f)-RR(f,i,iflag,w)*sigma(i,i)
        enddo
      endif
    sigma_dot(i,j,iflag,w)=sigma_dot(i,j,iflag,w)-RR_pd(i,j)*sigma(i,j)
      if(i.ne.j) then
        do f=1,nquant
          sigma_dot(i,j,iflag,w)=sigma_dot(i,j,iflag,w)-0.5*(RR(f,i,iflag,w)+RR(f,j,iflag,w))*sigma(i,j)
        enddo
      endif
    enddo
  enddo


end subroutine calculate_sigmadot
!-----------------------------------------------------------------  

subroutine calculate_RR(iflag,w)
  implicit none
  complex*16 AA(nquant),FF(nquant),XX(nquant,nquant),fac
  integer i,j
  real*8 tt
  integer,intent(in)::iflag,w

  tt=curr_time(iflag,w)

  do i=1,nquant
    AA(i) = exp(-iota*E_exc(i)*tt/hbar-gg(i,i,i,i))
    FF(i)=exp(-iota*(E_exc(i)-2*lambda_exc(i,i,i,i))*tt/hbar-conjg(gg(i,i,i,i)))
  enddo

  do i=1,nquant
    do j=1,nquant
      XX(i,j)=exp(2*(gg(i,i,j,j)+iota*lambda_exc(i,i,j,j)*tt/hbar))
      fac=(gdot(j,i,i,i)-gdot(j,i,j,j)-2*iota*lambda_exc(j,i,j,j)/hbar)
      fac=fac*(gdot(i,j,i,i)-gdot(i,j,j,j)-2*iota*lambda_exc(i,j,j,j)/hbar)
      XX(i,j)=XX(i,j)*(g2dot(j,i,i,j)-fac)
    enddo
  enddo

  !RR=0.d0

  do i=1,nquant
    do j=1,nquant
      if(i.ne.j)RR(i,j,iflag,w)=RR(i,j,iflag,w)+2*dt*real(conjg(FF(j))*AA(i)*XX(i,j))
    enddo
  enddo

end subroutine calculate_RR
!-----------------------------------------------------------------  
subroutine calculate_RR_pd(RR_pd)
  implicit none
  complex*16,intent(out)::RR_pd(nquant,nquant)
  integer i,j

  RR_pd=0.d0
  do i=1,nquant
    do j=1,nquant
      if(i.ne.j) then
        RR_pd(i,j) = real(gdot(i,i,i,i)+gdot(j,j,j,j)-2*gdot(i,i,j,j))
        RR_pd(i,j)=RR_pd(i,j)+iota*aimag(gdot(i,i,i,i)-gdot(j,j,j,j))
      endif
    enddo
  enddo

end subroutine calculate_RR_pd
!-----------------------------------------------------------------  
!-----------------------------------------------------------------  
subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and
  !eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are
  !re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 
subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and
  !eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are
  !re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  !real*8,intent(in) :: mat(n,n)
  complex*16,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO) 
   lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif


  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 
subroutine dipole_matrix_upper
  implicit none
  complex*16, allocatable :: mu_E(:)
  integer, allocatable :: double_pairs(:,:)
  integer :: i, j, idx, di, dj

  ! Allocate arrays
  allocate(mu_E(nsite))
  allocate(double_pairs(2, ndou_ex))

  ! initialize mu_E (transition dipole strengths)

mu_E = (/1.d0, 0.391d0, -0.312d0/)

  ! Initialize dipole matrix to zero
 dip_moment_neg = 0.0d0

  ! Ground ↔ Single excitations (only upper triangle)
  do i = 1,nsite
     dip_moment_neg(1, i+1) = mu_E(i)
  end do
  ! Generate all (i < j) double excitation pairs and store them
  idx = 1
  do i = 1, nsite-1
     do j = i+1,nsite
        double_pairs(1, idx) = i
        double_pairs(2, idx) = j
        idx = idx + 1
     end do
  end do

  ! Single ↔ Double excitations (only upper triangle)
  do idx = 1,ndou_ex 
     di = double_pairs(1, idx)  ! site i
     dj = double_pairs(2, idx)  ! site j

     ! Indices in mu:
     ! single excitation on site k → k+1
     ! double excitation (i,j) → idx + 1 + nsite
     if (di+1 < idx+1+nsite) dip_moment_neg(di+1, idx+1+nsite) = mu_E(dj)
     if (dj+1 < idx+1+nsite) dip_moment_neg(dj+1, idx+1+nsite) = mu_E(di)
     ! skip lower triangle
  end do

!do i=1,nquant
!do j=1,nquant
!write(*,*) i,j,dip_moment_neg(i,j)
!enddo
!enddo
!stop

dip_moment_pos=transpose(dip_moment_neg)
  ! Clean up
  deallocate(mu_E, double_pairs)
end subroutine dipole_matrix_upper

!-----------------------------------------------------------------------------
subroutine hamiltonian_full
  implicit none

  integer :: i, j, idx1, idx2
  integer :: d1_i, d1_j, d2_i, d2_j
  integer, allocatable :: double_pairs(:,:)

  ! Allocate pair list for double excitations
  allocate(double_pairs(2, ndou_ex))

  ! Zero out the full Hamiltonian
  Hamil_diab(:,:) = 0.0d0

  ! Ground state energy
  Hamil_diab(1,1) = 0.0d0

  ! --- Single-excitation block ---
  do i = 1, nsite
     do j = i, nsite
        Hamil_diab(i+1, j+1) = Hamil_s(i, j)
     end do
  end do

  ! --- Build double excitation index list ---
  idx1 = 1
  do i = 1, nsite - 1
     do j = i + 1, nsite
        double_pairs(1, idx1) = i
        double_pairs(2, idx1) = j
        idx1 = idx1 + 1
     end do
  end do

  ! --- Double-excitation block ---
  do idx1 = 1, ndou_ex
     d1_i = double_pairs(1, idx1)
     d1_j = double_pairs(2, idx1)

     do idx2 = idx1, ndou_ex
        d2_i = double_pairs(1, idx2)
        d2_j = double_pairs(2, idx2)

      if (idx1 == idx2) then
           ! Diagonal: sum of single-site energies
     Hamil_diab(idx1+1+nsite, idx1+1+nsite) = Hamil_s(d1_i, d1_i)+Hamil_s(d1_j, d1_j)
        else
           ! Off-diagonal: coupling terms between double excitations
           Hamil_diab(idx1+1+nsite, idx2+1+nsite) = 0.0d0

     if (d1_i == d2_i .and. d1_j /= d2_j) Hamil_diab(idx1+1+nsite,idx2+1+nsite) = Hamil_diab(idx1+1+nsite, idx2+1+nsite) + Hamil_s(d1_j, d2_j)
     if (d1_i == d2_j .and. d1_j /= d2_i) Hamil_diab(idx1+1+nsite,idx2+1+nsite) = Hamil_diab(idx1+1+nsite, idx2+1+nsite) + Hamil_s(d1_j, d2_i)
     if (d1_j == d2_i .and. d1_i /= d2_j) Hamil_diab(idx1+1+nsite,idx2+1+nsite) = Hamil_diab(idx1+1+nsite, idx2+1+nsite) + Hamil_s(d1_i, d2_j)
     if (d1_j == d2_j .and. d1_i /= d2_i) Hamil_diab(idx1+1+nsite,idx2+1+nsite) = Hamil_diab(idx1+1+nsite, idx2+1+nsite) + Hamil_s(d1_i, d2_i)
        end if
     end do
  end do

do i=1,nquant
  do j = i+1, nquant
    Hamil_diab(j,i) = Hamil_diab(i,j)
  enddo
enddo

!do i=1,nquant
!do j=1,nquant
!write(*,*) i,j,real(Hamil_diab(i,j))
!enddo
!enddo
  ! Clean up
  deallocate(double_pairs)
!stop

end subroutine hamiltonian_full

!-----------------------------------------------------------------------------

End Module mod_spectra_cmrt
