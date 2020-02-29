      subroutine isecnorm(x1,x2,targ,rho,efg)
implicit none
! Solve P[X1>=x1 and X2>=x2]=targ with (X1,X2) multivariate normal mean 0, unit variance, corrlation rho, for x2
      double precision x1,rho,bvnd,targ,temp,x2
!!    double precision dinvnr
      double precision xl,xu,eps,outl,outu,outt
      integer efg,count,maxcount
      logical debug
      debug=.false.
      maxcount=20
      eps=1.0d-4
      temp=eps
!!    call cumnor(x1,temp,xu)
      temp=1.0d0-temp
      if(temp.lt.targ) then
         efg=1
!        write(6,*) "Error",efg,"targ",targ,"x1",x1,"temp",temp
      else 
         efg=0
         xu=5.0; xl=-xu
!!       xu=dinvnr(1.0d0-eps,eps)
!!       xl=dinvnr(eps,1.0d0-eps)
! bvnd gives upper tail
         outl=bvnd(x1,xl,rho)
         outu=bvnd(x1,xu,rho)
         if(outl.lt.targ) efg=2
         if(outu.gt.targ) efg=3
!        write(6,'(10(a5,"=",f8.4,1x))') "targ",targ,"x1",x1, "xl",xl, "outl",outl, "xu",xu, "outu",outu
!        if(efg.gt.0) write(6,*) "Error",efg
         if(efg.eq.0) then
            outt=10.0d0
            count=0
10          if ((abs(outt-targ).gt.1.0d-8).and.(count.lt.maxcount)) then
               count=count+1
               if(((count/2)*2).eq.count) then
!  FortranForm[x/.Solve[{outu==a+b xu,outl==a+b xl,out==a+ b x},{x,a,b}][[1]]]
                  x2=(targ*xl-outu*xl-targ*xu+outl*xu)/(outl-outu)
               else
                  x2=.5d0*(xl+xu)
               end if
               outt=bvnd(x1,x2,rho)
!              write(6,'(10(a5,"=",f8.4,1x))') "x1",x1,"x2",x2,"rho",rho,"outt",outt
               if(outt.lt.targ) then 
                  xu=x2
                  outu=outt
                  if(debug) write(6,*) "In the branch resetting xu"
               end if
               if(outt.gt.targ) then 
                  xl=x2
                  outl=outt
!                 if(debug) write(6,*) "In the branch resetting xl"
               end if
               if(debug) then
!                write(6,'(10(a5,"=",f8.4,1x))') "targ",targ,"x1",&
!                  x1,"xl",xl,"outl",outl,"xu",xu,"outu",outu,&
!                  "x2",x2,"outt",outt
              end if
              goto 10
           end if
!          write(6,*) "count=",count
         else
!          write(6,'(10(a5,"=",f8.4,1x))') "targ",targ,"x1",x1, "xl",xl, "outl",outl, "xu",xu, "outu",outu
         end if
      end if
      return
      end
