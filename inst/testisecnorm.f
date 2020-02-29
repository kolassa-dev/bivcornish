      program test
      double precision x1,targ,rho,x2
      integer efg
      write(6,*) "x1,targ,rho"
      read(5,*) x1,targ,rho
      call isecnorm(x1,x2,targ,rho,efg)
      write(6,'(10(a5,"=",f8.4,1x))') "x1",x1,"rho",rho,"targ",targ,
     x   "x2",x2
      end
