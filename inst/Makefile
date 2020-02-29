bivcornish.so: terms.f bvnd.f isecnorm.f90
	R CMD SHLIB  terms.f bvnd.f isecnorm.f90 ; mv terms.so bivcornish.so
testisecnorm.o: testisecnorm.f
	gfortran -c -Wall testisecnorm.f
isecnorm.o: isecnorm.f90
	gfortran -c -Wall isecnorm.f90
isecnormtest.exe: isecnorm.o testisecnorm.o
	gfortran -Wall -o isecnormtest.exe testisecnorm.o isecnorm.o -L /home/kolassa/lib -ldcdf bvnd.f
clean:
	rm -f *.so *.o *.exe
