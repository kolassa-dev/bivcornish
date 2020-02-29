c
c bivnor.f
c Bivariate Normal Distribution
c CACM Algorithm 462
c Thomas G. Donnely (1973), CACM 16: 638
c
c
	double precision function bivnor (ah,ak,R)
	implicit double precision (a-h,o-z)
	gauss(t)=(1.0d0+erf(t/dsqrt(2.0d0)))/2.0d0
	twopi=6.283185307179587d0
	b=0.d0
	idig=15
	gh=gauss(-ah)/2.0d0
	gk=gauss(-ak)/2.0d0
	if (r) 10,30,10
   10	rr=1.0d0-r*r
	if (rr) 20,40,100
   20	continue
c       write(3,99999) r
99999	format(' bivnor r is',d26.16)
c       stop
   30	b=4.0d0*gh*gk
	go to 350
   40	if (r) 50,70,70
   50	if (ah+ak) 60,350,350
   60	b=2.0d0*(gh+gk)-1.0d0
	go to 350
   70	if (ah-ak) 80,90,90
   80	b=2.0d0*gk
	go to 350
   90	b=2.0d0*gh
	go to 350
  100	sqr=dsqrt(rr)
	if (idig-15) 120,110,120
  110	con=twopi*1.d-15/2.0d0
	go to 140
  120	con=twopi/2.0d0
	do 130 i=1,idig
	   con=con/10.0d0
  130	continue
  140	if (ah) 170,150,170
  150	if (ak) 190,160,190
  160	b=datan(r/sqr)/twopi + 0.25d0
	go to 350
  170 	b=gh
	if (ah*ak) 180,200,190
  180	b=b-0.5d0
  190	b=b+gk
	if (ah) 200,340,200
  200	wh=-ah
	wk=(ak/ah-r)/sqr
	gw=2.0d0*gh
	is=-1
  210	sgn=-1.0d0
	t=0.0d0
	if (wk) 220,320,220
  220	if (dabs(wk)-1.0d0)270,230,240
  230	t = wk*gw*(1.0d0-gw)/2.0d0
	go to 310
  240	sgn=-sgn
	wh=wh*wk
	g2=gauss(wh)
	wk=1.0d0/wk
	if (wk) 250,260,260
  250	b=b+0.5d0
  260	b=b-(gw+g2)/2.d0 +gw*g2
  270	h2=wh*wh
	a2=wk*wk
	h4=h2/2.0d0
	ex=dexp(-h4)
	w2=h4*ex
	ap=1.0d0
	s2=ap-ex
	sp=ap
	s1=0.0d0
	sn=s1
	conex=dabs(con/wk)
	goto 290
  280	sn=sp
	sp=sp+1.0d0
	s2=s2-w2
	w2=w2*h4/sp
	ap=-ap*a2
  290	cn=ap*s2/(sn+sp)
	s1=s1+cn
	if (dabs(cn)-conex) 300,300,280
  300	t=(datan(wk)-wk*s1)/twopi
  310	b=b+sgn*t
c       write(6,*) b
  320	if (is) 330,350,350
  330	if (ak) 340,350,340
  340	wh=-ak
	wk=(ah/ak-r)/sqr
	gw=2.0d0*gk
	is=1
	go to 210
  350	if (b) 360,370,370
  360	b=0.0d0
  370	if (b-1.0d0) 390,390,380
  380	b=1.0d0
  390	bivnor=b
	return
	end
        subroutine pbivnor(p,r,ah,ak)
        double precision p,r,ah,ak,bivnor
        external bivnor
        p=bivnor(ah,ak,r)
        return
        end
