tv={th1,th2}
(* Express the cumulant generating function with cumulants to order 4.*)
Ko[th1_,th2_]=th1^2*sigma1^2/2+th2^2*sigma2^2/2+th1*th2*sigma1*sigma2*rho+
   Sum[tv[[i]]*tv[[j]]*tv[[k]]*kappa[i,j,k]/6,{i,1,2},{j,1,2},{k,1,2}]+
   Sum[tv[[i]]*tv[[j]]*tv[[k]]*tv[[l]]*kappa[i,j,k,l]/24,{i,1,2},{j,1,2},{k,1,2},{l,1,2}]

(* Without loss of generality, take the case with both variances 1.*)
sigma1=sigma2=1
(* Enforce symmetry of the cumulants.*)
kappa[1,2,1]=kappa[1,1,2]
kappa[2,1,1]=kappa[1,1,2]
kappa[2,2,1]=kappa[1,2,2]
kappa[2,1,2]=kappa[1,2,2]
kappa[2,2,1,1]=kappa[1,1,2,2]
kappa[2,1,2,1]=kappa[1,1,2,2]
kappa[2,1,1,2]=kappa[1,1,2,2]
kappa[1,2,2,1]=kappa[1,1,2,2]
kappa[1,2,1,2]=kappa[1,1,2,2]
kappa[2,1,1,1]=kappa[1,1,1,2]
kappa[1,2,1,1]=kappa[1,1,1,2]
kappa[1,1,2,1]=kappa[1,1,1,2]
kappa[2,2,2,1]=kappa[1,1,1,2]
kappa[2,2,1,2]=kappa[1,1,1,2]
kappa[2,1,2,2]=kappa[1,1,1,2]
kappa[1,2,2,2]=kappa[1,1,1,2]
(* Let rtni be the inverse of the square root of the sample size, and let
K2 represent the cumulant generating function of the sum of n=1/rtni^2 
independent copies of the same distribution, multiplied by rtni to keep 
unit variances.*)
K2[th1_,th2_]=ExpandAll[Ko[th1 rtni,th2 rtni] /rtni^2]
(* Calculate the part of the cumulant generating functions with quadratic 
parts removed.*)
Q2[th1_,th2_]=K2[th1,th2]-(th1^2 sigma1^2/2+th2^2 sigma2^2/2+th1 th2 sigma1 sigma2 rho)
(* Exponentiate the remainder to give terms in the Edgeworth expansion
beyond the multivariate Gaussian.*)
P2[th1_,th2_]=Series[Exp[Q2[th1,th2]],{rtni,0,2}]

(* Let Phibar1 and phi be the univariate Gaussian survivor and density functions
respectively.  Define derivatives.*)
Phibar1'[x_]=-phi[x]
Phibar1''[x_]=x phi[x]
phi'[x_]=-x phi[x]
phi''[x_]=(x^2-1) phi[x]

(* Define bivariate versions of what would otherwise be Hermite polynomials.
h[0,0] is the bivariate Gaussian density.  h[-1,-1] is the bivariate survivor
function.  Note here that I am NOT dividing by the density, since for h with
at least one index -1, the quantity before division does not have a factor of
the Gaussian density.  Gaussian distribution quantities will be substituted
in later.*)
h[-1,-1][x_,y_]=Phibar2[x,y,rho]
h[0,-1][x_,y_]=Phibar1[(y-rho x)/Sqrt[1-rho^2]]*phi[x]
h[-1,0][x_,y_]=h[0,-1][y,x]
Do[
   h[j,-1][x_,y_]=D[-h[j-1,-1][x,y],x];
   h[-1,j][x_,y_]=D[-h[-1,j-1][x,y],y];
  ,{j,1,5}]
Do[h[j,k][x_,y_]=D[-h[j-1,k][x,y],x],{j,0,5},{k,0,5}]
jrule=phi[Times[Power[1-rho^2, -1/2],x_-rho y_]]->Sqrt[1-rho^2] phi2[x,y,rho]/phi[y]
Do[h[j,k][x_,y_]=Simplify[h[j,k][x,y]/.jrule],{j,-1,5},{k,-1,5}]

(* Now define derivatives for bivariate quantities.*)
Derivative[1,0,0][phi2][x_,y_,rho_]=-h[1,0][x,y]
Derivative[0,1,0][phi2][x_,y_,rho_]=-h[0,1][x,y]
Derivative[1,0,0][Phibar2][x_,y_,rho_]=-h[0,-1][x,y]
Derivative[2,0,0][Phibar2][x_,y_,rho_]=h[1,-1][x,y]
Derivative[0,1,0][Phibar2][x_,y_,rho_]=-h[-1,0][x,y]
Derivative[0,2,0][Phibar2][x_,y_,rho_]=h[-1,1][x,y]
Derivative[1,1,0][Phibar2][x_,y_,rho_]=h[0,0][x,y]

(* When the h quantities above are substituted into P2, the result will not
be correct, because the arguments of the normal density and distribution
functions need a factor of square root of n.  Also, the result will be less
consise.  So the Gaussian density and tail probability quantities are 
calculated separately in the attached code, and the following rule makes
the appropriate substitution.*)

erules={Phibar2[x1,x2,rho]->pb2,phi2[x1,x2,rho]->p2,
   phi[x1]->px,phi[x2]->py, Phibar1[x1]->pxm,Phibar1[x2]->pym,
   Phibar1[(x1-rho x2)/Sqrt[1-rho^2]]->p1c,
   Phibar1[(x2-rho x1)/Sqrt[1-rho^2]]->p2c}

(* The following set of rules substitutes the h functions into P2.  Rules
need to be applied in order, to make sure, for example, that terms of order
six are hit by the order six rule and not twice by the order three rule.*)
rules[k_]:=Table[th1^j th2^(k-j)->h[j-1,k-j-1][x1,x2],{j,0,k}]
S2[x1_,x2_]=ExpandAll[h[-1,-1][x1,x2]-1+P2[th1,th2]]
R2[x1,x2]=S2[x1,x2]
Do[ R2[x1,x2]=(R2[x1,x2]/.rules[k]),{k,6,2,-1}]
(* Enforce symmetry of the of the Gaussian density.*)
R2[x1,x2]=R2[x1,x2]//.{phi2[x2,x1,rho]->phi2[x1,x2,rho]}
(* Here's the bivariate Edgeworth series.*)
bivariatetail=Simplify[R2[x1,x2]]
univariatetail=Normal[bivariatetail/.{phi[x2]->0,
   Phibar2[x1,x2,rho]->Phibar1[x1],
   phi2[x1,x2,rho]->0,Phibar1[(x2-rho x1)/Sqrt[1-rho^2]]->1}]
(* Calculate the traditional Cornish-Fisher expansion.  Express the ordinate
as a function of sample size (here, as a function of the inverse square 
root of sample size, expand this as a power series in rtni, and invert
the series, as in pp. 47f. of Kolassa (2006).*)
sr1=Solve[Series[univariatetail//.x1->x1f[rtni],{rtni,0,2}]==alpha1,{alpha1,x1f'[0],x1f''[0]}]

(* The bivariate expansion will involve the survivor function for x1
conditional on x2.  The expression will be simpler if we express this using
this ratio to the conditional density.  Call the ratio m.*)
mrule=Phibar1[(x1f[0]-rho x2f[0])/Sqrt[1-rho^2]]-> phi[(1-rho^2)^(-1/2) (x1f[0]-rho x2f[0])] mratio
(* Substitute the univariate Cornish-Fisher expansion into the bivariate tail
probability, and expand in rtni.*)
t2=((Simplify[Series[Normal[bivariatetail]//.{x1->x1f[rtni],x2->x2f[rtni]},
   {rtni,0,2}]/.sr1[[1]]]/.mrule)/.jrule)/.{
      phi2[x2f[0],x1f[0],rho]->phi2[x1f[0],x2f[0],rho]}
(* Equate the resulting expansion to the desired target, and solve.*)
sr2=Simplify[Solve[t2==alpha2,{alpha2,x2f'[0],x2f''[0]}]]
(* Because of symmetry in the cumulants, there's no reason to pass these
as arrays. Make the required elements univariate constants.  Also, simplify
expressions for the leading terms of the ordinates.*)
krule={x1f[0]->x1,x2f[0]->x2,
   kappa[1,1,1]->k111,kappa[1,1,2]->k112,kappa[1,2,2]->k122, 
   kappa[2,2,2]->k222, kappa[1,1,1,1]->k1111, kappa[1,1,1,2]->k1112, 
   kappa[1,1,2,2]->k1122, kappa[1,2,2,2]->k1222, kappa[2,2,2,2]->k2222}
(* Write the important results to the Fortran output.*)
Splice["terms.fm","terms.f",FormatType->FortranForm,PageWidth->62]
