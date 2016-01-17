R=ZZ/31991[s,t,u,v, Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
S=ZZ/31991[X,Y,Z,W]

grobnerb = (M)->(
changeR();
m=gens M;
f=map(R,S,m);
eqgb=kernel f;
)


bgb = (M)->( FR = res M;
      PD = pdim M;
      apply(PD+1, i-> tally(degrees source FR.dd_i))
      )



--extalk=blist(3,2);
T=R**S
I=ideal(t^2*u,t^2*v,s^2*v,s*t*v);
J=ideal(t^2*u,t^2*v,s*t*u+2*s*t*v,s^2*v);
K=ideal(s^2*v+s*t*u,s*t*v+t^2*u,s*t*u+2*t^2*u+t^2*v,t^2*u+t^2*v);
txt=ideal(t^2*u,t^2*v,s*t*v,s^2*v);
expaper=ideal(t^2*u^2+s^2*u*v,t^2*u*v+s^2*v^2,t^2*v^2,s^2*u^2);

--- this code computes the maps of the approximation complex using the induced maps
--- the other codes use less sofisticated commands
zcomplex = (M,a,b)->(
changeT();
cpx=koszul gens M;
L=apply(5,i->kernel cpx.dd_i);
cpx2=koszul gens ideal(X,Y,Z,W);
C=new ChainComplex;
C#0=L_0;
C#1=L_1;
C#2=L_2;
C#3=L_3;
C#4=L_4;
C.dd#0=cpx2.dd_0;
C.dd#1=cpx2.dd_1;
C.dd#2=cpx2.dd_2;
C.dd#3=cpx2.dd_3;
C.dd#4=cpx2.dd_4;
B=apply(5,i-> super basis({2*a-1+i*a,b-1+i*b},C#i));
C2strand=new ChainComplex;
C2strand#0=B_0;
C2strand#1=B_1;
C2strand#2=B_2;
C2strand#3=B_3;
C2strand#4=B_4;
C2strand.dd#0=cpx2.dd_0;
C2strand.dd#1=cpx2.dd_1;
C2strand.dd#2=cpx2.dd_2;
C2strand.dd#3=cpx2.dd_3;
C2strand.dd#4=cpx2.dd_4;
matmaps=apply(5,i-> C2strand.dd_i*C2strand#i);
rels=apply(5,i-> C.dd_i*B_i);
indmaps=apply(4,i->inducedMap(image(B_i),image(B_(i+1)),C.dd_(i+1)));
mds=apply(4,k->map(T^(numgens target indmaps_k),T^(numgens source indmaps_k),(i,j)->
(indmaps_k)_(i,j)));
Cvst=chainComplex mds;
D={};
mns=maxCol(Cvst.dd_1);
D=D|{mns_0};
k=numgens source Cvst.dd_1;
rows=for i from 0 to k-1 list (if not member(i,mns_1) then continue i);
rows=delete(,rows);
D=D|{Cvst.dd_2^rows};
time eqzcpx= (det D_0)/(det D_1); --this line computes the determinant of the complex
)

----SOME DOCUMENTATION FOR THIS CODE
--cpx is the Koszul complex associated to the generators of the ideal
--cpx2 is the koszul complex associated to the sequence (X,Y,Z,W)
--C is the full Z-complex which was computed straight from the koszul
-- complex cpx by taking kernels of the maps
-- C2strand is the approximation complex in which the modules are
-- in the desired degree but the maps are still the ones from the kozxul complex
-- Cvst. gives the desired matrices in the right bidegree, ignore the modules.


Hmatdet = (M,a,b)->(
cpx=koszul gens M;
L=apply(5,i->kernel cpx.dd_i);
cpx2=koszul gens ideal(X,Y,Z,W);
C=new ChainComplex;
C#0=L_0;
C#1=L_1;
C#2=L_2;
C#3=L_3;
C#4=L_4;
C.dd#0=cpx2.dd_0;
C.dd#1=cpx2.dd_1;
C.dd#2=cpx2.dd_2;
C.dd#3=cpx2.dd_3;
C.dd#4=cpx2.dd_4;
B=apply(5,i-> super basis({2*a-1+i*a,b-1+i*b},C#i));
C2strand=new ChainComplex;
C2strand#0=B_0;
C2strand#1=B_1;
C2strand#2=B_2;
C2strand#3=B_3;
C2strand#4=B_4;
C2strand.dd#0=cpx2.dd_0;
C2strand.dd#1=cpx2.dd_1;
C2strand.dd#2=cpx2.dd_2;
C2strand.dd#3=cpx2.dd_3;
C2strand.dd#4=cpx2.dd_4;
matmaps=apply(5,i-> C2strand.dd_i*C2strand#i);
rels=apply(5,i-> C.dd_i*B_i);
indmaps=apply(4,i->inducedMap(image(B_i),image(B_(i+1)),C.dd_(i+1)));
mds=apply(4,k->map(T^(numgens target indmaps_k),T^(numgens source indmaps_k),(i,j)->
(indmaps_k)_(i,j)));
Cvst=chainComplex mds;
)



