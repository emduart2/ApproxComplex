R=ZZ/31991[s,t,u,v, Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
S=ZZ/31991[X,Y,Z,W]



----- this code generates a random ideal with a linearsyzygy
blin=(i,j)->(
M=super basis({i,j},R);
n=numgens source M;
ges= apply(2,i->(random(R^2,R^n)*transpose M)_(i,0));
p2=ges#0;
p3=ges#1;
Mp=super basis({i,j-1},R);
np=numgens source Mp;
p= random(R^1,R^np)*transpose Mp;
test=ideal(p*u,p*v,p2,p3);
)


------ Then function uplusv takes a polynomial f and returns
------ a tuple (f2,g2) with the property that f=f2*u+g2*v

uplusv = f->(
f2=contract(u,f);
G=f-u*f2;
g2=contract(v,G);
return (f2,g2);
)


---- linsyz computes the equation of the image 
---- knowing that we have a linear syzygy and that
---- the space U is basepoint free we truncate the computation of
---- the syzygies in the desired degree

linsyz = (M,a,b)->(
changeR();
time SZ = syz(gens M,DegreeLimit=>{2*a+2*b});
time bzz= super basis({3*a-1,2*b-1}, image SZ);
T=R**S;
g=map(T,R,matrix{{s,t,u,v}});
f=map(T,S,matrix{{X,Y,Z,W}});
rels=f(vars S)*g(bzz);
N=super basis({2*a-1,b-1},R);
lst=contract(g(N),rels);
rws=apply(rank source lst,i->lst_(0,i));
almost=pack(rank source rels,rws);
time MAT=matrix almost;
time eqtsyz = det MAT;
)


---- The code optimus uses the linear syzygy to compute
----- the other two syzygies of the desired degree, we use several contract
----- commands to extract coefficients and write the polynomials
----- as in the proof of the theorem
optimus = (M,a,b)->(
changeR();
time linSZ = syz(gens M,DegreeLimit=>{a+b+2});
A=contract(u,linSZ);
B=contract(v,linSZ);
C=gens M;
pv=C*A;
pu=C*B;
p=contract(v,pv);
--C_(0,0)=p*u;
--C_(0,1)=p*v;
newgens=matrix{{p*u,p*v,C_(0,2),C_(0,3)}};
s1=uplusv(C_(0,2));
s2=uplusv(C_(0,3));
SZ=matrix{{-v,s1_0,s2_0},{u,s1_1,s2_1},{0,-p,0},{0,0,-p}};
time bzz= super basis({2*a-1,b-1}, image SZ);
T=R**S;
g=map(T,R,matrix{{s,t,u,v}});
f=map(T,S,matrix{{X,Y,Z,W}});
rels=f(vars S)*g(bzz);
N=super basis({2*a-1,b-1},R);
lst=contract(g(N),rels);
rws=apply(rank source lst,i->lst_(0,i));
almost=pack(rank source rels,rws);
time MAT=matrix almost;
--time eqtopt = det MAT;
)


----this syzgs code seems to be the same thing as the linsyz code
---- not sure what it did at the moment

syzgs = (M,a,b)->(
changeR();
SZ = syz gens M;
bzz= super basis({3*a-1,2*b-1}, image SZ);
T=R**S;
g=map(T,R,matrix{{s,t,u,v}});
f=map(T,S,matrix{{X,Y,Z,W}});
rels=f(vars S)*g(bzz);
N=super basis({2*a-1,b-1},R);
lst=contract(g(N),rels);
rws=apply(rank source lst,i->lst_(0,i));
almost=pack(rank source rels,rws);
MAT=matrix almost;
pwr=det MAT;
eqtt = radical ideal(pwr);
--use S;
--eqtt=substitute(eqtt,S)
)



