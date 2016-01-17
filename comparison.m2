--Documentation---
--r= number of basepoints
--(a,b)=bidegree we want to look at
-- B = basis of the ideal in bidegree (a,b)
-- k= number of elements in the basis B
-- C = matrix of coefficients
-- test ideal of generic forms

gexample=(r,a,b)->(
I=kbpts(r);
B=super basis({a,b},I);
k=numgens source B;
C=random(R^4,R^k);
f=C*transpose B;
test=ideal f;
szf=syz gens test;
Cszf=transpose C * szf;
szB= syz gens ideal B;
)


compzero=(a,b)->(
I=intersect(ideal(s,t),ideal(u,v));
B=super basis({a,b},I);
k=numgens source B;
C=random(R^4,R^k);
f=C*transpose B;
test=ideal f;
szf=syz gens test;
Cszf=transpose C * szf;
szB= syz gens ideal B;
)
