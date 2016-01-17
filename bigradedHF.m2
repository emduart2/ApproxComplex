R=ZZ/31991[s,t,u,v, Degrees=>{{1,0},{1,0},{0,1},{0,1}}]

L={{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1}};

------ CODE FROM COMPUTING A BIGRADEDHF FROM THE LIST OF POINTS-----------
--- the function bhf computes the value of the hilbert function of the points at the
--- bidegree (i,j)
bhf=(i,j,L)->(
n=#L; --- number of points in the list
M=super basis({i,j},R); -- number of elements in the basis in degree (i,j)
f= for i from 0 to n-1 list map(R,R,L#i);
E=numgens source M;
-- Now we are going to evaluate each monomial in the basis at a point and save this
-- to be a row of a matrix.
rows={};
i=0;
while i<n
do(
lst= for j from 0 to E-1 list f_i(M_(0,j));
rows=rows|{lst};
i=i+1;
);
mat= matrix rows;
--- the hilbert function of the points in degree i,j is the rank of the matrix mat
hfij=rank mat;
return hfij;
)

hpmatrix=(i,j,L)->(
k=0;
hmat={};
while k<i
do(
rows=for q from 0 to j list bhf(k,q,L);
hmat=hmat|{rows};
k=k+1;
);
return matrix hmat;
)


----- CODE FOR COMPUTING A BIGRADED HF FROM THE IDEAL DEFINING THE POINTS---
----- The function bhfi computes the HF of the ideal at the value (i,j)
bhfq=(i,j,I)->(
M=super basis({i,j},I);
N=super basis({i,j},R);
m=numgens source M;
n=numgens source N;
hfq=n-m;
return hfq;
)

------------
idealbm=(i,j,I)->(
k=0;
hmat={};
while k<i
do(
rows= for q from 0 to j list numgens source basis({k,q},I);
hmat=hmat|{rows};
k=k+1;
);
return matrix hmat;
)



bmatrix=(i,j,I)->(
k=0;
hmat={};
while k<i
do(
rows=for q from 0 to j list bhfq(k,q,I);
hmat=hmat|{rows};
k=k+1;
);
return matrix hmat;
)

----- Code to compute the delta functions for a matrix
delrows=(m)->(
n=rank target m;
del=matrix m_0;
i=0;
while i<n-1
do(
c=m_(i+1)-m_i;
del=del| matrix c;
i=i+1;
);
return del;
)

delM=(m)->(
k=delrows(m);
kk= delrows(transpose k);
return transpose kk;
)

dM=(M)->(
m=rank target M;
n=rank source M;
j=1;
rows={};
while j<n-1
do(
row=for i from 1 to m-1 list (M_(i,j)-M_(i-1,j)-M_(i,j-1)+M_(i-1,j-1));
rows=rows|{row};
j=j+1;
)
)

