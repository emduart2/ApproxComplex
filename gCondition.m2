
installPackage "Depth"
installPackage "EliminationMatrices"
load "/home/eliana/Dropbox/UIUC/Research/BasePoints/macau2/zcomplex.m2"
load "/home/eliana/Dropbox/UIUC/Research/BasePoints/macau2/examples.m2"
load "/home/eliana/Dropbox/UIUC/Research/BasePoints/macau2/comparison.m2"
load "/home/eliana/Dropbox/UIUC/Research/BasePoints/macau2/reesexample.m2"
load "/home/eliana/Dropbox/UIUC/Research/BasePoints/macau2/bigradedHF.m2"


use R

I=ideal super basis({2,1},kbpts(3))

I=ideal super basis({2,1},kbpts(3)) --good
I=ideal super basis({3,1},kbpts(3)) -- bad

I=ideal super basis({3,1},kbpts(4))-- good
I=ideal super basis({4,1},kbpts(4))-- bad

I=ideal super basis({4,1},kbpts(5))-- good

I=ideal super basis({5,1},kbpts(5))-- bad

I=ideal super basis({9,1},kbpts(10))-- good

--- This bunch of code checks the MU1.2 conditions for an ideal I
--(1) Is the ideal perfect of grade 2 ?
pdm   = pdim coker gens I; --computes projective dimension of I and should be 2
dpth  = depth(I,R); -- computes the depth of I in R and should be 2
pdm == dpth -- checks if I is perfect 
pdm == 2 -- checks perfect grade 2
-- if we get both true then yes I is perfect grade 2


--(2) Does the ideal satisfy G_l, where l=analyticSpread of I ?
l=analyticSpread I -- compute analytic spread
rsltn = res coker gens I;
M=rsltn.dd_2;
r=rank M;
gCon  = for i from 0 to l-1 list (i, i+1<=codim minors(r-i,M)) -- this checks G_l condition
--- if the list has only true values then yes I satisfies G_l.

--(3) Is the min# of gens of I_1(phi) less than or equal to the analyticSpread ?
rsltn = res coker gens I; -- computes resolution of R/I
I1    = minors(1,rsltn.dd_2);
muI1  = rank source mingens I1
muI1<=l --- if this is true then yes
-- If conditions (1)-(3) hold then we are in good shape for using MU1.2

-- TFAE
-- (a) Q has the expected form
compare I; -- we compare the defining equations of R_I with the equations obtained from the Jacobian dual
eqnsdual==com(eqnsrees) -- if true then the defining equations of R_I have the expected form
-- (b) R_I is CM and I_(n-l)(phi)=I_1(phi)^(n-l) n=muI l=analyticSpread
muI= rank source mingens I
I1 = minors(1,rsltn.dd_2);
J=minors(muI-l, rsltn.dd_2)
JJ=I1^(muI-l)
J==JJ
-- (c) r(I)< l=analyticSpread and I_(n-l)(phi)=I_1(phi)^(n-l)
minred=minimalReduction(I) 
rednumber=reductionNumber(I,minred)
l=analyticSpread I
rednumber <= l
--- (d) After elementary row operations on phi I_(n-l)(phi')=I_1(phi)^(n-l)
-- not sure how to check this computationally.
---------------------------------------------------------------------
-- This next bunch of code checks the MU1.3 conditions for an ideal I
l=analyticSpread I -- compute analytic spread
d=dim R
gCon  = for i from 0 to d-1 list (i, i+1<=codim fittingIdeal(i,coker gens I))--checks G_d conditions
muI=rank source mingens I; -- #of minimal gens for I
muI>=d
l==d -- if true then consequence of MU1.3 is true
J=minimalReduction I
reductionNumber(I,J)
-------------------------------------------------------------------


-- Example from Lan's paper "On Rees algebras of linearly presented ideals"
-- of an ideal that satisfies G_2 but not G_3. This example doesn't seem to work
R=ZZ/31991[x,y,z]
phi=matrix{{x,0,z},{y,x,0},{0,y,x},{0,0,y}}
I=minors(3,phi)
F0=minors(3,phi)
F1=minors(3-1,phi)
F2=minors(3-2,phi)
F3=minors(3-3,phi)
codim F0
codim F1
codim F2
codim F3
fittingIdeal(2,coker gens I)

