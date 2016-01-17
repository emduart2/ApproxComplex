checkD=(M)-> (
primcomponents=primaryDecomposition M;
k=#primcomponents;
asprimes= apply(k,i-> radical primcomponents_i)
)

changeT=()-> (
use T;
test=substitute(test,T);
)

changeR=()-> (
use R;
test=substitute(test,R);
)


compare= (eq1,eq2)->(
use S;
I = substitute(eq2,S);
eq1 == ideal I
)

uplusv = f->(
f2=contract(u,f);
G=f-u*f2;
g2=contract(v,G);
return (f2,g2);
)
