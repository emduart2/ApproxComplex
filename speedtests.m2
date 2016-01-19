-- In this session I would like to test how far up in degree can we go with 
-- the theorem from Tensor product surfaces and linear syzygies. I also
-- want to know what is the advantage in computation time of having a linear syzygy.
-- write several examples using just the code zcomplex vs using linearsyz
-- write examples of how much time does a computation take if there are no linear
-- syzygies
-- how far up in degree can we go if there is a linear syzygy vs when there is no linear 
-- syzygy

use R
I=blin(5,2)
bgb coker gens I
I=substitute(I,T)
zcomplex(I,3,2)
