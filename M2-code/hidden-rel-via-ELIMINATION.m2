--------------------------------------------------------------------------------
-- MACAULAY2 CODE TO COMPUTE THE HIDDEN RELATION VIA ELIMINATION OF VARIABLES --
--------------------------------------------------------------------------------

-- Hidden relation in terms of t_0..t_2,d_0..d_2, and L=\Lambda=\prod \lambda_k 
-- The final polynomial is H = H_0+H_1*L+H_2*L^2

printWidth=120;
gbTrace=3;

kk = QQ
R = kk[a_0..a_5,t_0..t_2,d_0..d_2,L, MonomialOrder => Eliminate 6]

-- Formulas for the traces
T_0 = t_0 - (-a_0-a_5) ;
T_1 = t_1 - (a_0+a_4-a_5) ;
T_2 = t_2 - (-a_0+a_1+a_5) ;

-- Formulas for the determinants
D_0 = d_0 - (-a_2*a_3+a_0*a_5) ;
D_1 = d_1 - (-a_1*a_3+a_2*a_3+a_0*a_4-a_0*a_5) ;
D_2 = d_2 - (a_2*a_3-a_2*a_4-a_0*a_5+a_1*a_5) ;

-- Explicit formula for X
x = L*(-a_0^2*a_1^2 + 4*a_0^3*a_2 + 27*a_2^2*a_3^2 - 4*a_2*a_4^3 + 4*a_3*a_5^3 - (a_1^2 - 12*a_0*a_2)*a_4^2 - (a_0^2 + 12*a_1*a_3 - 2*a_0*a_4 + a_4^2)*a_5^2 - 2*(2*a_1^3 - 9*a_0*a_1*a_2)*a_3 + 2*(a_0*a_1^2 - 6*a_0^2*a_2 - 9*a_1*a_2*a_3)*a_4 + 2*(a_0^2*a_1 + a_1*a_4^2 + 3*(2*a_1^2 - 3*a_0*a_2)*a_3 - (2*a_0*a_1 - 9*a_2*a_3)*a_4)*a_5) - (a_2^2*a_3^2 - a_1*a_2*a_3*a_4 + a_0*a_2*a_4^2 + a_0^2*a_5^2 - (a_0*a_1*a_4 - (a_1^2 - 2*a_0*a_2)*a_3)*a_5);

-- Ideal of relations
J = ideal(T_0,T_1,T_2,D_0,D_1,D_2,x);

-- Main computation
time g = gb(J, SubringLimit => 1)
time G = gens g;

-- The hidden relation
h = G_(0,0)

-- h = H_2*L^2 + H_1*L + H_0
H_0 = sub(h,{L=>0});
H_1 = sub(diff(L,h),{L=>0});
H_2 = sub(diff(L,diff(L,h))/2,{L=>0});
sub(diff(L,diff(L,diff(L,h)))/6,{L=>0}) == 0

-- Work with a sub-family (eg. Hamiltonian vector fields)
subVals = {t_0 => 0, t_1 => 0, t_2 => 0}
hVals = sub(h, subVals)
factor hVals

--------------------------
-- Check irreducibility --
--------------------------

isPrime h  --true

--------------------
-- Export the H_i --
--------------------
printWidth=0;

toString factor H_0
toString factor H_1
toString factor H_2
