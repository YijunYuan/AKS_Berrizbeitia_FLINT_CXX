AKS_Berrizbeitia_FLINT_CXX
==========================
Well...This is the first version of my implementation of the Berrizbeitia improved version of AKS Algorithm.
I finished the part of n=4k+1, and the another part is under coding.
This program is completely based on Bill Hart's FLINT, which can be found at http://www.flintlib.org
I used to write it in pure C, but I finally find it hard to generate a Set in C. So I turn to C++. The set in STL is used.

Some details are here:
	1.I modified some polynomials in the algorithm to make the program faster. In my code, it is called poly2. It's not difficult to see that (1 + m*x)^n=1 + m*(a^[n / 2^s] mod n)*x^[n % 2^s]  (mod x^2^s-a, n) . This arrangement transfers the calculation of the mod of polys to the powermods of some "small" number. And it extends the max virtual value (MVV) of this implementation, for fmpz_mod_poly_set_coeff_fmpz only receive a ulong type for its second parameter.
	2.The exact MVV is hard to evaluate, but if is believed that the MVV is bigger than 2^31 bit. It's enough for using.
	3.Test is needed. I think there may be some BUGs are in this code. If someone find some thing wrong here, I'm glad to hear from you.

Yijun Yuan
2014 July 21th
