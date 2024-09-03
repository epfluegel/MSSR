convertSecret := proc(s, p)
# INPUT: s - an integer
#        p - a prime number < s
# OUTPUT: A t by t matrix which is s converted to a matrix, using its base p digits
	local t;

	t := ceil(sqrt(1+floor(log(s)/log(p))));
	Matrix(t, t, convert(s, base, p));

end:

computeCompanionMatrix := proc(S, p)
# INPUT: S - a matrix
#        p - a prime number
# OUTPUT: Either fail or the characteristic polynomial of S and a similarity transformation W bringing
# 	  S into companion form
	local C, i, t, v, w, W;
	
	t := LinearAlgebra[Dimension](S)[1];
	v := LinearAlgebra[Modular][Create](p, 1, t, random, integer);	
	w[1] := v;
	for i from 2 to t do
		w[i] := LinearAlgebra[Modular][Multiply](p, w[i-1], S)
	od;
	W := Matrix([seq(convert(w[i], 'list'), i=1..t)]);
	if LinearAlgebra[Modular][Determinant](p, W) = 0 then
		return false;
	fi;
	C := LinearAlgebra[Modular][Multiply](p, W, LinearAlgebra[Modular][Multiply](p, S, LinearAlgebra[Modular][Inverse](p, W)));
	[W, LinearAlgebra[Row](C, t), C];

end:

reduceShareSize := proc(s)
# INPUT: S - a matrix
# OUTPUT: P - a public matrix
#         Q - a private row vector
#         info - additional test info
	local l, p, t0, success;

	t0 := time();
	success := false;
	p := 1;
	l := 0;
	while (success=false) and (p <= evalf(sqrt(s))) do
		p := nextprime(p);
		l := l+1;
		success := computeCompanionMatrix(convertSecret(s, p), p);
	od;
	[l, p, ceil(sqrt(1+floor(log(s)/log(p)))), time()-t0, success]
end:

createShares := proc(s, p)
# INPUT: s - an integer
#        p - a prime number < s
# OUTPUT: A t by t matrix which is s converted to a matrix, using its base p digits
	
	convert(s, base, p);

end:

reconstructSecret := proc(s, p)
# INPUT: s - an integer
#        p - a prime number < s
# OUTPUT: A t by t matrix which is s converted to a matrix, using its base p digits
	
	convert(s, base, p);

end:
