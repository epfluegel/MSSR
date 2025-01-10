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
# OUTPUT: Either false or the characteristic polynomial of S and a similarity transformation W bringing
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

computeBlockCompanionMatrix := proc(S, p)
# INPUT: S - a matrix
#        p - a prime number
# OUTPUT: The characteristic polynomial of S and a similarity transformation W bringing S into block-companion form
	local C, i, j, k, t, v, w, W, WW, R, r;
	
	t := LinearAlgebra[Dimension](S)[1];
	v := LinearAlgebra[Modular][Create](p, 1, t, random, integer);	
	userinfo(5, computeBlockCompanionMatrix, "Vector is", convert(v, 'list'));
	w[1] := LinearAlgebra[Modular][Create](p, 1, t, random, integer);
    w[1] := Matrix(1, 8, {(1, 1) = 0, (1, 2) = 1, (1, 3) = 1, (1, 4) = 0, (1, 5) = 0, (1, 6) = 1, (1, 7) = 1, (1, 8) = 0});
	W[1] := Matrix(convert(w[1], 'list'));
	#userinfo(5, computeBlockCompanionMatrix, W[1]);
	r := 1;
	j := 1;
	userinfo(5, computeBlockCompanionMatrix, "Starting new block: ",j);
	for i from 2 to t do
		print(i,j,W[j]);
		w[i] := LinearAlgebra[Modular][Multiply](p, w[i-1], S);
        userinfo(5, computeBlockCompanionMatrix, "Iterated vector is", convert(w[i], 'list'), "augmented rank is", LinearAlgebra[Modular][Rank](p, ArrayTools[Concatenate](1, W[j], w[i])));
		if LinearAlgebra[Modular][Rank](p, ArrayTools[Concatenate](1, W[j], w[i])) > r then
			userinfo(5, computeBlockCompanionMatrix, "Increasing size", r);
			W[j] := ArrayTools[Concatenate](1, W[j], w[i]);
			r := r+1
		else
			if j=1 then WW := W[1] else WW := ArrayTools[Concatenate](1, WW, W[j]) fi;
            userinfo(5, computeBlockCompanionMatrix, "Finish block of dimension", r);
            R := LinearAlgebra[Modular][RowEchelonTransform](p, LinearAlgebra[Modular][Copy](p, WW), true, true, true, false);
			userinfo(5, computeBlockCompanionMatrix, "Row-echelon indices: ", R[2]);
            k := 1; while R[2][k]=k do k := k+1 od; 
            userinfo(5, computeBlockCompanionMatrix, "Pick index: ", k);
			j := j+1; 
			
			w[i] := LinearAlgebra[Transpose](LinearAlgebra[UnitVector](k, t));
            #userinfo(5, computeBlockCompanionMatrix, w[i]);
            #LinearAlgebra[Modular][Create](p, 1, t, random, integer);
			W[j] := Matrix(convert(w[i], 'list'));
			userinfo(5, computeBlockCompanionMatrix, "Starting new block: # = ",j);
			userinfo(5, computeBlockCompanionMatrix, "Vector is", convert(w[i], 'list'));
			r := 1;
		fi;
	od;
    
    # Concatenate last block 
	userinfo(5, computeBlockCompanionMatrix, "WW is", WW);
    WW := ArrayTools[Concatenate](1, WW, W[j]);
	userinfo(5, computeBlockCompanionMatrix, WW);
	#W := Matrix([seq(convert(w[i], 'list'), i=1..t)]);
	#if LinearAlgebra[Modular][Determinant](p, W) = 0 then
	#	R := LinearAlgebra[Modular][RowEchelonTransform](p, W, true, true, true, false);
	#	userinfo(5, computeBlockCompanionMatrix, R);
	#	return(R)
	#fi;
	C := LinearAlgebra[Modular][Multiply](p, WW, LinearAlgebra[Modular][Multiply](p, S, LinearAlgebra[Modular][Inverse](p, WW)));
	#[W, LinearAlgebra[Row](C, t), C];
	#[seq(W[i], i=1..j)]
    C

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
