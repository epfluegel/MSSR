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

FrobeniusNormalForm := proc(A, p)
    local F, Q;

    # Check if inputs are valid
    if not type(A, Matrix) then
        error "The input A must be a valid matrix.";
    end if;
    if not type(p, prime) then
        error "The input p must be a prime number.";
    end if;

    # Reduce the matrix A modulo p
    A := LinearAlgebra[Mod](A, p);

    # Compute the Frobenius normal form and the similarity transformation
    F, Q := LinearAlgebra[FrobeniusForm](A, 'transform' = true);

    # Return both the Frobenius form and the transformation matrix
    return F, Q;
end proc;


computeShiftedHessenbergForm := proc(S, p)
# INPUT: S - a matrix
#        p - a prime number
# OUTPUT: The characteristic polynomial of S and a similarity transformation W bringing S into block-triagonal companion form
	local C, i, j, k, t, v, vv, w, W, WW, R, r, rr, d, dlist, h;
	
	t := LinearAlgebra[Dimension](S)[1];
	v := LinearAlgebra[Modular][Create](p, 1, t, random, integer);	
    v := Matrix(1, 8, {(1, 1) = 0, (1, 2) = 1, (1, 3) = 1, (1, 4) = 0, (1, 5) = 0, (1, 6) = 1, (1, 7) = 1, (1, 8) = 0});
    WW := LinearAlgebra[Modular][Copy](p, v);
    W := LinearAlgebra[Modular][Copy](p, v);
	#v := Matrix(convert(v, 'list'));
	userinfo(5, computeShiftedHessenbergForm, "Starting with row vector", convert(v, 'list'));

	r := 1;
	#j := 1;
    d[0]:= 0;
    h := 1;
    #dlist := [];
	userinfo(5, computeShiftedHessenbergForm, "Starting new block: ", h);
	for i from 2 to t do
		#print(i,j,WW);
		vv := LinearAlgebra[Modular][Multiply](p, v, S);
        userinfo(5, computeShiftedHessenbergForm, "Iterated vector is", convert(vv, 'list'));
        
        WW:= ArrayTools[Concatenate](1, W, vv);
        rr := LinearAlgebra[Modular][Rank](p, WW);
        userinfo(5, computeShiftedHessenbergForm, "New rank is", rr);
        
		if rr > r then
			userinfo(5, computeShiftedHessenbergForm, "Increasing size", rr);
            v := vv;
            W := LinearAlgebra[Modular][Copy](p, WW);
			r := rr;
            #d := d+1;
		else
            userinfo(5, computeShiftedHessenbergForm, "Finish block of dimension", i-1-d[h-1]);
            
            R := LinearAlgebra[Modular][RowEchelonTransform](p, LinearAlgebra[Modular][Copy](p, W), true, true, true, false);
			userinfo(5, computeShiftedHessenbergForm, "Row-echelon indices: ", R[2],nops(R[2]));
            k := 1; while R[2][k]=k do k := k+1 od; 
            #k := 7;
            userinfo(5, computeShiftedHessenbergForm, "Pick index: ", k); 			
			v := LinearAlgebra[Transpose](LinearAlgebra[UnitVector](k, t));
			userinfo(5, computeShiftedHessenbergForm, "Vector is", convert(v, 'list'));            
            W := ArrayTools[Concatenate](1, W, v);
            
            #dlist := [op(dlist), d];
            d[h] := i-1;
            h := h+1;
			r := rr + 1;
            #d := 1;
            userinfo(5, computeShiftedHessenbergForm, "Starting new block: ", h);
            userinfo(5, computeShiftedHessenbergForm, "New rank is", LinearAlgebra[Modular][Rank](p, W));
  

		fi;
	od;
    
    d[h] := t;
    print(d[h], d[h-1]);
    userinfo(5, computeShiftedHessenbergForm, "Finish block of dimension", t-d[h-1]);
    #dlist := [op(dlist), d];    
	C := LinearAlgebra[Modular][Multiply](p, W, LinearAlgebra[Modular][Multiply](p, S, LinearAlgebra[Modular][Inverse](p, W)));
    [d, S, C]

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
