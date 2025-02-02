MSSR_eval := proc()

	local N, M, step, samples, ssample, size, sizes, r, k, success, result_p, result_d, tavg, lmax, lmin, lavg, p;
	
	N := 10: M := 500: step := 50: samples := 100: sizes := [10,20,30,40,50,100,200,300,400,500,1000]:
	for size in sizes do
		r := rand(10^(size-1)..10^(size)-1):
		for k from 1 to samples do
			ssample[size,k] := r():
			#print(ssample[size,k]):
		od:
		#print(size):
	od:
	result_p := []:
	for size in sizes do
		tavg := 0:
		lmin := infinity: lmax := - infinity: lavg := 0:
		for k from 1 to samples do
			success := reduceShareSize_p(ssample[size,k]):
			if success = false then break fi:
			lmin := min(lmin, success[1]):
			lmax := max(lmax, success[1]):
			lavg := lavg + success[1]:
			tavg := tavg + success[4]:
		od:
		lavg := lavg/samples:
		result_p := [op(result_p), [lmin, lmax, evalf(lavg), evalf(tavg)]]:
		
	od:


	result_d := [];
	for size in sizes do
		for p in [2,3,5,7,11,13] do
			tavg := 0;
			for k to samples do
					success := reduceShareSize_d(ssample[size, k], p);
					if success = false then
							break;
					end if;
					tavg := tavg + success[1];
			end do;
			result_d := [op(result_d), [p, evalf(tavg)]];
		od;
	end do;
	
	return(result_p, result_d);

end:

#	result_p;
#		[[1, 5, 2.290000000, 0.547], [1, 4, 2.520000000, 0.750], 
#
#			[1, 5, 2.300000000, 0.797], [1, 6, 2.460000000, 0.984], 
#
#			[1, 5, 2.280000000, 1.234], [1, 5, 2.700000000, 2.375], 
#
#			[1, 9, 2.320000000, 5.234], [1, 5, 2.420000000, 8.704], 
#
#			[1, 5, 2.320000000, 12.656], [1, 6, 2.320000000, 16.625], 
#
#			[1, 5, 2.060000000, 43.640]]
#	result_d;
#	[0.344, 0.812, 0.750, 1.266, 1.484, 3.516, 7.750, 14.000, 19.781, 
#
#		26.922, 67.641]
#