def main():
	
	global leakPattern
	global Threshold
	global mode
	global sampling
    global n, q, sigma

	R = Integers(q)
	Q = PolynomialRing(R, 'y')
	S = Q.quo(Q.gen() ** (n) + 1, 'x')

	w = find_omega()
	Omega = R(w[0])

    # Load NTT Matrix for given instance
	nttMatrix = load_ntt_matrix(challID)

    # Create basis for CVP solver
	[mat, B] = create_basis()

    sk = S(s)

	uTotal = []
	aTotal = []

	count = 0

	for measidx in range(numSample):
		flag = True

		skNTT = nttMatrix * vector(sk)

		e = b - a * sk
		
		# Create the matrix form of public key a
		aMatrix = create_a(a, n)

		
		# Apply NTT transform to error e
		eNTT = nttMatrix * vector(e)

		# Compute u = a * sk + e (not going to use it)
		u = a * sk + e

		# uPrime is same as u but in NTT form
		uPrime = aMatrix * vector(sk) + vector(e)
		uTarget = uPrime[:]

		ePoints = leak(eNTT, leakRate)

		ePoly1 = Q.lagrange_polynomial(ePoints)
			eCoeff1 = coeff_fix(ePoly1.list(), leakRate / 2)

			ePoly2 = Q.lagrange_polynomial(ePoints[1])
			eCoeff2 = coeff_fix(ePoly2.list(), leakRate / 2)

		# Start_Time = time.time()

		if leakPattern == '1-7mod16' or leakPattern == '1-15mod16' or leakPattern == '1mod16':
			numEqns = n / 8
		elif leakPattern == '1mod8':
			numEqns = n / 4

		p = multiprocessing.Process(target=worker2, args=(
			[status, aMat, uMat, countStatus]sysIdx, status, aMat, uMat, countStatus, mat, B, eCoeff1, eCoeff2, eList, uTarget, aMatrix))
			jobs.append(p)
			p.start()

		
		status = status.values()
		aMat = aMat.values()
		uMat = uMat.values()
		countStatus = countStatus.values()

		idx = 0
		for st in status:
			if st == False:
				flag = False
				break
			countTemp = countStatus[idx]
			if countTemp > 0:
				count += countTemp
				aTotal.extend(aMat[idx])
				uTotal.extend(uMat[idx])
			idx += 1

			# If we have enough rows
			if count == n:
				break;

		measidx += 1
		if flag == False:
			break;

		if count == n:
			break;

	if flag == True:
		aTotal = matrix(aTotal)
		uTotal = vector(uTotal)
		
		try:
			KEY = aTotal.solve_right(uTotal)
		except:
			print("couldn't solve")
			return 0
		found_Key
		
return