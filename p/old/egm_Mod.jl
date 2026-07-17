module egm

export numdev, aubdev, solveD

function numdev(EV, agrid)
	nl, na, np = size(EV)
	DEV = zeros(size(EV))

	for ip in 1:np
		for ia in 1:na
			evia = EV[:, ia, ip]
			a = agrid[ia]

			if ia > 1 && ia < na
				evia_0 = EV[:, ia-1, ip]
				evia_1 = EV[:, ia+1, ip]
				a1 = agrid[ia+1]
				a0 = agrid[ia-1]

				dr = (evia_1 - evia) / (a1 - a)
				dl = (evia - evia_0) / (a - a0)

				DEV[:, ia, ip] = (dr + dl) / 2

			elseif ia == 1
				evia_1 = EV[:, ia+1, ip]
				a1 = agrid[ia+1]
				dr = (evia_1 - evia) / (a1 - a)

				DEV[:, ia, ip] = dr

			else
				evia_0 = EV[:, ia-1, ip]
				a0 = agrid[ia-1]
				dl = (evia - evia_0) / (a - a0)

				DEV[:, ia, ip] = dl
			end
		end
	end

	return DEV
end

function aubdev(ev, agrid)
	m, anum = size(ev)
	Dev = zeros(m, anum)

	for ia in 1:anum
		if ia == anum
			Dev[:, ia] = Devr[:, 1]
		elseif ia == 1
			Dev[:, ia] = (ev[:, ia+1] - ev[:, ia]) / (agrid[ia+1] - agrid[ia])
			Devl[:, 1] = Dev[:, ia]
		else
			Devr[:, 1] = (ev[:, ia+1] - ev[:, ia]) / (agrid[ia+1] - agrid[ia])
			Dev[:, ia] = (Devl + Devr) / 2.0
			Devl = Devr
		end
	end

	return Dev
end

function solveD(EV0, ik, kchgrid)
	nk = length(kchgrid)
	k = kchgrid[ik]

	if ik == nk
		d = (EV0[ik] - EV0[ik-1]) / (k - kchgrid[ik-1])
	else
		d = (EV0[ik+1] - EV0[ik]) / (kchgrid[ik+1] - k)
	end

	return d
end

end # module