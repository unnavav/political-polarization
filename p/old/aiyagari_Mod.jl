module aiyagari

export prices, uc, getW

function prices(K, L, alpha, delta)
	w = (1-alpha)*(K^alpha)*(L^(-alpha))
	r = (alpha)*((K/L)^(alpha-1)) - delta
	return r, w
end

function uc(apr, params)
	y = params[1]
	beta = params[2]
	sigma = params[3]
	
	cons = y - apr
	if cons <= 0
		u = 2000
	elseif sigma == 1
		u = log(cons)
	else
		u = (cons)^(1-sigma)/(1-sigma)
	end
	
	return u
end

function getW(r, alpha, delta)
	w = r .+ delta
	w = w ./ alpha
	w = w .^ (alpha / (alpha - 1))
	w = w .* (1 - alpha)
	return w
end

function util(apr, params, vc)
	y = params[1]
	beta = params[2]
	sigma = params[3]

	if sigma == 1
		utils = log(y - apr) + beta * vc
	else
		utils = ((y - apr)^(1 - sigma)) / (1 - sigma) + beta * vc
	end

	return utils

end

end