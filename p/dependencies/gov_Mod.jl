module gov

export tax

function tax(gross, lambda, tau)
	return lambda * (gross ^ (1 - tau))
end

end # module