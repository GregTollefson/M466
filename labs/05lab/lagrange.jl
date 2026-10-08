ellf(k,t)=prod([x[1:k-1];x[k+1:length(x)]].-t)/
          prod([x[1:k-1];x[k+1:length(x)]].-x[k])

function ellp(k,t)
	p=1
	q=1
	for i=1:length(x)
		if i!=k
			p*=x[i]-t
			q*=x[i]-x[k]
		end
	end
	return p/q
end

pf(t)=sum((k->y[k]*ellf(k,t)).(1:length(x)))

function pp(t)
	s=0
	for k=1:length(x)
		s+=y[k]*ellp(k,t)
	end
	return s
end
