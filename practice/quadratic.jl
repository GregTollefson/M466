function quadratic(a, b, c)
	a == 0 && throw(DomainError(a, "a must be nonzero for a quadratic polynomial"))

	discriminant = b^2 - 4a*c
	discriminant > 0 || throw(DomainError(discriminant, "the discriminant must be positive"))

	return ((-b + sqrt(discriminant)) / (2a), (-b - sqrt(discriminant)) / (2a))
end