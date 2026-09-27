using LinearAlgebra
using SparseArray


A = [
     3.0   6.0  -1.0
    -6.0  -6.0   1.0
     2.0   1.0  -1.0
]

println("A =")
display(A)

# Question 1a)

# Gram-Schmidt
a1 = A[:, 1]
a2 = A[:, 2]
a3 = A[:, 3]

u1 = a1
e1 = u1 / norm(u1)

u2 = a2 - dot(a2, e1) * e1
e2 = u2 / norm(u2)

u3 = a3 - dot(a3, e1) * e1 - dot(a3, e2) * e2
e3 = u3 / norm(u3)

Q = [e1 e2 e3]
R = Q' * A

# Clean tiny roundoff values
R[abs.(R) .< 1e-12] .= 0.0

println("Q =")
display(Q)

println("R =")
display(R)

println("Q'Q =")
display(Q' * Q)

println("QR =")
display(Q * R)

println("Q is orthogonal: ", Q' * Q ≈ I)
println("QR = A: ", Q * R ≈ A)


function givens_qr(A; count_ops=false)

    # Convert A to a Float64 matrix so all arithmetic is floating-point
    A = Matrix{Float64}(A)

    # m = number of rows, n = number of columns
    m, n = size(A)

    # This implementation assumes A is square or tall
    if m < n
        error("wide matrices not supported!")
    end

    # Q starts as the m x m identity matrix.
    # As we apply Givens rotations to A, we accumulate their transposes in Q.
    Q = Matrix{Float64}(I, m, m)

    # Create an operation counter only if requested
    if count_ops
        ops = OpCount()
    end

    # Work column by column
    for p = 1:n

        # Move upward from the bottom of the matrix,
        # eliminating entries below A[p,p]
        for q = m:-1:p+1

            # r1 is the current pivot entry
            # r2 is the entry below the pivot that we want to eliminate
            r1 = A[p,p]
            r2 = A[q,p]

            # Compute
            #
            #     z = sqrt(r1^2 + r2^2)
            #
            # This is the norm of the 2-vector [r1, r2].
            z = sqrt(r1*r1 + r2*r2)

            # Count the work used to compute z
            if count_ops
                ops.mul += 2      # r1*r1 and r2*r2
                ops.add += 1      # r1^2 + r2^2
                ops.sqrt += 1     # sqrt(...)
            end

            # Choose c and s so the Givens rotation sends
            #
            #     [r1]
            #     [r2]
            #
            # to
            #
            #     [z]
            #     [0]
            #
            if z > 0
                c = r1 / z
                s = r2 / z

                if count_ops
                    ops.div += 2
                end
            else
                # If r1 = r2 = 0, no rotation is needed
                c = 1.0
                s = 0.0
            end

            # Construct the Givens rotation matrix.
            #
            # RotM is the identity except for rows/columns p and q:
            #
            #     [ c   s ]
            #     [-s   c ]
            #
            M = RotM(m, p, q, c, s)

            # Apply the rotation to A.
            #
            # This zeros A[q,p].
            A = M * A

            # Accumulate Q.
            #
            # Since R = M_k ... M_2 M_1 A,
            # we have
            #
            #     Q = M_1' M_2' ... M_k'
            #
            # so that A_original = Q*R.
            Q = Q * M'

            # Count the arithmetic required to apply the Givens rotation
            # to the active columns of A.
            #
            # For each column, the two updated entries are
            #
            #     a_new = c*a + s*b
            #     b_new = c*b - s*a
            #
            # giving:
            #     4 multiplications
            #     1 addition
            #     1 subtraction
            #
            # We do NOT count the work used to construct Q.
            if count_ops
                active_cols = n - p + 1

                ops.mul += 4 * active_cols
                ops.add += active_cols
                ops.sub += active_cols
            end
        end
    end

    # Return Q and the upper-triangular matrix R.
    # If operation counting is enabled, also return ops.
    if count_ops
        return Q, A, ops
    else
        return Q, A
    end
end