# Gram-Schmidt QR Factorization

The goal of the Gram-Schmidt process is to factor a matrix

$$
A
$$

into

$$
A = QR,
$$

where the columns of $Q$ are orthonormal and $R$ is upper triangular.

For

$$
A =
\begin{bmatrix}
3 & 6 & -1 \\
-6 & -6 & 1 \\
2 & 1 & -1
\end{bmatrix},
$$

write the columns of $A$ as

$$
A =
\begin{bmatrix}
a_1 & a_2 & a_3
\end{bmatrix}.
$$

The Gram-Schmidt process takes the original column vectors

$$
a_1,\;a_2,\;a_3
$$

and constructs an orthonormal set

$$
e_1,\;e_2,\;e_3.
$$

These vectors become the columns of $Q$.

---

# 1. Julia Function

```julia
function gramSchmidt_qr(A; count_ops=false)

    A = Matrix{Float64}(A)
    m, n = size(A)

    if n != 3
        error("This implementation expects exactly 3 columns.")
    end

    if count_ops
        ops = OpCount()
    end

    a1 = A[:,1]
    a2 = A[:,2]
    a3 = A[:,3]

    u1 = a1
    e1 = u1 / norm(u1)

    if count_ops
        ops.mul += m
        ops.add += m - 1
        ops.sqrt += 1
        ops.div += m
    end

    u2 = a2 - dot(a2,e1)*e1

    if count_ops
        ops.mul += 2*m
        ops.add += m - 1
        ops.sub += m
    end

    e2 = u2 / norm(u2)

    if count_ops
        ops.mul += m
        ops.add += m - 1
        ops.sqrt += 1
        ops.div += m
    end

    u3 = a3 - dot(a3,e1)*e1 - dot(a3,e2)*e2

    if count_ops
        ops.mul += 4*m
        ops.add += 2*(m - 1)
        ops.sub += 2*m
    end

    e3 = u3 / norm(u3)

    if count_ops
        ops.mul += m
        ops.add += m - 1
        ops.sqrt += 1
        ops.div += m
    end

    Q = [e1 e2 e3]

    R = Q' * A

    if count_ops
        # Optional literal implementation count for forming R
        ops.mul += 9*m
        ops.add += 9*(m - 1)
    end

    if count_ops
        return Q, R, ops
    else
        return Q, R
    end

end
```

---

# 2. Extract the Columns of A

The lines

```julia
a1 = A[:,1]
a2 = A[:,2]
a3 = A[:,3]
```

extract the three columns of $A$.

Thus,

$$
a_1 =
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix},
\qquad
a_2 =
\begin{bmatrix}
6\\
-6\\
1
\end{bmatrix},
\qquad
a_3 =
\begin{bmatrix}
-1\\
1\\
-1
\end{bmatrix}.
$$

The purpose of Gram-Schmidt is to construct orthonormal vectors from these original columns.

---

# 3. Construct the First Orthonormal Vector

The first step is

```julia
u1 = a1
```

so

$$
u_1=a_1.
$$

For this example,

$$
u_1=
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix}.
$$

The vector is normalized using

```julia
e1 = u1 / norm(u1)
```

Mathematically,

$$
e_1=\frac{u_1}{\|u_1\|}.
$$

The norm is

$$
\|u_1\|
=
\sqrt{3^2+(-6)^2+2^2}
=
\sqrt{49}
=
7.
$$

Therefore,

$$
e_1
=
\frac{1}{7}
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix}.
$$

Since $e_1$ has unit length,

$$
e_1^Te_1=1.
$$

---

# 4. Operation Count for Normalization

For a vector of length $m$,

$$
\|u\|
=
\sqrt{u_1^2+\cdots+u_m^2}
$$

requires

$$
m
$$

multiplications,

$$
m-1
$$

additions, and one square root.

Dividing every component by $\|u\|$ requires

$$
m
$$

divisions.

For the present problem,

$$
m=3.
$$

Thus each normalization requires

$$
3\text{ multiplications},
$$

$$
2\text{ additions},
$$

$$
3\text{ divisions},
$$

and

$$
1\text{ square root}.
$$

---

# 5. Remove the Component of a2 in the e1 Direction

The second vector begins with

```julia
u2 = a2 - dot(a2,e1)*e1
```

Mathematically,

$$
u_2
=
a_2-(a_2^Te_1)e_1.
$$

The scalar

$$
a_2^Te_1
$$

measures the component of $a_2$ in the direction of $e_1$.

The vector

$$
(a_2^Te_1)e_1
$$

is therefore the projection of $a_2$ onto $e_1$.

Thus

$$
u_2
=
a_2-\operatorname{proj}_{e_1}(a_2).
$$

This removes the component of $a_2$ parallel to $e_1$.

As a result,

$$
u_2^Te_1=0.
$$

Therefore $u_2$ is orthogonal to $e_1$.

---

# 6. Projection Formula

For a unit vector $e$, the projection of $a$ onto $e$ is

$$
\operatorname{proj}_e(a)
=
(a^Te)e.
$$

Because the vectors $e_i$ produced by Gram-Schmidt are normalized, no denominator is required.

For a general non-unit vector $u$, the projection would instead be

$$
\operatorname{proj}_u(a)
=
\frac{a^Tu}{u^Tu}u.
$$

Using normalized vectors simplifies the Gram-Schmidt formulas.

---

# 7. Normalize the Second Vector

Once the component parallel to $e_1$ is removed,

```julia
e2 = u2 / norm(u2)
```

gives

$$
e_2
=
\frac{u_2}{\|u_2\|}.
$$

Now

$$
e_2^Te_2=1
$$

and

$$
e_1^Te_2=0.
$$

Thus

$$
e_1
\quad\text{and}\quad
e_2
$$

are orthonormal.

---

# 8. Construct the Third Orthogonal Vector

The third vector must be orthogonal to both $e_1$ and $e_2$.

The code is

```julia
u3 = a3 - dot(a3,e1)*e1 - dot(a3,e2)*e2
```

which corresponds to

$$
u_3
=
a_3
-
(a_3^Te_1)e_1
-
(a_3^Te_2)e_2.
$$

This removes the component of $a_3$ in the $e_1$ direction and the component in the $e_2$ direction.

Equivalently,

$$
u_3
=
a_3
-
\operatorname{proj}_{e_1}(a_3)
-
\operatorname{proj}_{e_2}(a_3).
$$

The resulting vector satisfies

$$
u_3^Te_1=0
$$

and

$$
u_3^Te_2=0.
$$

---

# 9. Normalize the Third Vector

The final orthonormal vector is

```julia
e3 = u3 / norm(u3)
```

or

$$
e_3
=
\frac{u_3}{\|u_3\|}.
$$

At this point,

$$
e_i^Te_j
=
\begin{cases}
1, & i=j,\\
0, & i\ne j.
\end{cases}
$$

Thus the three vectors form an orthonormal basis.

---

# 10. Construct Q

The code

```julia
Q = [e1 e2 e3]
```

places the orthonormal vectors into the columns of $Q$:

$$
Q=
\begin{bmatrix}
| & | & |\\
e_1 & e_2 & e_3\\
| & | & |
\end{bmatrix}.
$$

Since the columns are orthonormal,

$$
Q^TQ=I.
$$

Therefore $Q$ is orthogonal.

For a square matrix,

$$
Q^{-1}=Q^T.
$$

---

# 11. Check Orthogonality

The orthogonality of $Q$ can be verified in Julia with

```julia
CheckOrth = Q' * Q
```

The result should be approximately

$$
I=
\begin{bmatrix}
1&0&0\\
0&1&0\\
0&0&1
\end{bmatrix}.
$$

Because of floating-point arithmetic, values that should be zero may instead appear as numbers such as

$$
10^{-16}.
$$

These can be cleaned for display with

```julia
CheckOrth[abs.(CheckOrth) .< 1e-12] .= 0.0
```

This check is diagnostic and is not part of the operation count for the Gram-Schmidt algorithm.

---

# 12. Construct R

Once $Q$ is known and

$$
A=QR,
$$

multiply both sides on the left by $Q^T$:

$$
Q^TA
=
Q^TQR.
$$

Since

$$
Q^TQ=I,
$$

we obtain

$$
R=Q^TA.
$$

The Julia implementation uses

```julia
R = Q' * A
```

The result should be upper triangular:

$$
R=
\begin{bmatrix}
* & * & *\\
0 & * & *\\
0 & 0 & *
\end{bmatrix}.
$$

---

# 13. How R Appears Naturally in Gram-Schmidt

Although the Julia implementation calculates

```julia
R = Q' * A
```

after constructing $Q$, the entries of $R$ are already generated naturally during the Gram-Schmidt process.

For example,

$$
r_{11}=\|a_1\|,
$$

and the projection coefficient

$$
r_{12}=e_1^Ta_2.
$$

Similarly,

$$
r_{22}=\|u_2\|,
$$

while

$$
r_{13}=e_1^Ta_3
$$

and

$$
r_{23}=e_2^Ta_3.
$$

Finally,

$$
r_{33}=\|u_3\|.
$$

Therefore,

$$
R=
\begin{bmatrix}
r_{11} & r_{12} & r_{13}\\
0 & r_{22} & r_{23}\\
0 & 0 & r_{33}
\end{bmatrix}.
$$

An efficient implementation can store these values as they are calculated instead of recomputing

$$
R=Q^TA.
$$

---

# 14. Verify A = QR

The factorization can be checked with

```julia
CheckA = Q * R
```

and compared with the original matrix.

Numerically,

$$
QR\approx A.
$$

A useful numerical check is

```julia
norm(A - Q*R)
```

which should be close to machine precision.

---

# 15. Gram-Schmidt Algorithm Summary

Given

$$
A=
\begin{bmatrix}
a_1&a_2&a_3
\end{bmatrix},
$$

the classical Gram-Schmidt process is:

1. Set

   $$
   u_1=a_1.
   $$

2. Normalize:

   $$
   e_1=\frac{u_1}{\|u_1\|}.
   $$

3. Remove the $e_1$ component from $a_2$:

   $$
   u_2=a_2-(a_2^Te_1)e_1.
   $$

4. Normalize:

   $$
   e_2=\frac{u_2}{\|u_2\|}.
   $$

5. Remove the $e_1$ and $e_2$ components from $a_3$:

   $$
   u_3
   =
   a_3
   -
   (a_3^Te_1)e_1
   -
   (a_3^Te_2)e_2.
   $$

6. Normalize:

   $$
   e_3=\frac{u_3}{\|u_3\|}.
   $$

7. Form

   $$
   Q=
   \begin{bmatrix}
   e_1&e_2&e_3
   \end{bmatrix}.
   $$

8. Construct the upper-triangular factor $R$.

Then

$$
\boxed{A=QR}.
$$

---

# 16. Operation Counting

For the three-column Gram-Schmidt calculation, excluding any separate multiplication used to reconstruct $R$, the operations are counted directly from the orthogonalization process.

### Normalizations

There are three normalizations.

Each normalization of a 3-vector requires

$$
3\text{ multiplications}
+
2\text{ additions}
+
3\text{ divisions}
+
1\text{ square root}.
$$

Across all three normalizations:

$$
9\text{ multiplications},
$$

$$
6\text{ additions},
$$

$$
9\text{ divisions},
$$

and

$$
3\text{ square roots}.
$$

### Constructing u2

The expression

$$
u_2=a_2-(a_2^Te_1)e_1
$$

requires:

- one dot product:
  
  $$
  3\text{ multiplications}+2\text{ additions},
  $$

- one scalar-vector multiplication:

  $$
  3\text{ multiplications},
  $$

- one vector subtraction:

  $$
  3\text{ subtractions}.
  $$

Thus,

$$
u_2:
\qquad
6\text{ multiplications}
+
2\text{ additions}
+
3\text{ subtractions}.
$$

### Constructing u3

The expression

$$
u_3
=
a_3-(a_3^Te_1)e_1-(a_3^Te_2)e_2
$$

contains two projections.

Therefore,

$$
u_3:
\qquad
12\text{ multiplications}
+
4\text{ additions}
+
6\text{ subtractions}.
$$

Combining these gives the algorithmic Gram-Schmidt count

$$
\boxed{
12\text{ additions}
}
$$

$$
\boxed{
9\text{ subtractions}
}
$$

$$
\boxed{
27\text{ multiplications}
}
$$

$$
\boxed{
9\text{ divisions}
}
$$

and

$$
\boxed{
3\text{ square roots}.
}
$$

If the literal Julia calculation

```julia
R = Q' * A
```

is also counted, then for a $3\times3$ matrix this adds

$$
27\text{ multiplications}
+
18\text{ additions}.
$$

The corresponding literal implementation count becomes

$$
\boxed{
30\text{ additions},
\quad
9\text{ subtractions},
\quad
54\text{ multiplications},
\quad
9\text{ divisions},
\quad
3\text{ square roots}.
}
$$

For comparison between numerical methods, the algorithmic count is generally the more meaningful one because the entries of $R$ can be retained as they are computed during Gram-Schmidt.

# Givens QR Factorization

The goal of the Givens QR algorithm is to factor a matrix

$$
A
$$

into

$$
A = QR,
$$

where $Q$ is orthogonal,

$$
Q^TQ = I,
$$

and $R$ is upper triangular.

For the matrix

$$
A =
\begin{bmatrix}
3 & 6 & -1 \\
-6 & -6 & 1 \\
2 & 1 & -1
\end{bmatrix},
$$

Givens rotations eliminate the entries below the diagonal one at a time.

For a $3\times3$ matrix, three rotations are required:

$$
A_{31}\rightarrow0,
\qquad
A_{21}\rightarrow0,
\qquad
A_{32}\rightarrow0.
$$

---

## Julia Rotation Matrix

The helper function constructs a Givens rotation matrix.

```julia
function RotM(n, p, q, c, s)

    M = sparse(Float64, I, n, n)

    M[p,p] = c
    M[p,q] = s
    M[q,p] = -s
    M[q,q] = c

    return M
end
```

The matrix begins as the identity and replaces the entries associated with rows and columns $p$ and $q$ by

$$
\begin{bmatrix}
c & s\\
-s & c
\end{bmatrix}.
$$

For example, if

$$
p=1,
\qquad
q=3,
$$

then for a $3\times3$ matrix,

$$
M=
\begin{bmatrix}
c & 0 & s\\
0 & 1 & 0\\
-s & 0 & c
\end{bmatrix}.
$$

Only coordinates $p$ and $q$ are affected.

---

# 1. Julia QR Function

```julia
function givens_qr(A; count_ops=false)

    A = Matrix{Float64}(A)
    m, n = size(A)

    if m < n
        error("wide matrices not supported!")
    end

    Q = Matrix{Float64}(I, m, m)

    if count_ops
        ops = OpCount()
    end

    for p = 1:n

        for q = m:-1:p+1

            r1 = A[p,p]
            r2 = A[q,p]

            z = sqrt(r1*r1 + r2*r2)

            if count_ops
                ops.mul += 2
                ops.add += 1
                ops.sqrt += 1
            end

            if z > 0

                c = r1 / z
                s = r2 / z

                if count_ops
                    ops.div += 2
                end

            else

                c = 1.0
                s = 0.0

            end

            M = RotM(m, p, q, c, s)

            A = M * A
            Q = Q * M'

            if count_ops
                active_cols = n - p + 1

                ops.mul += 4 * active_cols
                ops.add += active_cols
                ops.sub += active_cols
            end
        end
    end

    if count_ops
        return Q, A, ops
    else
        return Q, A
    end
end
```

---

# 2. Moving Through the Matrix

The outer loop is

```julia
for p = 1:n
```

where $p$ identifies the current column.

The inner loop is

```julia
for q = m:-1:p+1
```

which starts at the bottom of the current column and moves upward.

For a $3\times3$ matrix:

- When $p=1$,

  $$
  q=3,2.
  $$

  The algorithm eliminates

  $$
  A_{31}
  \quad\text{and}\quad
  A_{21}.
  $$

- When $p=2$,

  $$
  q=3.
  $$

  The algorithm eliminates

  $$
  A_{32}.
  $$

Thus the entries below the diagonal are removed one at a time.

---

# 3. Select the Two Entries to Rotate

For each rotation,

```julia
r1 = A[p,p]
r2 = A[q,p]
```

selects two elements from the same column.

They form the two-element vector

$$
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}.
$$

The goal is to choose a rotation such that

$$
\begin{bmatrix}
c & s\\
-s & c
\end{bmatrix}
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}
=
\begin{bmatrix}
z\\
0
\end{bmatrix}.
$$

Thus the lower element $r_2$ is eliminated.

---

# 4. Compute the Length of the Two-Vector

The code calculates

```julia
z = sqrt(r1*r1 + r2*r2)
```

which is

$$
z=\sqrt{r_1^2+r_2^2}.
$$

This is the Euclidean norm of

$$
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}.
$$

The operation count is:

$$
2
$$

multiplications,

$$
1
$$

addition, and

$$
1
$$

square root.

---

# 5. Calculate the Rotation Parameters

If

$$
z>0,
$$

the code calculates

```julia
c = r1 / z
s = r2 / z
```

so

$$
c=\frac{r_1}{z},
\qquad
s=\frac{r_2}{z}.
$$

This requires two divisions.

These values satisfy

$$
c^2+s^2=1.
$$

This is what makes the Givens matrix orthogonal.

---

# 6. Why the Rotation Produces a Zero

Consider

$$
M
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}
=
\begin{bmatrix}
c & s\\
-s & c
\end{bmatrix}
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}.
$$

The first component is

$$
cr_1+sr_2.
$$

Substituting

$$
c=\frac{r_1}{z},
\qquad
s=\frac{r_2}{z},
$$

gives

$$
cr_1+sr_2
=
\frac{r_1^2+r_2^2}{z}.
$$

Since

$$
z^2=r_1^2+r_2^2,
$$

we obtain

$$
cr_1+sr_2=z.
$$

The second component is

$$
-sr_1+cr_2.
$$

Substituting $c$ and $s$ gives

$$
-\frac{r_2r_1}{z}
+
\frac{r_1r_2}{z}
=
0.
$$

Therefore,

$$
\boxed{
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}
\longrightarrow
\begin{bmatrix}
z\\
0
\end{bmatrix}
}
$$

under the Givens rotation.

---

# 7. Special Case When Both Entries Are Zero

If

$$
z=0,
$$

then

$$
r_1=r_2=0.
$$

There is nothing to eliminate.

The code therefore uses

```julia
c = 1.0
s = 0.0
```

which gives

$$
\begin{bmatrix}
1 & 0\\
0 & 1
\end{bmatrix}.
$$

Thus the transformation acts as the identity.

---

# 8. Construct the Givens Matrix

The line

```julia
M = RotM(m, p, q, c, s)
```

constructs the full $m\times m$ matrix.

It is equal to the identity except for the $p,q$ coordinate plane.

For example,

$$
M=
\begin{bmatrix}
c & s & 0\\
-s & c & 0\\
0 & 0 & 1
\end{bmatrix}
$$

when $p=1$ and $q=2$.

A Givens rotation satisfies

$$
M^TM=I.
$$

Therefore,

$$
M^{-1}=M^T.
$$

---

# 9. Apply the Rotation to A

The key step is

```julia
A = M * A
```

which applies the row rotation to the matrix.

If the current rotation operates on rows $p$ and $q$, then for each active column $j$,

$$
A_{pj}^{\text{new}}
=
cA_{pj}+sA_{qj},
$$

and

$$
A_{qj}^{\text{new}}
=
-sA_{pj}+cA_{qj}.
$$

Each pair of updated entries requires:

$$
4
$$

multiplications,

$$
1
$$

addition, and

$$
1
$$

subtraction.

This structured update is what is counted in the operation-count code.

---

# 10. First Column Elimination

For the first column,

$$
p=1.
$$

The algorithm begins at the bottom:

$$
q=3.
$$

The first rotation eliminates

$$
A_{31}.
$$

After this rotation, the first column has the form

$$
\begin{bmatrix}
*\\
*\\
0
\end{bmatrix}.
$$

The next rotation uses

$$
q=2
$$

and eliminates

$$
A_{21}.
$$

The first column then becomes

$$
\begin{bmatrix}
*\\
0\\
0
\end{bmatrix}.
$$

---

# 11. Second Column Elimination

Next,

$$
p=2.
$$

There is only one entry below the diagonal:

$$
A_{32}.
$$

The algorithm constructs one more Givens rotation and eliminates it.

The resulting matrix has the form

$$
R=
\begin{bmatrix}
* & * & *\\
0 & * & *\\
0 & 0 & *
\end{bmatrix}.
$$

Thus $R$ is upper triangular.

---

# 12. Accumulating Q

The code also performs

```julia
Q = Q * M'
```

but the assignment specifies that the operations needed to construct $Q$ from the rotations are not counted.

Suppose the rotations are

$$
M_1,M_2,\ldots,M_k.
$$

The transformation of $A$ is

$$
R=M_k\cdots M_2M_1A.
$$

Since each $M_i$ is orthogonal,

$$
M_i^{-1}=M_i^T.
$$

Therefore,

$$
A
=
M_1^TM_2^T\cdots M_k^TR.
$$

Thus

$$
Q=
M_1^TM_2^T\cdots M_k^T
$$

and

$$
\boxed{A=QR}.
$$

---

# 13. Verify the QR Factorization

The factorization can be checked in Julia with

```julia
Q, R = givens_qr(A)

println("Q =")
display(Q)

println("R =")
display(R)

println("Factorization error:")
display(norm(A - Q*R))

println("Orthogonality error:")
display(norm(Q'*Q - I))
```

The errors

$$
\|A-QR\|
$$

and

$$
\|Q^TQ-I\|
$$

should both be close to machine precision.

---

# 14. Givens Algorithm Summary

For each column $p$:

1. Start at the bottom of the column.

2. Select

   $$
   r_1=A_{pp},
   \qquad
   r_2=A_{qp}.
   $$

3. Calculate

   $$
   z=\sqrt{r_1^2+r_2^2}.
   $$

4. Set

   $$
   c=\frac{r_1}{z},
   \qquad
   s=\frac{r_2}{z}.
   $$

5. Construct the Givens rotation.

6. Apply the rotation to rows $p$ and $q$.

7. The entry

   $$
   A_{qp}
   $$

   becomes zero.

8. Continue upward through the column.

9. Move to the next column.

After all rotations,

$$
A\rightarrow R.
$$

The accumulated transposed rotations form $Q$.

---

# 15. Operation Counting

For one Givens rotation, computing

$$
z=\sqrt{r_1^2+r_2^2}
$$

requires

$$
2\text{ multiplications}
+
1\text{ addition}
+
1\text{ square root}.
$$

Computing

$$
c=\frac{r_1}{z},
\qquad
s=\frac{r_2}{z}
$$

requires

$$
2\text{ divisions}.
$$

For every active column, the row update

$$
\begin{aligned}
a' &= ca+sb,\\
b' &= cb-sa
\end{aligned}
$$

requires

$$
4\text{ multiplications}
+
1\text{ addition}
+
1\text{ subtraction}.
$$

For the $3\times3$ example, the rotations act on:

- 3 active columns for the first rotation,
- 3 active columns for the second rotation,
- 2 active columns for the third rotation.

Thus there are

$$
3+3+2=8
$$

two-entry updates.

The final operation count is therefore

$$
\boxed{
11\text{ additions}
}
$$

$$
\boxed{
8\text{ subtractions}
}
$$

$$
\boxed{
38\text{ multiplications}
}
$$

$$
\boxed{
6\text{ divisions}
}
$$

and

$$
\boxed{
3\text{ square roots}.
}
$$

The cost of constructing $Q$ from the rotations is not included.

# Householder QR Factorization

The goal of the Householder QR algorithm is to factor a matrix

$$
A
$$

into

$$
A = QR,
$$

where $Q$ is orthogonal,

$$
Q^TQ = I,
$$

and $R$ is upper triangular.

For the matrix

$$
A =
\begin{bmatrix}
3 & 6 & -1 \\
-6 & -6 & 1 \\
2 & 1 & -1
\end{bmatrix},
$$

two Householder reflections are required.

The first reflection eliminates the entries below $A_{11}$, and the second reflection eliminates the remaining entry below $A_{22}$.

---

## Julia Function

```julia
function householder_qr(A; count_ops=false)

    A = Matrix{Float64}(A)
    m, n = size(A)

    if m < n
        error("wide matrices not supported!")
    end

    Q = Matrix{Float64}(I, m, m)

    if count_ops
        ops = OpCount()
    end

    for p = 1:min(n, m-1)

        x = A[p:m, p]
        z = norm(x)

        if count_ops
            k = length(x)
            ops.mul += k
            ops.add += k - 1
            ops.sqrt += 1
        end

        if z > 0

            α = -copysign(z, x[1])

            v = copy(x)
            v[1] = v[1] - α

            if count_ops
                ops.sub += 1
            end

            v = v / norm(v)

            if count_ops
                k = length(v)

                ops.mul += k
                ops.add += k - 1
                ops.sqrt += 1
                ops.div += k
            end

            H = Matrix{Float64}(I, m, m)
            H[p:m, p:m] -= 2.0 * (v * v')

            if count_ops
                k = length(v)

                ops.mul += 2 * k^2
                ops.sub += k^2
            end

            A = H * A

            if count_ops
                ops.mul += m^2 * n
                ops.add += m * n * (m - 1)
            end

            # Do not count work used to construct Q
            Q = Q * H
        end
    end

    if count_ops
        return Q, A, ops
    else
        return Q, A
    end
end
```

---

# 1. Initialization

The function begins with

```julia
A = Matrix{Float64}(A)
m, n = size(A)
```

The matrix is converted to floating-point form and its dimensions are stored.

For the example,

$$
m=n=3.
$$

The matrix $Q$ initially starts as the identity:

```julia
Q = Matrix{Float64}(I, m, m)
```

so

$$
Q=I.
$$

As the Householder reflectors are generated, they will be accumulated into $Q$.

---

# 2. Moving Through the Columns

The main loop is

```julia
for p = 1:min(n, m-1)
```

For a $3\times3$ matrix,

$$
p=1,2.
$$

Therefore, two Householder reflections are constructed.

For $p=1$, entries below $A_{11}$ are eliminated.

For $p=2$, the entry below $A_{22}$ is eliminated.

---

# 3. Select the Active Part of the Column

The line

```julia
x = A[p:m, p]
```

selects the part of column $p$ beginning at the diagonal.

For the first reflection,

$$
p=1,
$$

so

$$
x =
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix}.
$$

For the second reflection, only rows 2 and 3 are used, so $x$ has length 2.

If

$$
k=\operatorname{length}(x),
$$

then the Euclidean norm is

$$
\|x\|_2
=
\sqrt{x_1^2+x_2^2+\cdots+x_k^2}.
$$

This corresponds to

```julia
z = norm(x)
```

and requires

$$
k
$$

multiplications,

$$
k-1
$$

additions, and one square root.

For the first column,

$$
x =
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix},
$$

so

$$
z
=
\sqrt{3^2+(-6)^2+2^2}
=
\sqrt{49}
=
7.
$$

---

# 4. Choose the Target Value

The next line is

```julia
α = -copysign(z, x[1])
```

This chooses the sign of $\alpha$ so that the Householder calculation is numerically stable.

Conceptually,

$$
\alpha
=
-\operatorname{sign}(x_1)\|x\|.
$$

For the first reflection,

$$
x_1=3>0,
$$

so

$$
\alpha=-7.
$$

The desired transformation is therefore

$$
Hx =
\begin{bmatrix}
-7\\
0\\
0
\end{bmatrix}.
$$

The purpose of the reflector is to preserve the length of $x$ while rotating or reflecting it onto the first coordinate axis.

---

# 5. Construct the Householder Vector

The code starts with a copy of $x$:

```julia
v = copy(x)
```

and then changes its first element:

```julia
v[1] = v[1] - α
```

Mathematically,

$$
v=x-\alpha e_1,
$$

where

$$
e_1=
\begin{bmatrix}
1\\
0\\
\vdots\\
0
\end{bmatrix}.
$$

For the first reflection,

$$
x=
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix},
\qquad
\alpha=-7,
$$

so

$$
v=
\begin{bmatrix}
3-(-7)\\
-6\\
2
\end{bmatrix}
=
\begin{bmatrix}
10\\
-6\\
2
\end{bmatrix}.
$$

The line

```julia
v[1] = v[1] - α
```

requires one subtraction.

---

# 6. Normalize the Householder Vector

The next line is

```julia
v = v / norm(v)
```

which produces a unit vector:

$$
v
\leftarrow
\frac{v}{\|v\|}.
$$

If $v$ contains $k$ elements, calculating its norm requires

$$
k
$$

multiplications,

$$
k-1
$$

additions, and one square root.

Dividing each element of $v$ by the scalar $\|v\|$ requires

$$
k
$$

divisions.

After this step,

$$
v^Tv=1.
$$

---

# 7. Construct the Householder Reflector

The Householder reflector is

$$
H=I-2vv^T.
$$

The code first constructs the identity matrix:

```julia
H = Matrix{Float64}(I, m, m)
```

and then replaces the active portion with

```julia
H[p:m, p:m] -= 2.0 * (v * v')
```

The expression

```julia
v * v'
```

is an **outer product**, not a dot product.

If

$$
v=
\begin{bmatrix}
v_1\\
v_2\\
\vdots\\
v_k
\end{bmatrix},
$$

then

$$
vv^T=
\begin{bmatrix}
v_1v_1 & v_1v_2 & \cdots & v_1v_k\\
v_2v_1 & v_2v_2 & \cdots & v_2v_k\\
\vdots & \vdots & \ddots & \vdots\\
v_kv_1 & v_kv_2 & \cdots & v_kv_k
\end{bmatrix}.
$$

Thus $vv^T$ is a $k\times k$ matrix.

The important properties of a Householder matrix are

$$
H^T=H
$$

and

$$
H^TH=I.
$$

Therefore $H$ is both symmetric and orthogonal.

It also satisfies

$$
H^{-1}=H.
$$

Applying the same Householder reflection twice returns the original vector.

---

# 8. Apply the Reflector to A

The key elimination step is

```julia
A = H * A
```

Mathematically,

$$
A_{\text{new}}=HA.
$$

For the first reflector,

$$
A_1=H_1A.
$$

The first column becomes

$$
\begin{bmatrix}
3\\
-6\\
2
\end{bmatrix}
\longrightarrow
\begin{bmatrix}
-7\\
0\\
0
\end{bmatrix}.
$$

Thus, after the first reflector, the matrix has the form

$$
A_1=
\begin{bmatrix}
* & * & *\\
0 & * & *\\
0 & * & *
\end{bmatrix}.
$$

The first column is now finished.

---

# 9. Second Householder Reflection

On the next iteration,

$$
p=2.
$$

The code

```julia
x = A[p:m, p]
```

now extracts only rows 2 through 3 of column 2.

Therefore,

$$
x=
\begin{bmatrix}
A_{22}\\
A_{32}
\end{bmatrix}.
$$

A second Householder reflector $H_2$ is constructed so that

$$
H_2x=
\begin{bmatrix}
*\\
0
\end{bmatrix}.
$$

After applying this reflector,

$$
R=H_2H_1A
$$

has the form

$$
R=
\begin{bmatrix}
* & * & *\\
0 & * & *\\
0 & 0 & *
\end{bmatrix}.
$$

Thus $R$ is upper triangular.

---

# 10. Accumulating Q

The code also performs

```julia
Q = Q * H
```

but the assignment specifies that the work required to construct $Q$ from the reflectors should not be included in the operation count.

Starting with

$$
Q=I,
$$

after the first reflector,

$$
Q=H_1.
$$

After the second reflector,

$$
Q=H_1H_2.
$$

Since

$$
R=H_2H_1A,
$$

and

$$
H_i^{-1}=H_i,
$$

we can solve for $A$:

$$
A=H_1H_2R.
$$

Therefore,

$$
Q=H_1H_2
$$

and

$$
\boxed{A=QR}.
$$

---

# 11. Verify the QR Factorization

The factorization can be checked in Julia with

```julia
Q, R = householder_qr(A)

println("Q =")
display(Q)

println("R =")
display(R)

println("Factorization error:")
display(norm(A - Q*R))

println("Orthogonality error:")
display(norm(Q'*Q - I))
```

The quantities

$$
\|A-QR\|
$$

and

$$
\|Q^TQ-I\|
$$

should both be close to machine precision.

---

# 12. Householder Algorithm Summary

For each column $p$:

1. Select the portion of the column below and including the diagonal:

   $$
   x=A_{p:m,p}.
   $$

2. Calculate

   $$
   z=\|x\|.
   $$

3. Choose

   $$
   \alpha
   =
   -\operatorname{sign}(x_1)z.
   $$

4. Form

   $$
   v=x-\alpha e_1.
   $$

5. Normalize $v$:

   $$
   v\leftarrow\frac{v}{\|v\|}.
   $$

6. Construct the reflector:

   $$
   H=I-2vv^T.
   $$

7. Apply it to the matrix:

   $$
   A\leftarrow HA.
   $$

8. Repeat for the next column.

After the final reflection,

$$
A\rightarrow R.
$$

The product of the reflectors gives the orthogonal factor $Q$.

---

# 13. Operation-Counting Note

There are two possible ways to count the Householder operations.

The current Julia implementation explicitly constructs

$$
H=I-2vv^T
$$

and then evaluates

```julia
A = H * A
```

as a matrix multiplication.

A more efficient numerical implementation generally does **not** explicitly construct $H$. Instead, it uses

$$
HA
=
(I-2vv^T)A
$$

and therefore

$$
HA
=
A-2v(v^TA).
$$

The latter uses the special structure of the Householder reflector and requires fewer operations.

For the assignment, the counting convention should therefore be stated explicitly. The current instrumentation can count the literal implementation first, while the structured Householder calculation can be considered separately.

# Testing Procedure (fror all three factorizations)

We begin with the known solution

$$
x_{\text{true}} =
\begin{bmatrix}
1 \\
2 \\
3
\end{bmatrix}.
$$

Using

$$
A =
\begin{bmatrix}
3 & 6 & -1 \\
-6 & -6 & 1 \\
2 & 1 & -1
\end{bmatrix},
$$

we construct the right-hand side

$$
b = Ax_{\text{true}}.
$$

Each QR method—Gram-Schmidt, Givens rotations, and Householder reflections—computes a factorization

$$
A = QR,
$$

where $Q$ is orthogonal and $R$ is upper triangular.

Starting from

$$
Ax=b,
$$

substituting $A=QR$ gives

$$
QRx=b.
$$

Multiplying both sides by $Q^T$ gives

$$
Q^TQRx = Q^Tb.
$$

Since $Q$ is orthogonal,

$$
Q^TQ=I,
$$

so

$$
Rx=Q^Tb.
$$

Define

$$
y=Q^Tb.
$$

Then the system becomes

$$
Rx=y.
$$

Finally, solve the upper-triangular system

$$
Rx=y
$$

for $x$ using back substitution.

The computed solution can then be compared with the known solution $x_{\text{true}}$. A correct QR factorization should produce a value of $x$ very close to

$$
\begin{bmatrix}
1 \\
2 \\
3
\end{bmatrix}.
$$