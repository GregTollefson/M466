## Running `matnorm.jl`

From the directory containing `matnorm.jl`, run:

```bash
julia matnorm.jl > matnorm.out
```

This runs the program and redirects the terminal output to:

```text
matnorm.out
```

To view the saved output:

```bash
cat matnorm.out
```

After running the program, the directory should contain:

```text
matnorm.jl
matnorm.out
```

## Creating the PDF

Start the Julia REPL from the directory containing `matnorm.jl` and `matnorm.out`:

```bash
julia
```

Load the `j2pdf` utility:

```julia
include("../../utils/j2pdf.jl")
```

Then create the PDF:

```julia
j2pdf("matnorm")
```

This generates:

```text
matnorm.pdf
```
The pdf can be viewed on WSL/Linux using `evince`

```
evince matnorm.pdf
```

The PDF contains:

- the Julia source code from `matnorm.jl`
- the program output from `matnorm.out`

If the relative path to `j2pdf.jl` is different, adjust the `include(...)` path accordingly.

# Lab 04 — The Spectral Norm

## Overview

The goal of this lab is to compute the spectral norm of a matrix

$$
A \in \mathbb{R}^{m \times n}
$$

using three different methods:

1. A statistical approach using random unit vectors.
2. A spectral approach using the scaled power method.
3. Julia's built-in `opnorm` function.

The three methods are implemented as

```julia
spectral_stat(A, sim_size)
spectral_power(A)
spectral_opnorm(A)
```

and are called from the main program.

---

## Matrix Norm

The matrix norm induced by a vector norm is defined as

$$
\|A\|
=
\max
\left\{
\|Ax\| : \|x\| \leq 1
\right\}.
$$

For this lab, the vector norm is the Euclidean or $L_2$ norm,

$$
\|x\|_2
=
\sqrt{x^T x}
=
\sqrt{\sum_{i=1}^{n}x_i^2}.
$$

The corresponding matrix norm is therefore

$$
\|A\|_2
=
\max
\left\{
\|Ax\|_2 : \|x\|_2 \leq 1
\right\}.
$$

This matrix norm is called the **spectral norm**.

---

# Statistical Approach

## Restricting the Search to Unit Vectors

Suppose

$$
x \neq 0
$$

and

$$
\|x\| < 1.
$$

Define the corresponding normalized vector

$$
z = \frac{x}{\|x\|}.
$$

Then

$$
\|z\| = 1.
$$

Now consider

$$
\|Az\|.
$$

Substituting

$$
z = \frac{x}{\|x\|}
$$

gives

$$
\|Az\|
=
\left\|
A\frac{x}{\|x\|}
\right\|.
$$

Because $1/\|x\|$ is a scalar,

$$
\|Az\|
=
\frac{1}{\|x\|}\|Ax\|.
$$

Since

$$
\|x\| < 1,
$$

we have

$$
\frac{1}{\|x\|} > 1.
$$

Therefore,

$$
\|Az\| \geq \|Ax\|.
$$

This means the maximum stretch occurs on the boundary of the unit ball, so the spectral norm can be written as

$$
\boxed{
\|A\|_2
=
\max
\left\{
\|Az\|_2 : \|z\|_2 = 1
\right\}.
}
$$

Therefore, it is sufficient to search over unit vectors.

---

## Generating Random Directions

Julia's

```julia
rand(n)
```

generates a vector

$$
V \in \mathbb{R}^n
$$

whose components satisfy

$$
V_i \in [0,1).
$$

In two dimensions, these vectors would always lie in the first quadrant.

To generate positive and negative components, use the transformation

$$
X_i = 2V_i - 1.
$$

Since

$$
0 \leq V_i < 1,
$$

multiplying by $2$ gives

$$
0 \leq 2V_i < 2,
$$

and subtracting $1$ gives

$$
-1 \leq 2V_i - 1 < 1.
$$

Therefore,

$$
X_i \in [-1,1).
$$

In Julia this is

```julia
x = 2 * rand(4) .- 1
```

for a vector in $\mathbb{R}^4$.

The vector is then normalized:

$$
z = \frac{x}{\|x\|}.
$$

By construction,

$$
\|z\| = 1.
$$

The zero vector must be avoided because

$$
\frac{x}{\|x\|}
$$

would be undefined when

$$
\|x\| = 0.
$$

---

## Statistical Estimate

For each random unit vector $z$, calculate

$$
t = \|Az\|_2.
$$

The statistical method repeats this calculation for many different random directions and keeps the largest observed value:

$$
A_{\text{bound}}
=
\max_{1\leq i\leq N}
\|Az_i\|_2.
$$

For this lab,

$$
N = 1,000,000.
$$

Because only a finite collection of directions is tested,

$$
A_{\text{bound}}
\leq
\|A\|_2.
$$

Thus, the statistical result provides an approximation to the spectral norm and normally approaches the true value from below.

Because the vectors are generated randomly, the statistical result can be slightly different each time the program is run.

---

# Spectral Approach

The second approach uses the relationship between the spectral norm and the eigenvalues of

$$
B=A^TA.
$$

The spectral norm satisfies

$$
\boxed{
\|A\|_2
=
\sqrt{\lambda_{\max}(A^TA)}
}
$$

where

$$
\lambda_{\max}(A^TA)
$$

is the largest eigenvalue of $A^TA$.

For a matrix

$$
A \in \mathbb{R}^{5\times4},
$$

the transpose has dimensions

$$
A^T \in \mathbb{R}^{4\times5}.
$$

Therefore,

$$
B=A^TA
$$

has dimensions

$$
(4\times5)(5\times4)
=
4\times4.
$$

Thus,

$$
B\in\mathbb{R}^{4\times4}.
$$

---

## Why $A^TA$?

The matrix

$$
B=A^TA
$$

is symmetric because

$$
(A^TA)^T
=
A^TA.
$$

It is also positive semidefinite because for any vector $x$,

$$
x^TA^TAx
=
(Ax)^T(Ax)
=
\|Ax\|_2^2
\geq 0.
$$

Therefore, the eigenvalues of $A^TA$ are real and nonnegative.

If

$$
\lambda_1
=
\lambda_{\max}(A^TA),
$$

then

$$
\boxed{
\|A\|_2=\sqrt{\lambda_1}.
}
$$

---

# Scaled Power Method

The largest eigenvalue of

$$
B=A^TA
$$

is approximated using the scaled power method.

Begin with a nonzero initial vector

$$
w^{(0)} \in \mathbb{R}^4.
$$

A random normalized vector is used in the program.

Set

$$
y^{(0)}=w^{(0)}.
$$

For each iteration

$$
j=1,2,\ldots,30,
$$

perform the following steps.

---

## Step 1 — Matrix-Vector Multiplication

Calculate

$$
w^{(j)}
=
By^{(j-1)}.
$$

In Julia:

```julia
w = B * y
```

Repeated multiplication by $B$ causes the vector to increasingly align with the eigenvector corresponding to the largest eigenvalue.

To see why, suppose the initial vector can be expressed as a linear combination of the eigenvectors of $B$:

$$
y^{(0)}
=
c_1v_1+c_2v_2+\cdots+c_nv_n.
$$

Because

$$
Bv_i=\lambda_i v_i,
$$

multiplying once by $B$ gives

$$
By^{(0)}
=
c_1\lambda_1v_1+
c_2\lambda_2v_2+
\cdots+
c_n\lambda_nv_n.
$$

After $k$ multiplications,

$$
B^ky^{(0)}
=
c_1\lambda_1^kv_1+
c_2\lambda_2^kv_2+
\cdots+
c_n\lambda_n^kv_n.
$$

If

$$
|\lambda_1|>|\lambda_2|,
$$

then as $k$ increases, the term containing

$$
\lambda_1^k
$$

dominates.

Therefore, the direction of the vector approaches the dominant eigenvector

$$
v_1.
$$

---

## Step 2 — Find the Largest Component

Find the index $p$ such that

$$
|w_p^{(j)}|
=
\max_i |w_i^{(j)}|.
$$

In Julia:

```julia
p = argmax(abs.(w))
```

The expression

```julia
abs.(w)
```

takes the absolute value of every component of $w$, while

```julia
argmax(...)
```

returns the index of the largest value.

---

## Step 3 — Estimate the Largest Eigenvalue

The eigenvalue approximation is calculated as

$$
\lambda_1^{(j)}
\approx
\frac{
(w^{(j)})^T y^{(j-1)}
}{
(y^{(j-1)})^T y^{(j-1)}
}.
$$

Since

$$
w^{(j)}
=
By^{(j-1)},
$$

this can be written as

$$
\lambda_1^{(j)}
\approx
\frac{
(By)^Ty
}{
y^Ty
}.
$$

Because $B$ is symmetric,

$$
(By)^Ty
=
y^TBy.
$$

Therefore,

$$
\boxed{
\lambda_1^{(j)}
\approx
\frac{y^TBy}{y^Ty}
}
$$

which is the Rayleigh quotient.

In Julia:

```julia
lambda1 = (w' * y) / (y' * y)
```

As $y$ approaches the dominant eigenvector, this estimate approaches the largest eigenvalue

$$
\lambda_1.
$$

---

## Step 4 — Scale the Vector

The vector is rescaled using its largest-magnitude component:

$$
y^{(j)}
=
\frac{w^{(j)}}{w_p^{(j)}}.
$$

In Julia:

```julia
y = w / w[p]
```

After scaling,

$$
|y_i|\leq1
$$

for every component, and the component with the largest absolute value has magnitude $1$.

This scaling prevents the entries in the power iteration from becoming extremely large or extremely small while preserving the direction of the vector.

---

# Recovering the Spectral Norm

After 30 iterations, the power method produces an approximation to

$$
\lambda_{\max}(A^TA).
$$

The spectral norm is then

$$
\boxed{
\|A\|_2
=
\sqrt{\lambda_1}
}
$$

so the function returns

```julia
sqrt(lambda1)
```

rather than `lambda1` itself.

The vector $y$ approximates the **direction of the dominant eigenvector**, while `lambda1` approximates the corresponding **dominant eigenvalue**.

These are different quantities:

$$
y
\approx
v_1
$$

and

$$
\lambda_1
\approx
\lambda_{\max}(A^TA).
$$

The spectral norm is obtained from the eigenvalue:

$$
\|A\|_2
=
\sqrt{\lambda_1}.
$$

---

# Julia Built-In Spectral Norm

Julia provides the function

```julia
opnorm(A)
```

which directly computes the operator norm of the matrix.

For the Euclidean norm, this is the spectral norm:

$$
\boxed{
\operatorname{opnorm}(A)
=
\|A\|_2.
}
$$

The built-in implementation uses a more sophisticated numerical algorithm than the simple power method.

The `opnorm` result therefore provides a useful comparison for checking the statistical and power-method implementations.

---

# Comparison of the Three Methods

For the example matrix used during development, the three methods produced approximately

$$
\begin{aligned}
\text{Statistical approach}
&\approx 9.2023,\\
\text{Scaled power method}
&\approx 9.2047527188,\\
\text{Julia opnorm}
&\approx 9.2047527188.
\end{aligned}
$$

The statistical result is close but slightly smaller because only a finite number of randomly selected directions were tested.

After 30 iterations, the scaled power method converges very closely to the same result as `opnorm`.

Thus,

$$
\boxed{
\|A\|_2
\approx
9.2047527188
}
$$

for the example matrix.

The agreement among the three methods provides a useful consistency check on the calculations.

---

# Program Structure

The program is organized around three functions:

```julia
spectral_stat(A, sim_size)
spectral_power(A)
spectral_opnorm(A)
```

The `main()` function defines the matrix, sets the statistical simulation size, calls each method, and prints the results.

Conceptually,

$$
A
\longrightarrow
\begin{cases}
\texttt{spectral\_stat}\\
\texttt{spectral\_power}\\
\texttt{spectral\_opnorm}
\end{cases}
\longrightarrow
\texttt{main()}
\longrightarrow
\text{output}.
$$

The program can be executed from the shell with

```bash
julia MatrixNorms.jl
```

and the output can be redirected to a file with

```bash
julia MatrixNorms.jl > MatrixNorms.out
```

The program and its output can then be combined into the final submission PDF using the course `j2pdf` utility.