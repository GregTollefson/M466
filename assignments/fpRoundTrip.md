# Decimal ↔ Float64 Round-Trip Precision

IEEE 754 double-precision floating-point numbers (`Float64`) use a
53-bit significand. This corresponds to roughly

$$
53\log_{10}(2) \approx 15.95
$$

decimal digits of precision.

However, there are **two different round-trip questions**, and they have
different answers.

---

## 1. Float64 → Decimal → Float64

Start with a floating-point number

$$
x^* \in \mathbb{F}_{64}.
$$

Convert it to a decimal representation containing $m$ significant digits:

$$
x^*
\longrightarrow
\text{decimal with }m\text{ significant digits}.
$$

Then convert that decimal number back to the nearest `Float64`:

$$
\text{decimal}
\longrightarrow
y^* \in \mathbb{F}_{64}.
$$

The question is:

> What is the smallest number of decimal digits $m$ that guarantees

$$
x^* = y^*
$$

for every `Float64` number?

For IEEE double precision,

$$
\boxed{m=17}.
$$

Thus, **17 significant decimal digits are sufficient to uniquely identify
any `Float64` value**.

In other words,

$$
\boxed{
\text{Float64}
\rightarrow
17\text{-digit decimal}
\rightarrow
\text{Float64}
}
$$

is guaranteed to recover the original floating-point number.

---

## Why 17 Digits?

A `Float64` has 53 bits of significand precision:

$$
53\log_{10}(2)
\approx
15.95.
$$

Although this corresponds to roughly 16 decimal digits of numerical
precision, slightly more decimal digits are needed to distinguish every
possible neighboring `Float64` value.

The round-trip guarantee therefore requires

$$
\boxed{17\text{ significant decimal digits}}.
$$

This value is sometimes referred to as `max_digits10`.

---

## Julia's Printing Behavior

When Julia prints a `Float64`, it generally does **not** simply print
17 digits every time.

Instead, Julia prints the **shortest decimal representation that will
round back to exactly the same floating-point number**.

For example,

```julia
x = 0.1
println(x)
```

prints

```text
0.1
```

even though the actual binary `Float64` representation is not exactly
the mathematical decimal number $0.1$.

The decimal string `0.1` is sufficient because converting it back to
`Float64` produces the same floating-point value.

Thus Julia attempts to make the round trip

$$
x^*
\rightarrow
\text{decimal string}
\rightarrow
x^*
$$

an identity while printing as few digits as necessary.

---

## 2. Decimal → Float64 → Decimal

Now consider the opposite direction.

Start with a decimal number $x$ containing $m$ significant digits:

$$
x
\longrightarrow
\text{nearest Float64}.
$$

Then convert that floating-point number back to a decimal representation
with $m$ significant digits:

$$
\text{Float64}
\longrightarrow
y.
$$

The question is now:

> What is the largest value of $m$ such that

$$
x=y
$$

is guaranteed?

For IEEE double precision,

$$
\boxed{m=15}.
$$

Thus, a decimal number with at most 15 significant digits can be
converted to `Float64` and then back to a decimal number with 15
significant digits without changing the value, assuming the number is
within the appropriate representable range.

Therefore,

$$
\boxed{
15\text{-digit decimal}
\rightarrow
\text{Float64}
\rightarrow
15\text{-digit decimal}
}
$$

is guaranteed to preserve the original decimal value.

---

## Summary

The two directions give different guarantees:

| Round Trip | Guaranteed Significant Decimal Digits |
|---|---:|
| `Float64 → Decimal → Float64` | 17 |
| `Decimal → Float64 → Decimal` | 15 |

So,

$$
\boxed{
\text{Float64}
\rightarrow
17\text{ decimal digits}
\rightarrow
\text{same Float64}
}
$$

while

$$
\boxed{
15\text{-digit decimal}
\rightarrow
\text{Float64}
\rightarrow
\text{same 15-digit decimal}
}
$$

---

## Important Interpretation

The statement that 17 digits are sufficient for a `Float64` round trip
does **not** mean that double precision has 17 decimal digits of numerical
precision.

A `Float64` has approximately

$$
\boxed{15\text{--}16\text{ decimal digits of precision}}.
$$

The role of the 17th digit is different:

> It provides enough decimal information to uniquely distinguish one
> `Float64` number from neighboring `Float64` numbers.

So the useful distinction is

$$
\boxed{
\begin{aligned}
15\text{--}16 \text{ digits}
&\quad\text{: numerical precision},\\[4pt]
17 \text{ digits}
&\quad\text{: guaranteed Float64 round-trip identification}.
\end{aligned}
}
$$