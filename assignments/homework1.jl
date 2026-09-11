function decode_bfloat16(bits)
    # Remove spaces so either format works:
    # "0 10000000 1001001"
    # "0100000001001001"
    println("bfloat input:  ", bits)
    bits = replace(bits, " " => "")

    # Split the 16 bits into fields
    sign_bit     = bits[1]
    exponent_bits = bits[2:9]
    fraction_bits = bits[10:16]

    println("Sign bits:     ", sign_bit)
    println("Exponent bits: ", exponent_bits)
    println("Fraction bits: ", fraction_bits)

    # ----- Sign -----
    s = sign_bit == '0' ? 1.0 : -1.0
    println("\nSign = ", s)

    # ----- Exponent -----
    E = parse(Int, exponent_bits; base=2)
    e = E - 127

    println("Stored exponent E = ", E)
    println("Actual exponent e = ", E, " - 127 = ", e)

    # ----- Fraction -----
    f = 0.0

    println("\nFraction calculation:")

    for i in 1:length(fraction_bits)
        bit = fraction_bits[i]
        weight = 2.0^(-i)

        println("bit ", i,
                " = ", bit,
                ", weight = 2^(-", i, ") = ", weight)

        if bit == '1'
            f += weight
            println("    add ", weight, " -> f = ", f)
        end
    end

    # Normalized bfloat16 has an implicit leading 1
    significand = 1.0 + f

    println("\nFraction f = ", f)
    println("Significand 1 + f = ", significand)

    # ----- Put everything together -----
    value = s * significand * 2.0^e

    println("\nFinal calculation:")
    println("x = (-1)^s × (1 + f) × 2^e")
    println("x = ", s, " × ", significand, " × 2^", e)
    println("x = ", value)

    return value
end

decode_bfloat16("0 10000000 1001001")
println("\n,\n")
decode_bfloat16("0 01111011 1001101")