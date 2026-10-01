# Spectral_Norm.jl

using LinearAlgebra
using Plots

 

    function spectral_stat(A,sim_size)

        maxEigVec = 0.0
    
        function normalized_x() 
            while true            # use while loop to ensure x_norm is never zero
                x=2 * rand(4) .-1 # mult by 2 shifts interval to right, subtract 1 to generate interval on [-1,1)
                x_norm = norm(x)
                if x_norm > 0
                    return x/x_norm
                end
            end
        end

        for i=1:sim_size
            z = normalized_x()
            t=norm(A*z)
            if maxEigVec < t
                maxEigVec = t
            end
        end

        #println("Spectral Norm using Monte Carlo simulation")
        return maxEigVec
    
    end

    function spectral_power(A)

        function normalized_x() 
            while true            # use while loop to ensure x_norm is never zero
                x=2 * rand(4) .-1 # mult by 2 shifts interval to right, subtract 1 to generate interval on [-1,1)
                x_norm = norm(x)
                if x_norm > 0
                    return x/x_norm
                end
            end
        end

        B = A'A
        w = normalized_x()
        y = copy(w)    # copy creates a new/separate object
        lambda1 = 0.0
        for j = 1:30
            w=B*y
            p=argmax(abs.(w))
            lambda1=(w'*y)/(y'*y)
            #println("lambda1^($j)=$lambda1")
            y=w/w[p]
        end

        #println("Spectral Norm using power methon ")
        return sqrt(lambda1) 


    end

    function spectral_opnorm(A)

        #println("Spectral Norm Julia built-in opnorm(A)")
        return opnorm(A)
    
    end


using LinearAlgebra


function spectral_stat(A, sim_size)

    function normalized_x()
        while true
            x = 2 * rand(4) .- 1
            x_norm = norm(x)

            if x_norm > 0
                return x / x_norm
            end
        end
    end

    maxEigVec = 0.0

    for i = 1:sim_size
        z = normalized_x()
        t = norm(A * z)

        if maxEigVec < t
            maxEigVec = t
        end
    end

    return maxEigVec
end


function spectral_power(A)

    function normalized_x()
        while true
            x = 2 * rand(4) .- 1
            x_norm = norm(x)

            if x_norm > 0
                return x / x_norm
            end
        end
    end

    B = A' * A
    w = normalized_x()
    y = copy(w)

    lambda1 = 0.0

    for j = 1:30
        w = B * y
        p = argmax(abs.(w))
        lambda1 = (w' * y) / (y' * y)
        y = w / w[p]
    end

    return sqrt(lambda1)
end


function spectral_opnorm(A)
    return opnorm(A)
end


function main()

    A = [
         0.23   0.29  -1.85  -3.48;
         4.24   4.00   2.36   4.82;
         4.56  -2.97   1.90  -0.51;
        -1.18  -2.55  -4.86   0.50;
         3.31  -2.04  -2.58   0.35
    ]

    sim_size = 1_000_000

    stat_result   = spectral_stat(A, sim_size)
    power_result  = spectral_power(A)
    opnorm_result = spectral_opnorm(A)

    println("    Statistical approach = ", stat_result)
    println("    Power method         = ", power_result)
    println("    opnorm               = ", opnorm_result)

end

# added for debug using REPL and VSCode in parallel 
# runs the functions only if running in the REPL 
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

