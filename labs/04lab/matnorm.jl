# Spectral_Norm.jl

using LinearAlgebra
#using Plots

 

    function spectral_stat(A,sim_size)

        maxEigVec = 0.0
        sim_size = 1000000
    
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
            if maxEigVec < t    # checks of the new value is greater than the current max
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
        w = normalized_x()  # generate a random vector of length 4 with values in [-1,1)
        y = copy(w)         # copy creates a new/separate object
        lambda1 = 0.0
        for j = 1:30        # Tried a few different values for the number of iterations, 30 seems to have good accuracy.  
            w=B*y
            p=argmax(abs.(w))
            lambda1=(w'*y)/(y'*y)    # calculate the Rayleigh quotient
            #println("lambda1^($j)=$lambda1")
            y=w/w[p]                 # normalize w by the p-th element of w to avoid overflow and assign it to y for the next iteration
        end

        #println("Spectral Norm using power methon ")
        return sqrt(lambda1) 


    end

    function spectral_opnorm(A)

        #println("Spectral Norm Julia built-in opnorm(A)")
        return opnorm(A)
    
    end


function main()

    # A1 is used for testing the functions, A is used for the final results
    A1 = [
         0.23   0.29  -1.85  -3.48;
         4.24   4.00   2.36   4.82;
         4.56  -2.97   1.90  -0.51;
        -1.18  -2.55  -4.86   0.50;
         3.31  -2.04  -2.58   0.35
    ]

    # Genereatd from https://fractal.math.unr.edu/~ejolson/466-23/specnorm/snmatrix.cgi
    A= [ -1.88 -1.50 -2.77  0.81;
          2.21  3.92  1.70  2.96;
         -2.65 -3.62 -1.02 -1.62;
         -0.43 -0.65 -0.52  0.66;
          4.20 -4.80  3.29  1.84 ]


    sim_size = 1_000_000

    stat_result   = spectral_stat(A, sim_size)
    power_result  = spectral_power(A)
    opnorm_result = spectral_opnorm(A)

    println("The Matrix A is:")
    display(A)
    println("")

    println("    Statistical approach = ", stat_result)
    println("    Power method         = ", power_result)
    println("    opnorm               = ", opnorm_result)

end

# added for debug using REPL and VSCode in parallel 
# runs the functions only if running in the REPL 
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

