# ==============================================================================
# M466 Lab - System Information
#
# Machine:         Meridian1
# OS:              Windows 11 / Ubuntu WSL
#
# CPU:             Intel Core Ultra 5 245HX
# WSL CPUs:        14
# Architecture:    x86_64
# RAM:             31 GiB available to WSL
#
# GPU:             NVIDIA RTX PRO 1000 Blackwell Generation Laptop GPU
# Architecture:    NVIDIA Blackwell
# CUDA Cores:      2,560
# Tensor Cores:    80 (5th generation)
# RT Cores:        20 (4th generation)
# GPU Memory:      8 GB GDDR7 (~7.96 GiB usable)
# Memory Bus:      128-bit
# Memory Bandwidth: 384 GB/s
# Compute Cap.:    12.0 (sm_120)
# GPU Power Cap:   110 W
#
# NVIDIA Driver:   596.53
# Driver CUDA:     13.2
# CUDA.jl Runtime: 13.4.0
#
# Julia:           1.13.0
# LLVM:            20.1.8
# Julia Project:   ~/M466
#------------------------------------------------------------------------------
# Commands used to obtain system information
#
# WSL/Linux terminal:
#
#   lscpu
#       CPU model, architecture, and CPUs available to WSL
#
#   free -h
#       RAM and swap available to WSL
#
#   nvidia-smi
#       GPU model, VRAM, power limit, driver, and supported CUDA version
#
# Julia REPL:
#
#   versioninfo()
#       Julia, OS, CPU, LLVM, and thread information
#
#   using CUDA
#   dev = CUDA.device()
#   CUDA.capability(dev)
#       GPU device and compute capability
#
#   CUDA.versioninfo()
#       CUDA runtime, driver, libraries, Julia CUDA packages, and GPU information
#
# CUDA core count:
#   Obtained from published specifications for the exact GPU model;
#   nvidia-smi and CUDA.versioninfo() do not directly report CUDA core count.
# ==============================================================================



# flops.jl

using LinearAlgebra
using CUDA

function main()

    warmpup = 100
    loop = [1000, 2000, 3000, 4000]
    trials = 5
    T = Float64
    
    println("n = ", warmpup, " warm-up iterations")
    println("precision = ", T)
    println()

    # --------------------------------------------------
    # Warm-up - Initialize GPU and compile kernels
    # --------------------------------------------------

    nw = 100

    A = rand(T, nw, nw)
    b = A * ones(T, nw)
    A \ b

    Ag = CuArray(A)
    bg = CuArray(b)

    CUDA.@sync begin
        Ag \ bg
    end

    # --------------------------------------------------
    # CPU - LAPACK+BLAS
    # --------------------------------------------------

    println("CPU")
    BLAS.set_num_threads(7)  # note: CPU has 14 logical cores, but we limit BLAS to 7 threads.
                             # Setting BLAS threds to 7 had better performance than 14 threads on this CPU for this problem size.
    println("BLAS threads = ", BLAS.get_num_threads())

    cpu_flmax = 0.0

    for j in 1:length(loop)
        n = loop[j]
        println("Matrix size n = ", n)

        for i = 1:trials

            A = rand(T, n, n)
            xtrue = ones(T, n)
            b = A * xtrue

            Tn = @elapsed x = A \ b

            flop_count(n) = (2/3) * n^3 + 2 * n^2  # Estimate of floating-point operations for LU factorization and triangular solves
            flops = flop_count(n) / Tn
            cpu_flmax = max(cpu_flmax, flops)

            error = norm(x - xtrue) / norm(xtrue)

            
            println("Trial ", i)
            println("  Tn     = ", Tn)
            println("  Error  = ", error)
            println("  GFLOPS = ", flops / 1e9)
        end

        println("Maximum CPU GFLOPS = ", cpu_flmax / 1e9)

    end 

    # --------------------------------------------------
    # GPU - cuSOLVER+cuBLAS
    # --------------------------------------------------

    println()
println("GPU")

# Store the maximum GPU floating-point performance observed
# over all benchmark trials (FLOPS/second).
gpu_flmax = 0.0

# Repeat the experiment several times to observe performance variability.
for j in 1:length(loop)
        n = loop[j]
        println("Matrix size n = ", n)
    for i = 1:trials

        # --------------------------------------------------
        # Construct the linear system Ax = b on the CPU
        # --------------------------------------------------

        # Generate a random n x n matrix using precision T
        # (Float32 or Float64).
        A = rand(T, n, n)

        # Choose the exact solution to be a vector of ones.
        xtrue = ones(T, n)

        # Construct b so that xtrue is the known solution:
        #
        #       A*xtrue = b
        #
        b = A * xtrue


        # --------------------------------------------------
        # Transfer data from CPU memory to GPU memory
        # --------------------------------------------------

        # CuArray copies the CPU arrays into GPU VRAM.
        Ag = CuArray(A)
        bg = CuArray(b)
        xtrue_g = CuArray(xtrue)

        # Wait for any previously issued GPU operations to finish
        # before starting the timer.
        CUDA.synchronize()


        # --------------------------------------------------
        # Solve Ax = b on the GPU and measure execution time
        # --------------------------------------------------

        Tn = @elapsed begin

            # Solve the dense linear system on the GPU.
            #
            # For a general dense matrix this involves an LU
            # factorization followed by triangular solves.
            # CUDA.jl uses NVIDIA's GPU linear algebra stack
            # (e.g. cuSOLVER/cuBLAS) for these operations.
            xg = Ag \ bg

            # GPU operations are asynchronous: the CPU can continue
            # before the GPU has actually finished.
            #
            # synchronize() forces Julia to wait until the GPU solve
            # is complete, ensuring Tn measures the actual solve time.
            CUDA.synchronize()
        end


        # --------------------------------------------------
        # Calculate achieved floating-point performance
        # --------------------------------------------------

        # Estimate the number of floating-point operations required
        # for the LU factorization and triangular solves:
        #
        #       (2/3)n^3 + 2n^2
        #
        # Divide by elapsed time to obtain FLOPS/second.
        flop_count(n) = (2/3) * n^3 + 2 * n^2  # Estimate of floating-point operations for LU factorization and triangular solves
        flops = flop_count(n) / Tn

        # Keep the highest observed performance across all trials.
        gpu_flmax = max(gpu_flmax, flops)


        # --------------------------------------------------
        # Calculate numerical error
        # --------------------------------------------------

        # Relative 2-norm error:
        #
        #       ||x_computed - x_true||_2
        #       -------------------------
        #              ||x_true||_2
        #
        # Both vectors are CuArrays, so this calculation is also
        # performed using GPU data.
        error = norm(xg - xtrue_g) / norm(xtrue_g)


        # --------------------------------------------------
        # Report results for this trial
        # --------------------------------------------------


        println("Trial ", i)
        println("  Tn     = ", Tn)
        println("  Error  = ", error)

        # Convert FLOPS/sec to GFLOPS:
        #
        #       1 GFLOP = 10^9 FLOPs
        #
        println("  GFLOPS = ", flops / 1e9)
    end
    println("Maximum GPU GFLOPS = ", gpu_flmax / 1e9)
    end
end

main()