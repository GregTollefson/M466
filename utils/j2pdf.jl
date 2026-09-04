"""
    j2pdf(name)

Convert `name.jl` into `name.pdf`.

If `name.out` exists, include it as a "Program Output" section.
Otherwise, generate a PDF containing only the Julia source code.

The files are assumed to be in the current working directory.
"""
function j2pdf(name::String)

    jlfile  = name * ".jl"
    outfile = name * ".out"
    texfile = name * ".tex"

    # The Julia source file is required
    isfile(jlfile) || error("File not found: $jlfile")

    # Read the Julia source
    jltext = read(jlfile, String)

    # Start constructing the LaTeX document
    tex = raw"""
\documentclass[11pt]{article}
\usepackage[margin=0.75in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{listings}

\lstset{
    basicstyle=\ttfamily\small,
    breaklines=true,
    columns=fullflexible,
    frame=single
}

\begin{document}

\section*{Julia Source Code}

\begin{lstlisting}
""" * jltext * raw"""
\end{lstlisting}
"""

    # Add program output only if the .out file exists
    if isfile(outfile)
        outtext = read(outfile, String)

        tex *= raw"""

\section*{Program Output}

\begin{lstlisting}
""" * outtext * raw"""
\end{lstlisting}
"""
    end

    # Finish the LaTeX document
    tex *= raw"""

\end{document}
"""

    # Write temporary LaTeX file
    write(texfile, tex)

    # Compile PDF
    run(`pdflatex -interaction=nonstopmode $texfile`)

    println("Created $(name).pdf")
end