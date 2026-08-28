"""
    j2pdf(name)

Combine `name.jl` and `name.out` into `name.pdf`.
The files are assumed to be in the current working directory.
"""
function j2pdf(name::String)

    jlfile  = name * ".jl"
    outfile = name * ".out"
    texfile = name * ".tex"

    # Make sure the source files exist
    isfile(jlfile)  || error("File not found: $jlfile")
    isfile(outfile) || error("File not found: $outfile")

    # Read the Julia source and output
    jltext  = read(jlfile, String)
    outtext = read(outfile, String)

    # Construct the LaTeX document
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

\section*{Program Output}

\begin{lstlisting}
""" * outtext * raw"""
\end{lstlisting}

\end{document}
"""

    # Write temporary LaTeX file
    write(texfile, tex)

    # Compile PDF
    run(`pdflatex -interaction=nonstopmode $texfile`)

    println("Created $(name).pdf")
end