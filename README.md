# Alternating reciprocal cubes of Fibonacci and Lucas sequences

**twinheterodox · 18 November 2026 · Number theory**

The alternating sum of the reciprocal cubes of the Fibonacci numbers is irrational. The theorem holds for the first-kind Lucas sequence of every real quadratic unit greater than one, including the Pell sequence. If one retains the terms with indices $k,2k,3k,\ldots$ and alternates their signs with the new index, the sum is irrational for each fixed positive integer $k$.

[Read the paper (PDF)](output/pdf/alternating-lucas-cubes.pdf) · [LaTeX source](paper/alternating-lucas-cubes.tex)

## Reproduce the finite checks

From the repository root, run

~~~sh
python3 -m pip install -r checks/requirements.txt
python3 -B checks/verify.py
~~~

The independent audit programs recompute every coefficient of the formal determinant for sizes 2 through 5 and both signs of the norm. They also compare two evaluations at both real embeddings for the Fibonacci, Pell and $(P,Q)=(3,1)$ sequences. The exact norm margin is $296/11025$. The general irrationality result is proved in the manuscript. [Reference results](checks/reference-results.json) · [Independent audit scope](review/independent-audit.json)

## Build the manuscript

From the same directory, run

~~~sh
python3 -B tools/build_paper.py
~~~

This uses pdfLaTeX, Latin Modern fonts, and Poppler's pdfinfo. It runs pdfLaTeX twice, checks the final log, and writes the reader PDF under output/pdf. The recorded build environment is in [the build record](paper/build-record.json).

## Sources and licensing

The proof uses the formal cyclotomic divisor from Krisman's ordinary reciprocal-cube theorem. The signed conjugate pole moments are evaluations of one holomorphic function, which allows a contour integral to estimate their sum. The [source notes](sources/inputs.json) identify the exact arithmetic input, its antecedents, and the classical odd-index transcendence theorem.

The manuscript and accompanying explanatory text are offered under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The repository's computer programs are distributed under [Apache-2.0](LICENSE).
