# FEM
code matlab for FEM
The problem:

$\Omega$ is a open bounded domain with a regulier boundary. A is a tensor which is uniformly bounded
\begin{equation*}
    \exists C > 0, \forall(x,y)\in \Omega, \forall i,j, A_{i,j}(x,y) \leq C
\end{equation*}
and satisfies the hypothesis of uniform coerciveness.
\begin{equation*}
    \exists c>0, \forall(x,y)\in \Omega, \forall \xi \in \mathbb{R}^2, A(x,y)\xi \cdot \xi \ge|\xi|^2
\end{equation*}
and $f \in L^2(\Omega)$.

with Neumann boundary condition:

\\\textsl{Trouver $u \in H^1(\Omega)$  telle que:}
\begin{equation} 
  \left\{
    \begin{aligned}
      {&u-\nabla\cdot(A(x,y)\nabla u) = f\quad dans \quad \Omega \\
      &A(x,y)\nabla u \cdot n =0\quad sur \quad \partial\Omega\\
      \end{aligned}
    \right.
\end{equation}
