# Double mutation-response scanning

$$
\nonumber
\newcommand{\fvec}{\mathbf{f}}
\newcommand{\evec}{\mathbf{e}}
\newcommand{\avec}{\mathbf{a}}
\newcommand{\amat}{\mathbf{A}}
\newcommand{\umat}{\mathbf{U}}
\newcommand{\kmat}{\mathbf{K}}
\newcommand{\cmat}{\mathbf{C}}
\newcommand{\bmat}{\mathbf{B}}
\newcommand{\vvec}{\mathbf{v}}
\newcommand{\uvec}{\mathbf{u}}
\newcommand{\x}{\mathbf{x}}
\newcommand{\y}{\mathbf{y}}
$$
## Simpler problems
### Quadratic form

Let $\amat$ be a square matrix then find maximum of:
$$
\x'\amat\x = x_iA_{ij}x_j
$$
subject to the constraint: 
$$
\x'\x = \alpha^2
$$
The lagrange function is:
$$
L(\x,\lambda) = \x'\amat\x - \lambda (\x'\x -\alpha^2)
$$
then:
$$
\begin{array}
dL/d\x = \amat\x +(x'\amat)'-2\lambda\x  = 0\\
dL/d\lambda = -(\x'\x - \alpha^2) = 0
\end{array}
$$
the solution is:
$$
\frac{(\amat + \amat')}{2}\x=\lambda\x \\
\x'\x = \alpha^2
$$
Therefore, the solution os an eigenvector of a symmetric matrix. If we want the maximum, then it is the eigenvector with maximum eigenvalue. Plus, we normalize this vector so that $|\x|^2 = \alpha^2$. 

To finish, the value of the quadratic form at the maximum is:
$$
\x'A\x=\x'(\amat_s+\amat_u)\x = \alpha^2\lambda_{\max}+\x'\amat_u\x
$$
but the second term is zero because $\amat_u$ is anti-symmetric, thus:
$$
\max(\x'\amat\x) = \alpha^2\lambda_\max \label{xAx_max}
$$


### Maximum of bilinear form, given $\x$
Let $\amat be an $N \times M$ matrix, and $\x$ and $\y$ be column vectors of size $N$ and $M$ respectively. Then find $\y$ that satisfies:
$$
\x'\amat\y=\max \\
\y'\y = \alpha^2
$$
then, we need to maximize:
$$
L(\y,\lambda) = \x'\amat\y - \lambda (\y'\y -\alpha^2)
$$


Deriving with respect to $\y$ and $\lambda$:
$$
\amat'\x - 2\lambda\y = 0; \y'\y-\alpha^2=0
$$
with solutions:
$$
\y = \alpha \frac{\amat'\x}{|\amat'\x|}\\ \label{ymax}
\lambda =\frac{|\amat' \x|}{2\alpha}
$$
Replacing this into the original quadratic form:
$$
\max(\x'\amat\y) = \alpha |A'\x| \label{xAy_max_y}
$$
Note that if $\amat = \mathbf{I}$, $\y$ is colinear with $\x$ and $\max(\x'\y)$ is just the product of the length of $\y$ times the length of $\x$. 
### Maximum of maximum

Now, let's find $\x$ and $\y$ that maximise the bilinear form, subject to conditions $\x'\x = \alpha_x^2$ and $\y'\y = \alpha_y^2$.
For any given $\x$, the $\y$ that maximises the form is given by $\eqref{ymax}$ and the maximum is given by $\eqref{xAy_max_y}$. Thus, now, we are looking to maximise $\alpha_y \sqrt{\x' \amat \amat' x}$ over $\x$. Equivalently, we need to maximise the quadratic form $\alpha_y^2 \x' \amat \amat'º x$ with constraint $\x'\x = \alpha_x^2$. As shown above, the solution to this is the eigenvector of $\amat\amat'$ with eigenvalue $\lambda_\max$, and, from $\eqref{xAx_max}$ it follows:
$$
\begin{eqnarray}
\max(\x' \amat \y) & = & \alpha_y \sqrt{\max(\x'\amat'\amat\x)} \\
& = & \alpha_y \alpha_x \sqrt{\lambda_\max}
\end{eqnarray}
$$

### Mean of maximum

Easier: calculate the mean over $\x$ of the square of the maximum given by $\eqref{xAy_max_y}$. 
$$
\begin{eqnarray}
<\max_\y(\x'\amat \y)^2> & = &  \alpha_y^2 <\x'\amat\amat'\x> \\
& = & \alpha_y^2 \texttt{Tr}(\amat \amat' <\x \x'>) \\
& = & \alpha_x^2 \alpha_y^2 \frac{\texttt{Tr}(\amat \amat')}{n_x}
\end{eqnarray}
$$
To obtain the last equation, I assumed  $\x$ are uniformly distributed onto the  hypersphere $\x'\x = \alpha_x^2$. Then $<x_i x_j> = \alpha_x^2 \delta_{ij}/n_x$, so that $<\x \x'> = \frac{\alpha_x^2}{n_x} \mathbf{I}$.

## Minimum of sum of two vectors

What's the minimum of:
$$
\newcommand{\X}{\mathbf{X}}
\newcommand{\Y}{\mathbf{Y}}
F(\x,\y) = (\X\x + \Y \y)'(\X\x + \Y \y)-\x'\X'\X\x
$$
 subject to constraints: $\x'\x = \alpha_x^2$, $\y'\y = \alpha_y^2$? Equivalently:
$$
F(\x,\y) = \y'\Y'\Y\y + 2\x'\X'\Y\y \label{fxy_sum_vectors}
$$

Let's fist consider the problem of given $\x = \x_0$ find $\y$ of length $\alpha_y$ that minimises $F(\x,\y)$. The Lagrangian is:
$$
\newcommand{\lx}{\mathbf{l}_x}
L(\x,\y) =\y'\Y'\Y\y + 2\x'\X'\Y\y - \lx'(\x - \x_0)-\lambda_y(\y'\y-\alpha_y^2)
$$

where $\lx' = (\lambda_{x_1}, \lambda_{x_2},…)$. The minimum satisfies:
$$
\begin{eqnarray}
\y &=& 	(\Y'\Y - \lambda_y)^{-1}\Y'\X\x_0 \\ \label{ymin_sum_vectors}
	\y'\y &=& \alpha_y^2
\end{eqnarray}
$$
Note that $\lambda_y$ must not be equal to an eigenvalue of $\Y'\Y$ for the inverse of the factor in parenthesis to exist. Replacing $\eqref{ymin_sum_vectors}$ into $\eqref{fxy_sum_vectors}$ we find:
$$

$$

## using least squares?

What's the minimum of:
$$
\newcommand{\X}{\mathbf{X}}
\newcommand{\Y}{\mathbf{Y}}
F(\x,\y) = (\X\x + \Y \y)'(\X\x + \Y \y)
$$
 given $\x$?
$$
F(\x,\y) = \x'\X'\X\x + \y'\Y'\Y\y + 2\x'\X'\Y\y
$$


$$
\newcommand{\lx}{\mathbf{l}_x}
L(\x,\y) =\y'\Y'\Y\y + 2\x'\X'\Y\y - \lx'(\x - \x_0)-\lambda_y(\y'\y-\alpha_y^2)
$$

The minimum satisfies:
$$
\begin{eqnarray}
\y &=& 	-(\Y'\Y)^{-1}\Y'\X\x \\ 
\end{eqnarray}
$$
which is the least-square solution of the overdetermined system of equations $\Y \y = - \X x$. Replacing the last equation in $F$, we find:
$$
F = \x´\X'[1 - (\Y'\Y)^{-1}\Y']^2\X\x
$$

### mean of min

$$
<F> = \texttt{Tr}(\X'[1 - (\Y'\Y)^{-1}\Y']^2\X<\x\x'>)\\
=\alpha_x^2 \texttt{Tr}(\X'[1 - (\Y'\Y)^{-1}\Y']^2\X)/n_x
$$

