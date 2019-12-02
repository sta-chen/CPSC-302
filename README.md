# CPSC-302
* [Linear Algebra Basic](#linear-algebra-basic)
* [Linear Least Squares Problems](#linear-least-squares-problems)
* [Direct Method and Iterative Method](#direct-method-and-iterative-method)
    * [LU Decomposition](#lu-decomposition)
    * [Jacobi Method](#jacobi-method)
    * [Gauss-Seidel Method](#gauss-seidel-method)
# Linear Algebra Basic
1. An $n\times n$ matrix $A$ will have at most $n$ eigenvalues. Solve $\mathrm{det}(A - \lambda I) = 0$ to find eigenvalues. Complex eigenvalues of a real matrix occur in conjugate pairs.
2. $A$ is a square matrix, if $\mathrm{det}(A) = 0$, then $A$ is **singular** matrix. Nonsingular means invertible.
3. $A$ is **symmetric** if and only if $A = A^T$.
    * $A$ is also **positive definite** if and only if all its eigenvalues are positive.
    * $A$ is also **positive definite** if and only if $\mathbf{x}^TA\mathbf{x}>0$ for all $\mathbf{x} \neq \mathbf{0}$
4. **Condition number** of matrix $A$ is denoted as $\kappa (A) = \|A\|\|A^{-1}\|$
    * $\kappa(A) = \|A\|\|A^{-1}\| \ge \|A^{-1}A\| = \|I\| = 1$
    * Matrix $A$ is ideally conditioned when $\kappa(A) = 1$
    * $A$ is symmetric positive definite, $\kappa_2(A) = \frac{\lambda_1}{\lambda_n} = \frac{\mathrm{max}(\lambda)}{\mathrm{min}(\lambda)}$
    * $A$ is an $m\times n$ rectangular matrix, $B = A^TA$ with eigenvalues $\lambda$, $\kappa_2(A) = \frac{\sigma_1}{\sigma_n} = \sqrt{\frac{\lambda_1}{\lambda_n}}$ because $\kappa_2(A^TA) = \frac{\lambda_1}{\lambda_n} = [\kappa(A)]^2$
5. **Relative error estimation** is bounded by $\frac{\|x-\hat{x}\|}{\|x\|} \le \kappa(A)\frac{\|\hat{r}\|}{\|b\|}$
6. Matrix norm of $A$
    * $\|A\|_\infty$ is the maximum of the absolute row sum
    * $\|A\|_1$ is the maximum of the absolute column sum
    * $\|A\|_2 = \sqrt{\rho(A^T A)}$ is the largest singular value of $A$
7. **Spectral radius** of a square matrix is the maximum eigenvalue magnitue: $\rho = \mathrm{max}|\lambda_i|$
# Linear Least Squares Problems
Stardard methods for solving the linear least squares problem include the following:  
1. Normal equations: fast, simple, intuitive, but less robust in ill-conditioned situations.
2. QR decomposition: this is the stardard approach implemented in general-purpose software. It is more computationally expensive than the normal equations approach if $m\gg n$ but is more robust.
3. SVD: used mostly when $A$ is rank deficient or nearly rank deficient because the QR approach may not be sufficiently robust in this case. The SVD-based approach is very robust but is significantly more expensive in general and cannot be adapted to deal efficiently with sparse matrices.
## Gram-Schmidt Orthogonalization
It's numerically unstable if the columns of $A$ is nearly linearly dependent; Using $\mathbf{a}_2$ and $\mathbf{q}_1$ to compute $\mathbf{r}_{12}$: $\mathbf{r}_{ij} = <\mathbf{a}_j, \mathbf{q}_i>$.  
In Gram-Schmidt orthogonalization, it's assumed the columns of $A$ are linearly independent; if $\mathbf{a}_j$ is linearly dependent on $\mathbf{a}_1$ through $\mathbf{a}_{j-1}$, then at $j^{\mathrm{th}}$ stage of the algorithm we will get $\mathbf{r}_{jj} = \mathbf{0}$.  
**Modified Gram-Schmit Orthogonalization**: using $\mathbf{q}_1$ through $\mathbf{q}_{j-1}$ to construct $\mathbf{q}_j$. Because $\mathbf{q}_i$ are orthogonal to one another and less prone to damaging effects of roundoff errors. $\mathbf{r}_{ij} = <\mathbf{q}_j, \mathbf{q}_i>$.
# Direct Method and Iterative Method
## When to use iterative method:
* If $A$ is large and sparse, $LU$ decomposition / Gaussian elimination may introduce **fill-in**. If the amount of fill-in is significant, then applying the direct method may become costly.
* Sometimes we only need a *rough approximation* $\hat{x}$.
* Sometimes we have a pretty good idea of an approximate guess for the solution.
* Sometimes we only know the matrix-vector products.
## Direct Method
### LU Decomposition
Given a real nonsingular matrix $A$, apply LU decomposition first, $A = LU$;  
Given a right-hand-side vector $\mathbf{b}$  
1. Forward substitution: solve $L\mathbf{y} = \mathbf{b}$
2. Backward substitution: solve $U\mathbf{x} = \mathbf{y}$    

Cost of LU decomposition: $\mathcal{O}(n^3)$  
$\mathrm{det}(A) = \mathrm{det}(L)\cdot \mathrm{det}(U)$
### Cholesky Decomposition
Gaussian elimination without pivoting is a stable algorithm for symmetric positive definite matrics.  
Since $A$ is symmetric positive definite, we can stably write $A = LU$
![](https://i.imgur.com/TD3zYWi.png)
$D$ is a square matrix with positive element $u_{kk}$ on its diagonal.  
$A = LU = A^T = LDL^T$, $D^{\frac{1}{2}} = \mathrm{diag}\{\sqrt{u_{11}}, ..., \sqrt{u_{nn}}\}$, $A = GG^T$ where $G = LD^{\frac{1}{2}}$, and $G$ is a lower triangular matrix.  
Compute directly, using symmetry, in $\frac{1}{3}n^3+\mathcal{O}(n^2)$ flops.
```matlab
R = chol(A);
```
In MATLAB, `chol(A)` gives $R = G^T$  
The cost of Cholesky decomposition half the space and half the flops of LU decomposition.
## Iterative Method
### Stationary iteration and relaxation methods
Stationary iteratie method: $\mathbf{x}_{k+1} = \mathbf{x}_k + M^{-1} \mathbf{r}_k$  
It's called *stationary* because $M$ is independent of the iteration counter $k$
### Jacobi Method
Simultaneous relaxation  
We choose $M = D$ yielding the $k^{\mathrm{th}}$ iteration: $\mathbf{x}_{k+1} = \mathbf{x}_k+D^{-1}\mathbf{r}_k$, where $D$ is a diagonal matrix consisting the diagonal elements of $A$  
**e.g.** $A = \pmatrix{1&4&7 \cr 2&5&8\cr 3&6&9}$, then $D = \pmatrix{1&0&0\cr 0&5&0\cr 0&0&9}$  
``` matlab
Dv = diag(A);
[~,n] = size(A);
x = zeros(n,1);
r = b;
tol = 1.e-6;
for i = 1:100
    x = x + r./Dv;
    r = b - A*x;
    if norm(r)/norm(b) < tol, break, end
end
A*x, b
```
### Gauss-Seidel Method
We choose $M = E$ yielding the $k^{\mathrm{th}}$ iteration: $\mathbf{x}_{k+1} = \mathbf{x}_k+E^{-1}\mathbf{r}_k$, where $E$ is an $n\times n$ lower triangular matrix consisting of corresponding elements of $A$ and zero elsewhere  
**e.g.** $A = \pmatrix{1&4&7 \cr 2&5&8\cr 3&6&9}$, then $E = \pmatrix{1&0&0\cr 2&5&0\cr 3&6&9}$  
```matlab
E = tril(A);
[~,n] = size(A);
x = zeros(n,1);
r = b;
tol = 1.e-6;
for i = 1:100
x = x + E \ r;
r = b - A*x;
if norm(r)/norm(b) < tol, break, end
end
A*x, b
```
