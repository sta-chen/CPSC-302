# CPSC-302
* [Direct Method and Iterative Method](#direct-method-and-iterative-method)
* [Jacobi Method](#jacobi-method)
* [Gauss-Seidel Method](#gauss-seidel-method)
# Direct Method and Iterative Method
## When to use iterative method:
* If $A$ is large and sparse, $LU$ decomposition / Gaussian elimination may introduce **fill-in**. If the amount of fill-in is significant, then applying the direct method may become costly.
* Sometimes we only need a *rough approximation* $\hat{x}$.
* Sometimes we have a pretty good idea of an approximate guess for the solution.
* Sometimes we only know the matrix-vector products.
## Iterative Method
### Stationary iteration and relaxation methods
Stationary iteratie method: $\mathbf{x}_{k+1} = \mathbf{x}_k + M^{-1} \mathbf{r}_k$  
It's called *stationary* because $M$ is independent of the iteration counter $k$
### Jacobi Method
Simultaneous relaxation  
We choose $M = D$ yielding the $k^{\mathrm{th}}$ iteration: $\mathbf{x}_{k+1} = \mathbf{x}+k+D^{-1}\mathbf{r}_k$, where $D$ is a diagonal matrix consisting the diagonal elements of $A$  
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
We choose $M = E$ yielding the $k^{\mathrm{th}}$ iteration: $\mathbf{x}_{k+1} = \mathbf{x}+k+E^{-1}\mathbf{r}_k$, where $E$ is an $n\times n$ lower triangular matrix consisting of corresponding elements of $A$ and zero elsewhere  
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
