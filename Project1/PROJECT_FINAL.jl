using Plots
using LinearAlgebra

a, b = 1, 2 # given in the task
y1, y2 = 0, log(2) # y-axis values, min and max

function F(Y)
    N = length(Y) # number of equations(unknowns) in the system
    h = (b - a) / (N + 1) # formula given in the task, spacing bewtween consecutive points
 
    X = similar(Y) # ensures X to be same type and size as Y
    
    X[1] = (Y[2] - 2Y[1] + y1) / h^2 + ((Y[2] - y1) / (2h))^2 + Y[1] - log(a + h) # 'a + h' used instead of just 'a' = better numerical stability and accuracy, commom approach
 
    for i in 2:N-1 # it is to 'N-1' as n is calculated seperetely
        X[i] = (Y[i+1] - 2Y[i] + Y[i-1]) / h^2 + ((Y[i+1] - Y[i-1]) / (2h))^2 + Y[i] - log(a + i*h) # between boundaries, log(a + i*h), on table
    end
 
    X[N] = (y2 - 2Y[N] + Y[N-1]) / h^2 + ((y2 - Y[N-1]) / (2h))^2 + Y[N] - log(b - h) # 'b - h', same explanation as for 'a - h'

    return X
end

function Jacobian(F, X) # F - function, X - vector 
    n = length(X) # calculates size of vector, also size of Jacobian matrix
    J = zeros(n, n) # currently initialization is with 0's, here Jacobian will be stored
    tol = 1e-6

    # calculating each element of Jacobian
    for j in 1:n # column iteration
        xi = copy(X) # this is done so X is not effected in every itteration, instead copy is created
        xi[j] += tol # approximate derivatives at the perturbed points
        Fi = F(xi) # evaluates function F at current point
        for i in 1:n # row iteration
            J[i, j] = (Fi[i] - F(X)[i]) / tol # computing finite difference approximation, difference between 'i' vector and original point
        end
    end
    return J
end

function LU_decomposition(A)
    n = size(A, 1) # calculation of square matrix
    L = zeros(n, n) # initialization of empty matrix, will store LT part of LU decomposition
    U = copy(A) # UT part of LU decomposition

    for k in 1:n-1 # column matrix, exclude last row as it does not need to be eliminated
        for i in k+1:n # row itterates, below diagonal
            L[i, k] = U[i, k] / U[k, k] # calculating to eliminate entries below diagonal
            for j in k+1:n # iterates over columns right to 'k'--> 'k+1', update enetries for U
                U[i, j] -= L[i, k] * U[k, j] # eliminates the entries below diagonal in column 'k' of 'U'
            end
        end
    end

    for i in 1:n
        L[i, i] = 1.0 # set diagonal entries to 1
    end

    return L, U
end

function LU_solve(L, U, b)
    n = length(b) # lenght of vector b
    y = zeros(n)
    x = zeros(n)

    # forward substitution: Ly = b
    y[1] = b[1]/L[1,1]
    for i in 2:n # itterating over each eq in system
        y[i] = b[i] # from each eq y, to correspond b
        for j in 1:i-1
            y[i] -= L[i, j] * y[j] # updating 'y[i]', subtracting product of correspond element from 'L' and preivously calculated value 'y[j]'
        end
        y[i] = y[i]/L[i,i]
    end

    # backward substitution: Ux = y
    for i in n:-1:1 # in reverse order from last eq
        x[i] = y[i]
        for j in i+1:n # iterates over elements of 'x'
            x[i] -= U[i, j] * x[j] # updaitng 'x[i]', similarly like in forward substitution
        end
        x[i] = x[i] / U[i, i] # dividing by diagonal of 'U' at '[i,i]'
    end

    return x
end

function NewtonMethod(F, J, x0; maxIter=1000, tol=1e-10)
    x = copy(x0)
    n = length(x)
    for iter in 1:maxIter
        Jx = Jacobian(F, x)
        L, U = LU_decomposition(Jx)
        deltaX = LU_solve(L, U, -F(x))
        x += deltaX
        if norm(deltaX) < tol
            return x
        end
    end
    println("MAX number of iterations reached!")
end

N = 100
t = LinRange(a, b, N + 2)

J(x) = Jacobian(F, x)

solution = NewtonMethod(F, J, zeros(N); maxIter=1000, tol=1e-10)
    
yAnalysis = [y1; solution; y2]
plot(t, yAnalysis, label="Numerical Analysis Approximation", ylabel="y(t)", xlabel="t", color=:blue)
scatter!(t, yAnalysis)

savefig("PROJECT_FINAL.png")

println("Solution: ", solution)

