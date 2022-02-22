"""
    fdweights(z, x, m)

Calculates FD (finite difference) weights.

Return an array of size m+1 by length(x) containing, in successive rows, the
weights for derivatives 0, 1, ..., m.

# Arguments
- `z`: location where approximations are to be accurate
- `x`: vector of x-coordinates (grid points)
- `m`: highest derivative for which to find weights

# Examples
```julia-repl
julia> fdweights(0, -2:2, 6)
5×5 adjoint(::Matrix{Float64}) with eltype Float64:
 -0.0         0.0        1.0   0.0       -0.0
  0.0833333  -0.666667   0.0   0.666667  -0.0833333
 -0.0833333   1.33333   -2.5   1.33333   -0.0833333
 -0.5         1.0        0.0  -1.0        0.5
  1.0        -4.0        6.0  -4.0        1.0
```
"""
function fdweights(z, x, m)
    n = length(x)
    if m >= n
        m = n - 1 # set to highest possible derivative
    end

    c1 = 1
    c4 = x[1] - z
    C = zeros(n, m+1)
    C[1, 1] = 1
    
    for i = 1:n-1
        i1 = i + 1
        mn = min(i, m)
        c2 = 1
        c5 = c4
        c4 = x[i1] - z
        for j = 0:i-1
            j1 = j + 1
            c3 = x[i1] - x[j1]
            c2 = c2 * c3
            if j == i-1
                for s = mn:-1:1
                    s1 = s + 1
                    C[i1, s1] = c1 * (s * C[i1-1, s1-1] - c5 * C[i1-1, s1]) / c2
                end
                C[i1, 1] = -c1 * c5 * C[i1-1, 1] / c2
            end
        for s = mn:-1:1
            s1 = s + 1
            C[j1, s1] = (c4 * C[j1, s1] - s*C[j1, s1-1]) / c3
        end
        C[j1, 1] = c4 * C[j1,1] / c3
        end
        c1 = c2
    end
    return C'
end
