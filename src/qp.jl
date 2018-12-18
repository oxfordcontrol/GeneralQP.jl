mutable struct Data{T}
    """
    Data structure for the solution of
        min         ½x'Px + q'x
        subject to  Ax ≤ b
    with an active set algorithm.

    The most important element is F, which holds
    - F.QR:      an updatable QR factorization of the working constraints
    - F.Z:       a view on an orthonormal matrix spanning the nullspace of the working constraints
    - F.P:       the hessian of the problem
    - F.U & F.D: the ldl factorization of the projected hessian, i.e. F.U'*F.D*F.D = F.Z'*F.P*F.Z

    The Data structure also keeps matrices of the constraints not in the working set (A_ignored and b_ignored)
    stored in a continuous manner. This is done to increase use of BLAS.
    """
    x::Vector{T}

    n::Int
    m::Int

    F::NullspaceHessianLDL{T}
    q::Vector{T}
    A::Matrix{T}
    b::Vector{T}
    working_set::Vector{Int}
    ignored_set::Vector{Int}
    λ::Vector{T}
    residual::T

    iteration::Int
    done::Bool

    e::Vector{T}

    A_ignored::SubArray{T, 2, Matrix{T}, Tuple{UnitRange{Int}, Base.Slice{Base.OneTo{Int}}}, false}
    b_ignored::SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
    A_shuffled::Matrix{T}
    b_shuffled::Vector{T}

    verbosity::Int
    printing_interval::Int
    r_max::Int

    function Data(P::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T},
        x::Vector{T}; r_max=Inf, verbosity=1) where T

        m, n = size(A)
        working_set = findall((A*x - b)[:] .>= -1e-11)
        ignored_set = setdiff(1:m, working_set)

        F = NullspaceHessianLDL(P, Matrix(view(A, working_set, :)'))
        if F.m == 0
            remove_constraint!(F, 0)
        end
        A_shuffled = zeros(m, n)
        l = length(ignored_set)
        A_shuffled[end-l+1:end, :] .= view(A, ignored_set, :)
        b_shuffled = zeros(m)
        b_shuffled[end-l+1:end] .= view(b, ignored_set)

        e = zeros(n);
        λ = zeros(n); λ .= NaN

        new(x, n, m, F, q, A, b, working_set, ignored_set, λ,
            NaN, 0, false, e,
            view(A_shuffled, m-l+1:m, :),
            view(b_shuffled, m-l+1:m),
            A_shuffled, b_shuffled,
            2, 1, r_max)
    end
end

function solve(P::Array{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T},
    x::Vector{T}; kwargs...) where T

    data = Data(P, q, A, b, kwargs)

    if data.verbosity > 0
        print_header(data)
        print_info(data)
    end

    while !data.done && data.iteration <= Inf
        iterate!(data)

        if data.verbosity > 0
            mod(data.iteration, 10*data.printing_interval) == 0 && print_header(data)
            (mod(data.iteration, data.printing_interval) == 0 || data.done) && print_info(data)
        end
    end
    return x
end

function iterate!(data::Data{T}) where{T}
    direction, stepsize, new_constraints = calculate_step(data)
    data.x .+= stepsize*direction

    if !isempty(new_constraints)
        add_constraint!(data, new_constraints[1])
    end
    if isempty(new_constraints) || data.F.m == 0
        if data.F.artificial_constraints > 0
            remove_constraint!(data.F, 0)
        else
            idx = check_kkt!(data)
            !data.done && remove_constraint!(data, idx)
        end
    end
    data.iteration += 1
end

function calculate_step(data)
    gradient = data.F.P*data.x + data.q
    if data.F.D[end] >= data.F.indefinite_tolerance
        qw = data.F.Z'*(gradient)
        direction = -data.F.Z*reverse(data.F.U\(data.F.D\(data.F.U'\reverse(qw))))
        α_min = 1
    else
        e = view(data.e, 1:data.F.m); e[end] = 1
        direction = data.F.Z*reverse(data.F.U\e)
        direction .*= -sign(dot(direction, gradient))
        e .= 0
        α_min = Inf
    end

    ratios = abs.(data.b_ignored - data.A_ignored*data.x)./(data.A_ignored*direction)
    ratios[data.A_ignored*direction .<= 0] .= Inf

    idx = argmin(ratios)
    α_constraint = ratios[idx]

    if α_constraint <= α_min
        new_constraints = [idx]
    else
        new_constraints = []
    end

    α = min(α_min, α_constraint) 
    α_max = roots(Poly(
        [norm(data.x)^2 - data.r_max^2,
        2*dot(direction, data.x),
        norm(direction)^2]
        ))
    if isreal(α_max)        
        α_max = maximum(α_max)
    else
        α_max = Inf
    end

    return direction, min(α, α_max), new_constraints
end

function check_kkt!(data)
    grad = data.F.P*data.x + data.q
    λ = -data.F.QR.R1\data.F.QR.Q1'*grad
    data.λ[1:length(λ)] .= λ
    data.residual = norm(data.F.Z'*grad)
    # data.residual = norm(grad + data.A[data.working_set, :]'*λ)

    idx = NaN
    if all(λ .>= 0)
        data.done = true
    else
        data.done = false
        idx = argmin(λ)
    end

    return idx
end