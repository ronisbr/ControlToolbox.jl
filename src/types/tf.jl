export feedback, series, parallel, tf, tf2ss, zpk

# Transfer Function Constructors
# ==============================

"""
### function zpk(z::Array, p::Array, K::Number)

Compute the transfer function given the list of poles and zeros and the
gain.

##### Args

* z: Array with the zeros.
* p: Array with the poles.
* K: System gain.

##### Returns

* The transfer function.

"""

function zpk(z::Array, p::Array, K::Number; real_coefs = true)
    # Number of zeros.
    num_z = length(z)

    # Number of poles.
    num_p = length(p)

    # Numerator.
    num = K*poly(z)

    # Denominator.
    den = poly(p)

    # Check if the transfer function must only have real coefficients.
    if real_coefs
        num_r = Poly(floor.(real(num.a),10))
        den_r = Poly(floor.(real(den.a),10))

        num = num_r
        den = den_r
    end

    # Transfer function.
    TransferFunction(num, den)
end


"""
### function tf(num::Vector, den::Vector)

Create a new transfer function given a numerator and denominator.

##### Args

* num: Array with the numerator coefficients. Highest order first.
* den: Array with the denominator coefficients. Highest order first.

##### Returns

The new transfer function with numerator `num` and denominator `den`.

"""

function tf(num::Vector, den::Vector)
    # Poly receives a vector with lowest order first. Hence, we must reverse the
    # order of the array.
    TransferFunction(Poly(reverse(num)), Poly(reverse(den)))
end

tf(num::Vector) = TransferFunction(num, [1])

# Print Functions
# ===============

"""
### function Base.show(io::IO, G::TransferFunction)

Print a transfer function.

##### Args

* io: IO in which the transfer function will be printed.
* G: Transfer function that will be printed.

"""

function Base.show(io::IO, G::TransferFunction)
    const margin::Int64 = 4

    # HACK: Round the numerator and denominator to 4 decimal digits. This is
    # done by Octave. We should find a way to improve the precision in roots
    # computation to avoid such errors when computing the minimum realization of
    # the system.
    num = Poly(round.(G.num.a,4), :s)
    den = Poly(round.(G.den.a,4), :s)

    num_str = "$(num)"
    den_str = "$(den)"

    max_size = max(length(num_str)-6, length(den_str)-6)

    # Compute trailing spaces to center the transfer function.
    trailing_num = ceil(Int64,(max_size-( length(num_str)-6) )/2)+margin
    trailing_den = ceil(Int64,(max_size-( length(den_str)-6) )/2)+margin

    print(io, "\n")
    print(io, " "^trailing_num)
    printpoly(STDOUT, num; descending_powers=true)
    print("\n")
    print(io, " "^margin, "-"^max_size, "\n")
    print(io, " "^trailing_den)
    printpoly(STDOUT, den; descending_powers=true)
    print("\n")
end

Base.print(io::IO, G::TransferFunction) = show(io, G)

# Transfer Functions Operations
# =============================

function +(G1::TransferFunction, G2::TransferFunction)
    # Check if the denominators are equal.
    if (G1.den == G2.den)
        num_new = G1.num + G2.num
        den_new = G1.den
    else
        num_new = G1.num*G2.den + G2.num*G1.den
        den_new = G1.den*G2.den
   end

    TransferFunction(num_new, den_new)
end

+(G1::TransferFunction, k::Number) =
    TransferFunction(G1.num+k*G1.den, G1.den)

+(k::Number, G2::TransferFunction) =
    TransferFunction(G2.num+k*G2.den, G2.den)

-(G1::TransferFunction, G2::TransferFunction) =
    G1+(-1*G2)

-(G1::TransferFunction, k::Number) =
    G1+(-1*k)

-(k::Number, G1::TransferFunction) =
    k+(-1*G1)

*(G1::TransferFunction, G2::TransferFunction) =
    TransferFunction(G1.num*G2.num, G1.den*G2.den)

*(k::Number, G2::TransferFunction) =
    TransferFunction(k*G2.num, G2.den)

*(G1::TransferFunction, k::Number) =
    TransferFunction(G1.num*k, G1.den)

/(G1::TransferFunction, G2::TransferFunction) =
    TransferFunction(G1.num*G2.den, G1.den*G2.num)

/(k::Number, G2::TransferFunction) =
    TransferFunction(k*G2.den, G2.num)

/(G1::TransferFunction, k::Number) =
    TransferFunction(G1.num, G1.den*k)

"""
### function feedback(G::TransferFunction, H::TransferFunction, m::Int64 = -1)

Compute the transfer function of the following system:


                ┌────────┐
         +      │        │
    ─────>O────>│  G(s)  ├───┬────────>
          ^ (m) │        │   │
          │     └────────┘   │       (m) should be +1 or -1 indicating a
          │     ┌────────┐   │       positive or negative feedback. By
          │     │        │   │       default, it is set to -1.
          └─────│  H(s)  │<──┘
                │        │
                └────────┘

or:

        G(s)
    ------------ , m = -1
    1 + G(s)H(s)

        G(s)
    ------------ , m = +1
    1 - G(s)H(s)

Notice that `minreal()` function is called to return the minimal realization of
the system.

##### Args

* G: Transfer function in the direct path.
* H: Transfer function in the feedback path.
* m: Feedback sum sign. It must be `-1` or `+1`. By default, it is set to `-1`.

##### Returns

The minimal realization of the following transfer function:

        G(s)
    ------------ , m = -1
    1 + G(s)H(s)

        G(s)
    ------------ , m = +1
    1 - G(s)H(s)

"""

function feedback(G::TransferFunction, H::TransferFunction, m::Int64 = -1)
    if m == -1
        minreal(G/(1+G*H))
    elseif m == +1
        minreal(G/(1-G*H))
    else
        error("feedback(): `m` must be -1 or +1.")
    end
end

"""
### function series(G1::TransferFunction, G2::TransferFunction)

Compute the transfer function of the following system:


         ┌────────┐    ┌────────┐
         │        │    │        │
    ────>│  G(s)  ├───>│  H(s)  ├───>
         │        │    │        │
         └────────┘    └────────┘
or:

    G(s)H(s)

##### Args

* G1: First function.
* G2: Second function.

##### Returns

The following transfer function:

    G(s)H(s)

"""

function series(G1::TransferFunction, G2::TransferFunction)
    G1*G2
end

"""
### function parallel(G1::TransferFunction, G2::TransferFunction, m::Int64 = +1)

Compute the transfer function of the following system:


                 ┌─────────┐
                 │         │             (m) should be +1 or -1 indicating a
           ┌─────┤  G1(s)  ├─────┐       sum or subtraction. By default, it is
           │     │         │     │       set to +1.
           │     └─────────┘   + │
     ──────┤                     O─────>
           │     ┌─────────┐ (m) ^
           │     │         │     │
           └─────┤  G2(s)  │<────┘
                 │         │
                 └─────────┘

or:

    G(s) + H(s), m = +1

    G(s) - H(s), m = -1

##### Args

* G1: First function.
* G2: Second function.
* m: Sum signal. It must be `+1` or `-1`. By default, it is set to `+1`.

##### Returns

The following transfer function:

    G1(s) + G2(s), m = +1

    G1(s) - G2(s), m = -1

"""

function parallel(G1::TransferFunction, G2::TransferFunction, m::Int64 = +1)
    if m == +1
        G1 + G2
    elseif m == -1
        G1 - G2
    else
        error("parallel(): `m` must be -1 or +1.")
    end
end

# Conversion from Transfer Function to State Space
# ================================================

"""
### function tf2ss(G::TransferFunction, form=:controllable)

Convert a transfer function to state space.

##### Args

* G: The transfer function that will be converted.
* form: Representation form (only `:controllable` is currently implemented).

##### Returns

* The space state representation.

##### Remarks

The output will always be of type `Float64` due to the operations to convert
between the representations.

##### Example

```julia
tf2ss(G)
```

"""

function tf2ss(G::TransferFunction, form=:controllable)
    # Get system dimension.
    sys_dim = length(G.den.a)-1

    # Get the number of elements of numerator.
    num_dim = length(G.num.a)-1

    # Check if G is proper.
    ( sys_dim < num_dim ) && error("The transfer function must be proper.")

    den = G.den.a

    # Get the coeficient of s^sys_dim.
    a0 = den[end]

    # Check if G is strictly proper.
    if ( sys_dim == num_dim )
        # If not, then transform
        #
        #      num_0      num_1
        # G = -------- = ------- + K = G1 + K
        #      den_0      den_1
        #
        # in which K is scalar and G1 is strictly proper.
        #

        a, b = divrem(G.den, G.num)

        # In this case, there is a direct path input-output.
        D = a.a''

        num = -b.a
        num_dim -= 1
    else
        # In this case, there is not a direct path input-output.
        D = zeros(1, 1)

        num = G.num.a
    end

    # Controllable canonical form.
    if ( form == :controllable )
        # Fill the matrices.
        A = [ zeros(sys_dim-1) eye(sys_dim-1);
             [-den[i]/a0 for i = 1:sys_dim ]' ];

        B = zeros(sys_dim)
        B[end] = 1

        C = zeros(1, sys_dim)
        C[1,1:num_dim+1] = [ num[i]/a0 for i = 1:num_dim+1 ]'

    elseif ( form == :observable )
       error("Conversion from transfer function to state space observable form is not implemented yet.")
    else
        error("Invalid form!")
    end

    StateSpace(A,B,C,D)
end

