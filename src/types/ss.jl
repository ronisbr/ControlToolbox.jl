export ss, ss2tf

# State Space Constructor
# =======================

"""
### function StateSpace(A::Matrix, B::Matrix, C::Matrix, D::Matrix)

Initialize a state space with the following form:

x = Ax + Bu
y = Cx + Du

##### Args

* A: Matrix A.
* B: Matrix B.
* C: Matrix C.
* D: Matrix D.

##### Returns

The new state space.

"""

function ss(A::Array, B::Array, C::Array, D::Array)
    StateSpace(A,B,C,D)
end

# Print Functions
# ===============

"""
### function Base.show(io::IO, Sys::StateSpace)

Print space state.

##### Args

* io: IO in which the state space will be printed.
* Sys: Space state.

"""

function Base.show(io::IO, Sys::StateSpace)
    print(io, "\n")
    print(io, "A = ")
    print(io, Sys.A, "\n\n")
    print(io, "B = ")
    print(io, Sys.B, "\n\n")
    print(io, "C = ")
    print(io, Sys.C, "\n\n")
    print(io, "D = ")
    print(io, Sys.D, "\n\n")
end

Base.print(io::IO, Sys::StateSpace) = Base.show(io, Sys)

# Conversion from State Space to Transfer Function
# ================================================

"""
### function ss2tf(sys::StateSpace)

Convert a state space representation of a system to a transfer function.

##### Args

* sys: State space representation.

##### Returns

The transfer function.

##### Remarks

The output will always be of type `Float64` due to the operations to convert
between the representations.

"""

function ss2tf(sys::StateSpace)
    ( !issiso(sys) ) && error("ss2tf() is only supported for SISO systems.")

    # Get the invariant zeros of the system.
    z = izero(sys)

    # Get the poles of the system.
    p = pole(sys)

    # Get the gain of the system.
    K = gain(sys)

    # Compute the transfer function.
    zpk(z, p, K)
end
