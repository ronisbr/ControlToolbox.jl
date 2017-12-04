export gain, izero, pole

"""
### function gain(sys::StateSpace)

Compute the gain of a system in a State Space representation.

##### Args

* sys: State space.

##### Remarks

The idea was obtained from the routine TG04BX of Octave-Force project.

"""

function gain(sys::StateSpace)
    ( !issiso(sys) ) && error("gain() is only defined for SISO systems.")

    # To compute the gain, notice that:
    #
    #                                #zeros
    #                                 prod  (s-zi)
    #                  -1             i = 1
    #    G(s) = C.(sI-A) . B+D = K . --------------
    #                                #poles
    #                                 prod  (s-pi)
    #                                 i = 1
    #
    # Let s0 \in |R be a number that is neither a pole nor a zeros. Hence,
    #
    # #zeros                        #poles
    #  prod  (s0-zi) |= 0    and     prod  (s0-pi) != 0
    #  i = 1                         i = 1
    #
    # Finally, the system gain can be computed using:
    #
    #                        #poles
    #      _              _   prod  (s0-zi)
    #     |         -1     |  i = 1
    # K = |C.(s0.I-A) . B+D|.--------------
    #     |_              _| #zeros
    #                         prod  (s0-pi)
    #                         i = 1
    #

    # System poles.
    p = pole(sys)

    # System invariant zeros.
    z = izero(sys)

    # Find a real value for s0 that is neither a pole nor a zero.
    #
    # TODO: Improve this algorithm.
    #
    # This version was adapted from routine TG04BX of Octave-Forge project.
    s0 = 2*maximum(abs.([z;p]))

    # Compute:
    #
    #    #zeros
    #     prod  (s0-zi)
    #     i = 1
    #
    # Notice that this is a real number because s0 is real and if zi is complex
    # then its conjugate is also a zero, for all i.
    prod_zeros = prod(s0-z)

    # Compute:
    #
    #    #poles
    #     prod  (s0-pi)
    #     i = 1
    #
    # Notice that this is a real number because s0 is real and if pi is complex
    # then its conjugate is also a pole, for all i.
    prod_poles = prod(s0-p)

    # Compute:
    #              _              _
    #             |         -1     |
    #     G(s0) = |C.(s0.I-A) . B+D|
    #             |_              _|
    #
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    G_s0 = C*inv(s0*eye(A)-A)*B+D

    # Compute the gain.
    (G_s0*prod_poles/prod_zeros)[1]
end

"""
### function izero(sys::StateSpace)

Compute the invariant zeros of the system.

##### Args

* sys: Dynamic system in space state representation.

##### Returns

* The invariant zeros of the system.

"""

function izero(sys::StateSpace)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D

    # Number of states.
    num_states = size(A,1)

    # Compute the invariant zeros using generalized eigenvalues.
    M = [A -B; C -D]

    N = zeros(M)
    N[1:num_states, 1:num_states] = eye(A)

    D,V = eig(M,N)

    # Remove the eigenvalues at +Inf and -Inf.
    D[ -Inf .< real(D) .< +Inf ]
end

"""
### function izero(G::TransferFunction)

Compute the invariant zeros of the system.

##### Args

* G: Transfer function of the system.

##### Returns

* The zeros of the transfer function.

"""

function izero(G::TransferFunction)
    # Compute the roots of the denominator.
    roots(G.num)
end

"""
### function pole(sys::StateSpace)

Compute the poles of the system.

##### Args

* sys: Dynamic system in state space representation.

##### Returns

* The poles of the system.

"""

function pole(sys::StateSpace)
    # Compute the eigenvalues of matrix A.
    eig(sys.A)[1]
end

"""
### function pole(G::TransferFunction)

Compute the poles of the system.

##### Args

* G: Transfer function of the system.

##### Returns

* The poles of the transfer function.

"""

function pole(G::TransferFunction)
    # Compute the roots of the denominator.
    roots(G.den)
end
