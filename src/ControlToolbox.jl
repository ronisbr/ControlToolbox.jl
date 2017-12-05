VERSION >= v"0.4.0-dev+6521" && __precompile__()

module ControlToolbox

using PyCall
using Polynomials
using OrdinaryDiffEq

import Base: +, -, *, /
import Base: display, print, step, show

export StateSpace
export TransferFunction
export change_lsim_tolerances

################################################################################
#                                    Types
################################################################################

#                                 State-Space
# ==============================================================================

"""
Space state definition.
"""
immutable StateSpace
    A::Array
    B::Array
    C::Array
    D::Array

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

    function StateSpace(A::Array, B::Array, C::Array, D::Array)

        # Check if the dimensions are consistent.
        #
        # dx/dt = Ax + Bu
        #     y = Cx + Du
        #

        num_states = size(A,1)
        num_measurements = size(C,1)
        num_inputs = size(B,2)

        ( isempty(B) ) && (B = zeros(num_states); num_inputs = 1)
        ( isempty(D) ) && (D = zeros(num_measurements, num_inputs) )

        ( num_states != size(A,2) ) &&
        error("Matrix A must be square.")

        ( num_states != size(B,1) ) &&
        error("Inconsistent dimensions between matrices A and B.")

        ( num_states != size(C,2) ) &&
        error("Inconsistent dimensions between matrices A and C.")

        ( num_measurements != size(D,1) ) &&
        error("Inconsistent dimensions between matrices C and D.")

        ( (num_inputs != size(D,2)) ) &&
        error("Inconsistent dimensions between matrices B and D.")

        new(A,B,C,D)
    end
end

#                              Transfer Function
# ==============================================================================

"""
Transfer function definition.
"""

immutable TransferFunction
    num::Poly
    den::Poly

    """
    ### function TransferFunction(num::Poly, den::Poly)

    Initialize a new transfer function.

    ##### Args

    * num: Numerator (see Poly).
    * den: Denominator (see Poly).

    ##### Returns

    The new transfer function.

    """

    function TransferFunction(num::Poly, den::Poly)
        ( ( length(num.a) == 0 ) || ( length(den.a) == 0 ) ) &&
        error("Numerator or denominator cannot be empty.")

        if (num.a[end] == 0)

            # Check if the numerator is 0.
            if length(num.a) == 1
                new(Poly(0), den)
            else
                TransferFunction(Poly(num.a[1:end-1]), den)
            end
        elseif (den.a[end] == 0)
            TransferFunction(num, Poly(den.a[1:end-1]))
        else
            new(num, den)
        end
    end
end

################################################################################
#                                   Includes
################################################################################

include("aux.jl")

include("analysis/bode.jl")
include("analysis/gain_poles_zeros.jl")
include("analysis/margin.jl")
include("analysis/minreal.jl")
include("analysis/nyquist.jl")
include("analysis/pole_info.jl")
include("analysis/rlocus.jl")
include("simulation/impulse.jl")
include("simulation/lsim.jl")
include("simulation/ramp.jl")
include("simulation/step.jl")
include("plots/plot.jl")
include("types/ss.jl")
include("types/tf.jl")

end # module
