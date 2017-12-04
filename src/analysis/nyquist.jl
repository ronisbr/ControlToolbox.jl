export nyquist

#==#
#
# @brief Compute nyquist diagram.
# @param[in] sys System in state space representation.
# @param[in] w Array of frequencies in which the bode diagram will be computed.
# If empty, the frequencies will be chosen according to the system poles.
# @param[in] plot Plot bode diagram (default = true).
#
# @return Real(G(jw)), Imag(G(jw)), w (freq. array) for each pair input/output.
#
#==#

"""
### function nyquist(sys::StateSpace, w::Array=[]; plot=true)

Compute the Nyquist diagram.

##### Args

* sys: System in state space representation.
* w: Array of frequencies in which the bode diagram will be computed.  If empty,
  the frequencies will be chosen according to the system poles.


##### Returns

Real(G(jw)), Imag(G(jw)), w (freq. array) for each input/output pair.

"""

function nyquist(sys::StateSpace, w::Array=[]; plot=true)
    A, B, C, D = sys.A, sys.B, sys.C, sys.D

    # Number of inputs.
    num_u = size(B,2)

    # Number of outputs.
    num_y = size(C,1)

    # Number of states.
    num_states = size(A,1)

    # Identity matrix.
    eye_s = eye(num_states)

    # Check if the vector of frequencies was provided.
    if (isempty(w))
        # If not, compute the bode diagram from a frequency two decades lower
        # than the minimum frequency and two decades higher than the maximum
        # frequency.
        poles = pole(sys)

        # Get maximum frequency of the system.
        highest_freq  = maximum(abs.(real(poles[poles .!= 0])))

        w_final   = (highest_freq != 0)  ?  highest_freq*100.0 : 1e+2

        w = collect(-w_final:0.01:w_final)
    end

    # Frequency response for each pair input/output.
    real_fr = Vector{Matrix{Float64}}(length(w))
    imag_fr = Vector{Matrix{Float64}}(length(w))

    for i = 1:length(w)
        s = Complex(0,w[i])

        fr = C*inv(s*eye_s-A)*B+D

        real_fr[i]  = Float64[real(f) for f in fr]
        imag_fr[i]  = Float64[imag(f) for f in fr]
    end

    ( plot ) && plot_nyquist(real_fr, imag_fr, num_u, num_y)

    real_fr, imag_fr, w
end

#==#
#
# @brief Compute bode diagram.
# @param[in] G Transfer function of the dynamic system.
# @param[in] w Array of frequencies in which the bode diagram will be computed.
# If empty, the frequencies will be chosen according to the system poles.
# @param[in] plot Plot bode diagram (default = true).
#
# @return Real(G(jw)), Imag(G(jw)), w (freq. array) for each pair input/output.
#
#==#

function nyquist(G::TransferFunction, w::Array=[]; plot=true)
    # Convert the transfer function to state space and call bode.jl.
    nyquist(tf2ss(G), w; plot=plot)
end

