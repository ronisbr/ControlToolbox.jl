export bode

"""
### function bode(sys::StateSpace, w::Array=[]; plot=true)

Compute the bode diagram.

##### Args

* sys: System in state space representation.
* w: Array of frequencies in which the bode diagram will be computed. If empty,
  the frequencies will be chosen according to the system poles.

##### Parameters

* plot: Plot the bode diagram (default = `true`).

##### Returns

Gain (dB), Phase (deg), w (freq. array) for each input/output pair.

"""

function bode(sys::StateSpace, w::Array=[]; plot=true)
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

        # Get the minimum and maximum frequency of the system removing the poles
        # at origin.
        highest_freq  = maximum(abs.(real(poles[poles .!= 0])))
        smallest_freq = minimum(abs.(real(poles[poles .!= 0])))

        w_initial = (smallest_freq != 0) ? smallest_freq/100.0 : 1e-2
        w_final   = (highest_freq != 0)  ?  highest_freq*100.0 : 1e+2

        w = collect(w_initial:0.01:w_final)
    end

    # Frequency response for each pair input/output.
    gain  = Vector{Vector{Float64}}(length(w))
    phase = Vector{Vector{Float64}}(length(w))

    for i = 1:length(w)
        s = Complex(0,w[i])

        # `fr` is always a vector. Hence, we can squeeze one dimension here.
        fr = squeeze(C*inv(s*eye_s-A)*B+D,1)

        gain[i]  = 20*log10.([ norm(f) for f in fr ])
        phase[i] = [ angle(f)*180.0/pi for f in fr ]

        if i > 1
            # Make sure that the phase plot is continuous. This can be done by
            # analysing the difference between the current phase and the latest
            # one and checking if it is higher than +180.0 or lower than -180.0.
            Δ = phase[i] - phase[i-1]
            phase[i] += (Δ .> 180.0) * (-360.0) + (Δ .< -180.0) * (+360.0)
        end
    end

    ( plot ) && plot_bode(gain, phase, w, num_u, num_y)

    gain, phase, w
end

"""
### function bode(G::TransferFunction, w::Array=[]; plot=true)

Compute the bode diagram.

##### Args

* G: Transfer function of the dynamic system.
* w: Array of frequencies in which the bode diagram will be computed. If empty,
  the frequencies will be chosen according to the system poles.

##### Parameters

* plot: Plot the bode diagram (default = `true`).

##### Returns

Gain (dB), Phase (deg), w (freq. array) for each input/output pair.

"""

function bode(G::TransferFunction, w::Array=[]; plot=true)
    # Convert the transfer function to state space and compute the Bode diagram.
    bode(tf2ss(G), w; plot=plot)
end

