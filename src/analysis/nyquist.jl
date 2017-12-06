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
    w_user = !isempty(w)

    w_final = 0.0

    tol = 0.0

    # If the user did not provide a vector with frequencies, then they will be
    # computed automatically.
    if (!w_user)
        # If not, compute the bode diagram from a frequency two decades lower
        # than the minimum frequency and two decades higher than the maximum
        # frequency.
        poles = pole(sys)

        # Get maximum frequency of the system.
        highest_freq  = maximum(abs.(real(poles[poles .!= 0])))

        w_final   = (highest_freq != 0)  ?  highest_freq*100.0 : 1e+2

        # Compute the tolerance using the DC gain.
        dc_gain = maximum(D + C*pinv(eye(A)-A)*B)
        tol     = dc_gain/50.0

        # Don't let tolerance to be 0.
        (tol == 0) && (tol = 1e-2)

        # Initialize the vector `w` with 100 points between the minimum and
        # maximum frequencies.
        w = collect(linspace(0,w_final,100))
    end

    # Frequency response for each pair input/output.
    real_fr     = Vector{Matrix{Float64}}(0)
    imag_fr     = Vector{Matrix{Float64}}(0)
    mir_imag_fr = Vector{Matrix{Float64}}(0)

    # Compute the Nyquist plot.
    i  = 1
    wi = w[1]

    while i <= length(w)
        wi = w[i]

        s = Complex(0,wi)

        fr = C*pinv(s*eye_s-A)*B+D

        real_fr_i  = Float64[real(f) for f in fr]
        imag_fr_i  = Float64[imag(f) for f in fr]

        if w_user || (i == 1)
            i += 1
            push!(real_fr,     +real_fr_i)
            push!(imag_fr,     +imag_fr_i)

            (!w_user) && push!(mir_imag_fr, -imag_fr_i)
        else
            # Verify the maximum distance between the current point and the last
            # one.
            Δx = maximum(abs.(real_fr_i - real_fr[i-1]))
            Δy = maximum(abs.(imag_fr_i - imag_fr[i-1]))

            # Accept the point if the tolerance is higher than the distances.
            if (maximum([Δx;Δy]) > tol) && ( (w[i]-w[i-1]) > 1e-2 )
                if i == 2
                    w = [collect(linspace(w[1],w[2],10)); w[3:end]]
                elseif i == length(w)
                    w = [w[1:i-2]; collect(linspace(w[i-1],w[i],10))]
                else
                    w = [w[1:i-2]; collect(linspace(w[i-1],w[i],10)); w[i:end]]
                end
            else
                i += 1
                push!(real_fr,     +real_fr_i)
                push!(imag_fr,     +imag_fr_i)
                push!(mir_imag_fr, -imag_fr_i)
            end
        end
    end

    if !w_user
        real_fr = [reverse(real_fr);     real_fr]
        imag_fr = [reverse(mir_imag_fr); imag_fr]
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

