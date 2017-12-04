export rlocus

"""
### function rlocus(G::TransferFunction, k = []; plot = true)

Compute and, optionally, plot the root locus of the open-loop system represented
by the transfer function `G`.

##### Args

* G: Transfer function.
* k: Vector with the open-loop gains to compute the root locus.

##### Parameters

* plot: If true, then the root locus will be plotted.

##### Returns

* root_locus: An array with the closed-loop poles of the system.
* root_locus_k: An array with the gains used to compute the closed-loop gains.

"""

function rlocus(G::TransferFunction, k = []; plot = true)
    # Check the transfer function.
    if length(G.num.a) > length(G.den.a)
        error("rlocus() does not support improper systems.")
    end

    # Compute the poles and zeros of the transfer function.
    p = convert(Array{Complex128,1}, pole(G))
    sort!(p, by=real)

    z = convert(Array{Complex128,1}, izero(G))

    # Minimal realization of the transfer function.
    Gm = minreal(G)
    mp = pole(Gm)
    mz = izero(Gm)
    sort!(mp, by=real)

    # Array to store the root locus.
    root_locus   = Array{Array{Complex128,1}}(0)
    root_locus_k = Array{Float64,1}(0)

    # Select gains to compute the root locus.
    # =======================================
    #
    # The idea of some parts of this algorithm was obtained from the file
    # `rlocus.m` from Octave-Forge.

    k_user = true

    if isempty(k)
        # Indicate that the use did not provided a gain vector.
        k_user = false

        # In the following, there is an ad-hoc algorithm to estimate the highest
        # gain necessary to provide the entire information to the user.

        # Real-axis break-in and breakaway points.
        # ========================================

        # Compute the real-axis break-in and breakaway points.
        num_dKds = polyder(Gm.den)*Gm.num - Gm.den*polyder(Gm.num)

        # Compute the points.
        break_points = roots(num_dKds)

        # Compute the gains related to those points.
        k_crossing =
            real(map(s->polyval(-G.den,s)/polyval(G.num, s), break_points))

        # Get the maximum gain. If it is empty, then use 100.
        k_max = maximum( [k_crossing; 100.0] )

        # Asymptotes
        # ==========

        k_asym = 0.0

        # Number of asymptotes.
        num_asymptotes = length(mp) - length(mz)

        if num_asymptotes > 0
            # Angle and start point for the asymptotes.
            σ_a = (sum(mp) - sum(mz))/num_asymptotes
            θ_a = [ (2*m+1)*pi/num_asymptotes for m in 0:(num_asymptotes-1) ]

            # Get 2x the maximum distance from the center of the asymptotes to
            # the open-loop poles and zeros and the break points.
            aux = [mp; mz; break_points]

            Δ_max  = 0.0

            if !isempty(aux)
                Δ_max = 2*maximum(abs.( [mp; mz; break_points] - σ_a ) )
            end

            # Get the gain related to Δ_max for each asymptotes.
            s = σ_a + Δ_max*exp.(complex.(0.0, θ_a))
            k_asym = real(map(s->polyval(-G.den,s)/polyval(G.num, s), s))
        end

        # Maximum gain
        # ============

        # Update the maximum gain given the previous information.
        k_max = maximum( [k_max; k_asym] )

        # Compute 50 points between 0 and k_max.
        k = collect(linspace(0, k_max, 50))
    end

    # Compute the tolerance for smoothing.
    # ====================================

    # Compute the tolerance for smoothing the root locus.
    ϵ = 0.0

    if !isempty(mp)
        η = [abs.( maximum(real(mp)) - minimum(real(mp)) );
             abs.( maximum(imag(mp)) - minimum(imag(mp)) ); ]

        # Julia function `minimum` fails if the argument is empty.
        (!isempty(η)) && (ϵ = 0.01*minimum(η))
    end

    # The tolerance must not be 0.
    (ϵ == 0.0) && (ϵ = 1e-2)

    # Compute the root locus.
    # =======================

    # Get the poles of the minimal realization of the system.
    push!(root_locus,   mp)
    push!(root_locus_k, 0.0)

    # Compute closed-loop poles.
    pole_k_1 = mp

    i = 1
    ki = k[i]

    # Notice that, if the minimum realization of the system has no poles, then
    # nothing must be computed.
    while !isempty(mp) && (i <= length(k))
        ki = k[i]

        (ki == 0) && (i += 1; continue)

        # Closed-loop dynamic matrix.
        den_cl = Gm.den + ki*Gm.num

        pole_k = convert(Array{Complex128,1},roots(den_cl))

        # Sort the poles for ki to be as closest as possible to the poles
        # computed in the last iteration.
        for j = 1:length(pole_k_1)
            # Do not consider the past values, since they are already sorted.
            ind = indmin(abs.(pole_k[j:end]-pole_k_1[j]))

            aux             = pole_k[j]
            pole_k[j]       = pole_k[j+ind-1]
            pole_k[j+ind-1] = aux
        end

        # If the user did not provide a gain vector, then add more points to
        # smooth the root locus.
        if !k_user
            # Compute the maximum distance between to current poles and the
            # poles of the last iteration.
            Δ = maximum(norm.(pole_k[1:length(pole_k_1)] - pole_k_1))

            # If the difference if too high, add more gains.
            if (Δ > ϵ) && (k[i-1] < k[i])
                k_new = collect(linspace(k[i-1], k[i], 10))

                if i < length(k)
                    k = [k[1:i-2]; k_new; k[i+1:end]]
                else
                    k = [k[1:i-2]; k_new;]
                end

                continue
            end
        end

        # Poles of the closed-loop system. Notice that we should only select the
        # number of poles in `pole_k_1`. The other should be those that will be
        # cancelled with the zeros.
        pole_k_1 = copy(pole_k[1:length(pole_k_1)])
        push!(root_locus, pole_k_1)
        push!(root_locus_k, ki)
        i += 1
    end

    # Now, put the zeros into the root locus by comparing the distance between
    # them and the latest computed closed-loop poles.
    aux      = copy(pole_k_1)
    zeros_k  = NaN*ones(pole_k_1)

    for i = 1:length(mz)
        ind = indmin(abs.(mz[i]-aux))
        zeros_k[ind] = mz[i]

        # Remove the matched pole to avoid problems with zeros that have
        # multiplicity higher than 1.
        aux[ind] = Inf
    end

    push!(root_locus, zeros_k)
    push!(root_locus_k, Inf)

    # Plot the root locus.
    # ====================

    # Check if the root locus must be plotted.
    if plot
        mz = izero(Gm)
        plot_rlocus(root_locus, root_locus_k, p, z, mp, mz)
    end

    root_locus, root_locus_k
end
