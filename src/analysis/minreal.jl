export minreal

"""
### function minreal(G::TransferFunction)

Compute the minimal realization of a dynamic system represented by a transfer
function.  In other words, compute a transfer function without common poles and
zeros.

##### Args

* G: Transfer function.

##### Returns

The transfer function of the minimal realization of the dynamic system.

"""

function minreal(G::TransferFunction)
    # Get the list of poles and zeros.
    p = pole(G)
    z = izero(G)

    # Get the system gain.
    k = G.num[end]/G.den[end]

    # Get the poles and zeros of the minimal system.
    mp = p
    mz = Array{Complex128}(0)

    # Auxiliary variable to compute tolarance.
    #
    # This algorithm was adapted from file `__minreal__.m` of Octave Forge.
    sqrt_eps = sqrt(eps())

    # Remove the poles that can be cancelled with zeros.
    for i = 1:length(z)
        # Compute the tolerance.
        tol = (abs(z[i]) < sqrt_eps) ? 10000*sqrt_eps : 10000*abs(z[i])*sqrt_eps

        # Compute the distance from the poles to the zero.
        Δ = abs.(p-z[i])

        # Find if there are any poles within a distance of `tol` from the
        # analyzed zero.
        ind = find(Δ .< tol)

        if !isempty(ind)
            # Delete the closest pole to the current zero.
            deleteat!(mp, indmin(Δ))
        else
            push!(mz,round(z[i],6))
        end
    end

    # Return the minimal transfer function.
    zpk(mz, mp, k)
end
