export margin

function margin(G::TransferFunction; plot = true, tol = sqrt(eps()))
    # Gain Margin
    # ===========

    # Compute the phase crossover frequency, which is used to compute the gain
    # margin.
    jw  = Poly([0;im])
    p   = G.num(jw)*conj(G.den(jw))
    imp = Poly(imag(p.a))
    wv  = roots(imp)

    # Search for the real, positive frequency that yields the lowest gain
    # margin.
    wcg = NaN
    Gm  = Inf

    for w in wv
        # Discard imaginary solutions.
        (abs(imag(w)) > tol) && continue

        rw = real(w)

        # Discard negative solutions.
        (rw < 0.0) && continue

        Gw = G.num(rw*im)/G.den(rw*im)

        # The real part of the transfer function at `w` must be negative.
        ( real(Gw) > 0.0 ) && continue

        # Do not consider undefined solutions.
        (isnan(Gw)) && continue

        # Compute the new gain margin candidate.
        Gm_candidate = -real(G.den(rw*im)/G.num(rw*im))

        # Check if it is lower than the last one.
        (Gm_candidate < Gm) && (Gm = Gm_candidate; wcg = rw)
    end

    # Phase Margin
    # ============

    # Compute the gain crossover frequency, which is used to compute the phase
    # margin.
    p = polyval(G.num,jw)*conj(polyval(G.num,jw)) -
        polyval(G.den,jw)*conj(polyval(G.den,jw))

    wv = roots(p)

    # Search for the real, positive frequency that yields the lowest phase
    # margin.
    wcp  = NaN
    Pm   = Inf

    for w in wv
        # Discard imaginary solutions.
        (abs(imag(w)) > tol) && continue

        rw = real(w)

        # Discard negative solutions.
        (rw < 0.0) && continue

        # Compute the new phase margin candidate.
        Pm_candidate = 180.0 + angle(G.num(rw*im)/G.den(rw*im))*180.0/pi

        # If it is greater than 180 deg, then we must subtract 360 deg.
        (Pm_candidate > 180.0) && (Pm_candidate -= 360.0)

        # Check if it is lower than the last one.
        (Pm_candidate < Pm) && (Pm = Pm_candidate; wcp = rw)
    end

    if plot
        # Compute bode diagram and plot.
        gain, phase, w = bode(G; plot = false)
        plot_margin(gain, phase, w, Gm, Pm, wcg, wcp)
    end

    Gm, Pm, wcg, wcp
end




