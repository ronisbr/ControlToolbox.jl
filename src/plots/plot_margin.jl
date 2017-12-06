"""
### function plot_margin(gain, phase, w, Gm, Pm, wcg, wcp)

Plot the bode diagram with the margin information.

##### Args

* gain: Gain [dB].
* phase: Phase [deg].
* w: Frequencies [rad/s].
* Gm: Gain margin.
* Pm: Phase margin [deg/s].
* wcg: Phase cross frequency for the gain margin computation [rad/s].
* wcp: Gain cross frequency for the phase margin computation [rad/s].

##### Parameters

* tol: Tolerance to check if a value is grater than 0.

"""

function plot_margin(gain, phase, w, Gm, Pm, wcg, wcp; tol = sqrt(eps()))
    # Compute the gain margin in dB.
    Gm_dB = 20*log10(Gm)

    if (plot_system == :pyplot)
        fig, ax = plt.subplots(2, 1, sharex=true)

        # Set window title.
        fig[:canvas][:set_window_title]("Bode Diagram")

        # Plot data.
        ax[1][:semilogx](w, gain,  linewidth=plw, zorder=6)
        ax[2][:semilogx](w, phase, linewidth=plw, zorder=6)
        ax[1][:grid]("on", which="major", color="#AAAAAA", linewidth=glw)
        ax[1][:grid](zorder = 0)
        ax[2][:grid]("on", which="major", color="#AAAAAA", linewidth=glw)
        ax[2][:grid](zorder = 0)

        # Change figure border width.
        ax[1][:spines]["top"][:set_linewidth](blw)
        ax[1][:spines]["bottom"][:set_linewidth](blw)
        ax[1][:spines]["left"][:set_linewidth](blw)
        ax[1][:spines]["right"][:set_linewidth](blw)

        ax[2][:spines]["top"][:set_linewidth](blw)
        ax[2][:spines]["bottom"][:set_linewidth](blw)
        ax[2][:spines]["left"][:set_linewidth](blw)
        ax[2][:spines]["right"][:set_linewidth](blw)

        # Title of the figure.
        title = s"$\bf{GM} = $"

        if !isinf(Gm)
            title *= "$(round(Gm_dB,2)) dB at $(round(wcg,2)) rad/s, "
        else
            title *= "Inf, "
        end

        title *= s"$\bf{PM} = $"

        if !isinf(Pm)
            title *= "$(round(Pm,2)) deg at $(round(wcp,2)) rad/s"
        else
            title *= "Inf"
        end

        fig[:suptitle](title, fontsize=tifsz)

        # Axes labels.
        ax[2][:set_xlabel]("Freq. (rad/s)")

        ax[1][:set_ylabel]("Gain (dB)")
        ax[2][:set_ylabel]("Phase (deg)")

        # Plot a line indicating the margins if applicable.
        if (!isinf(Gm)) && (wcg > tol)
            ax[1][:arrow](wcg, 0, 0, -Gm_dB,
                          linewidth=3.0, length_includes_head=true,
                          head_width=0, head_length=0,
                          facecolor="#F0B300", edgecolor="#F0B300", zorder=4)

            ax[1][:axhline](y = 0.0, color="k", zorder=5)
            ax[2][:axvline](x = wcg, color="k", zorder=3, linestyle="--")
            ax[1][:axvline](x = wcg, color="k", zorder=3, linestyle="--")
        end

        if !isinf(Pm) && (wcp > tol)
            # Find the approximate phase value of the transfer function at the
            # crossing frequency.
            ind_wcp   = indmin(abs.(w-wcp))
            phase_wcp = phase[ind_wcp][1]

            # Find the reference that was used to compute the phase margin.
            pc_ref = round(phase_wcp/180)*180
            (pc_ref == 0.0) && (pc_ref = sign(phase_wcp)*180.0)

            # Avoid problems if phase margin is 0.
            if (Pm != 0.0)
                ax[2][:arrow](wcp, pc_ref, 0, +Pm,
                              linewidth=3.0, length_includes_head=true,
                              head_width=0, head_length=0,
                              facecolor="#F0B300", edgecolor="#F0B300", zorder=3)
            end

            ax[2][:axhline](y = pc_ref, color="k", zorder=5)
            ax[2][:axvline](x = wcp,    color="k", zorder=3, linestyle="--")
            ax[1][:axvline](x = wcp,    color="k", zorder=3, linestyle="--")
        end
    end

    nothing
end
