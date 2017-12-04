"""
### function plot_bode(gain, phase, w, num_u, num_y)

Plot the bode diagram.

##### Args

* gain: Gain [dB].
* phase: Phase [deg].
* w: Frequencies [rad/s].
* num_u: Number of inputs.
* num_y: Number of outputs.

"""

function plot_bode(gain, phase, w, num_u, num_y)
    if (plot_system == :pyplot)
        # Get the number of graphics.
        num_plots = length(gain[1])

        fig, ax = plt.subplots(num_y*2, num_u, sharex=true, sharey=true)

        for i = 1:num_plots
            gain_i  = map(x->x[i], gain)
            phase_i = map(x->x[i], phase)

            ax[2(i-1)+1][:semilogx](w, gain_i,  linewidth=plw)
            ax[2(i-1)+2][:semilogx](w, phase_i,  linewidth=plw)
            ax[2(i-1)+1][:grid]("on", which="major",
                         color="#AAAAAA", linewidth=glw)
            ax[2(i-1)+2][:grid]("on", which="major",
                         color="#AAAAAA", linewidth=glw)

            # Change figure border width.
            ax[2(i-1)+1][:spines]["top"][:set_linewidth](blw)
            ax[2(i-1)+1][:spines]["bottom"][:set_linewidth](blw)
            ax[2(i-1)+1][:spines]["left"][:set_linewidth](blw)
            ax[2(i-1)+1][:spines]["right"][:set_linewidth](blw)

            ax[2(i-1)+2][:spines]["top"][:set_linewidth](blw)
            ax[2(i-1)+2][:spines]["bottom"][:set_linewidth](blw)
            ax[2(i-1)+2][:spines]["left"][:set_linewidth](blw)
            ax[2(i-1)+2][:spines]["right"][:set_linewidth](blw)
        end

        for i = 1:num_u
            (num_u > 1) && ax[2*num_y*(i-1)+1][:set_title]("From Input $i")
            ax[2*num_y*i][:set_xlabel]("Freq. (rad/s)")
        end

        for i = 1:num_y
            ax[2(i-1)+1][:set_ylabel]("Gain (dB)")
            ax[2(i-1)+2][:set_ylabel]("Phase (deg)")

            ax_aux = ax[2(num_u-1)*num_y+2(i-1)+1][:twinx]()
            ax_aux[:yaxis][:set_label_position]("right")
            ax_aux[:set_yticks]([])
            ax_aux[:set_ylabel]("To Output $i")

            ax_aux = ax[2(num_u-1)*num_y+2(i-1)+2][:twinx]()
            ax_aux[:yaxis][:set_label_position]("right")
            ax_aux[:set_yticks]([])
            ax_aux[:set_ylabel]("To Output $i")
        end

        fig[:suptitle]("Bode Diagram", fontsize=tifsz)
    end
end
