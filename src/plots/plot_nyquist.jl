"""
### function plot_nyquist(real_fr, imag_fr, num_u, num_y)

Plot the Nyquist diagram.

##### Args

* real_fr: Real part of the frequency response.
* imag_fr: Imaginary part of the frequency response.
* num_u: Number of inputs.
* num_y: Number of outputs.

"""

function plot_nyquist(real_fr, imag_fr, num_u, num_y)
    if (plot_system == :pyplot)
        # Get the number of graphics.
        num_plots = length(real_fr[1])

        fig, ax = plt.subplots(num_y, num_u, sharex=true, sharey=true)

        ( num_plots == 1 ) && (ax = [ax])

        for i = 1:num_plots
            real_fr_i = map(x->x[i], real_fr)
            imag_fr_i = map(x->x[i], imag_fr)

            # Get the default color for the plots.
            color = pybuiltin("next")(ax[i][:_get_lines][:prop_cycler])["color"]

            ax[i][:plot](real_fr_i, +imag_fr_i, color=color, linewidth=plw)
            ax[i][:plot](real_fr_i, -imag_fr_i, color=color, linewidth=plw)

            ax[i][:grid]("on", which="major",
                         color="#AAAAAA", linewidth=glw, zorder=1)

            # Mark the -1 point but do not resize the plot.
            ax[i][:autoscale](false)
            ax[i][:scatter]([-1.0],[0.0],
                            s=mksz, marker="+", linewidths=plw,
                            edgecolor="red", color="red", zorder=3)

            # Create arrows at:
            #     1) 1/4 of the range of the real part;
            #     2) 1/2 of the range of the real part;
            #     3) 3/4 of the range of the real part.
            min_real_fr = minimum(real_fr_i[isfinite.(real_fr_i)])
            max_real_fr = maximum(real_fr_i[isfinite.(real_fr_i)])
            Δp_real     = max_real_fr - min_real_fr

            fq_x  = min_real_fr + 1*Δp_real/4
            mid_x = min_real_fr + 2*Δp_real/4
            tq_x  = min_real_fr + 3*Δp_real/4

            for p in [fq_x, mid_x, tq_x]
                # Search for the freq. that provides the closest real part to p.
                ind = indmin(abs.(real_fr_i-p))

                # Compute the tangent.
                Δ = [real_fr_i[ind+1] - real_fr_i[ind];
                     imag_fr_i[ind+1] - imag_fr_i[ind];]

                ax[i][:annotate]("",
                                 xytext = [real_fr_i[ind];
                                           imag_fr_i[ind]],
                                 xy     = [real_fr_i[ind]+Δ[1]*0.01;
                                           imag_fr_i[ind]+Δ[2]*0.01],
                                 arrowprops = Dict("arrowstyle"=>"simple"),
                                 size       = arrsz)

                ax[i][:annotate]("",
                                 xytext = [ real_fr_i[ind+1];
                                           -imag_fr_i[ind+1]],
                                 xy     = [ real_fr_i[ind+1]-Δ[1]*0.01;
                                           -imag_fr_i[ind+1]+Δ[2]*0.01],
                                 arrowprops = Dict("arrowstyle"=>"simple"),
                                 size       = arrsz)
            end

            # Change figure border width.
            ax[i][:spines]["top"][:set_linewidth](blw)
            ax[i][:spines]["bottom"][:set_linewidth](blw)
            ax[i][:spines]["left"][:set_linewidth](blw)
            ax[i][:spines]["right"][:set_linewidth](blw)
        end

        for i = 1:num_u
            (num_u > 1) && ax[num_y*(i-1)+1][:set_title]("From Input $i")
            ax[num_y*i][:set_xlabel]("Real")
        end

        for i = 1:num_y
            ax[i][:set_ylabel]("Imag.")

            ax_aux = ax[(num_u-1)*num_y+(i-1)+1][:twinx]()
            ax_aux[:yaxis][:set_label_position]("right")
            ax_aux[:set_yticks]([])
            ax_aux[:set_ylabel]("To Output $i")
        end

        fig[:suptitle]("Nyquist Diagram", fontsize=tifsz)
    end
end
