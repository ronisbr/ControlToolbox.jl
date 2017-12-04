# HACK: This global variable stores information about the root locus diagrams
# that are open, so that it is possible to interact with them. This was the only
# way I could think to handle such case. This dictionary is indexed by the
# pointer of the figure canvas. In python, they create a class for each open
# figure and the variables are stored in an instantiation of this class. In
# julia, I do not know any good way to do this.
global rlocus_dict = Dict{Any, Any}()

"""
### function plot_rlocus(root_locus, root_locus_k, p, z, mp, mz)

Plot the root locus.

##### Args

* root_locus: Array with the closed-loop poles.
* root_locus_k: Array with the gains used to compute the closed-loop poles.
* p: Poles of the system.
* z: Zeros of the system.
* mp: Poles of the minimum realization of the system.
* mz: Zeros of the minimum realization of the system.

"""

function plot_rlocus(root_locus, root_locus_k, p, z, mp, mz)
    if plot_system == :pyplot
        # Create the figure.
        # ==================
        fig = plt.figure()
        ax  = plt.gca()

        # Plot axes.
        plt.axhline(0.0, color="black", linewidth=alw, zorder = 1)
        plt.axvline(0.0, color="black", linewidth=alw, zorder = 1)

        # Configure title and labels.
        plt.title("Root Locus")
        plt.xlabel("Real Axis")
        plt.ylabel("Imaginary Axis")

        # Add the root locus information to the global dictionary.
        rlocus_dict[fig[:canvas]] =
            (root_locus, root_locus_k)

        fig[:canvas][:mpl_connect]("button_press_event", rlocus_on_mouse_press)

        # Get real and imaginary parts of root locus.
        rl_real = map(x->real(x), root_locus)
        rl_imag = map(y->imag(y), root_locus)

        for i = 1:length(root_locus[1])
            rl_real_i = map(x->x[i], rl_real)
            rl_imag_i = map(x->x[i], rl_imag)

            # Get the next color.
            color = pybuiltin("next")(ax[:_get_lines][:prop_cycler])["color"]

            ax[:plot](rl_real_i, rl_imag_i, linewidth=plw, color=color)
        end

        # Mark poles and zeros.
        ( length(p) > 0 ) &&
        ax[:scatter](map(x->real(x), p), map(y->imag(y), p),
                     s=mksz, marker="x", linewidths=plw,
                     edgecolor="red", color="red")

        ( length(z) > 0 ) &&
        ax[:scatter](map(x->real(x), z), map(y->imag(y), z),
                     s=mksz, marker="o", linewidths=plw,
                     edgecolor="blue", facecolor="none")

        # Asymptotes
        # ==========

        # Compute the number of asymptotes.
        num_asymptotes = length(mp) - length(mz)

        σ_a = (sum(mp) - sum(mz))/num_asymptotes
        θ_a = [ (2*m+1)*pi/num_asymptotes for m in 0:(num_asymptotes-1) ]

        # Get the X and Y axis limits.
        xl = ax[:get_xlim]()
        yl = ax[:get_ylim]()

        # Get the limits in polar coordinates.
        θl_1 = mod(atan2(yl[2], real(xl[2]-σ_a)), 2*pi)
        θl_2 = mod(atan2(yl[2], real(xl[1]-σ_a)), 2*pi)
        θl_3 = mod(atan2(yl[1], real(xl[1]-σ_a)), 2*pi)
        θl_4 = mod(atan2(yl[1], real(xl[2]-σ_a)), 2*pi)

        for i = 1:length(θ_a)
            # Start point.
            x_0 = real(σ_a)
            y_0 = imag(σ_a)

            if (θ_a[i] <= θl_1) || (θ_a[i] >= θl_4)
                x_f = xl[2]
                y_f = y_0 + (x_f - x_0)*tan(θ_a[i])
            elseif (θ_a[i] >= θl_2) && (θ_a[i] <= θl_3)
                x_f = xl[1]
                y_f = y_0 + (x_f - x_0)*tan(θ_a[i])
            elseif (θ_a[i] > θl_1) && (θ_a[i] < θl_2)
                y_f = yl[2]
                x_f = x_0 + (y_f - y_0)/tan(θ_a[i])
            else
                y_f = yl[1]
                x_f = x_0 + (y_f - y_0)/tan(θ_a[i])
            end

            ax[:plot]([x_0;x_f], [y_0;y_f], linestyle="--", color="#BBBBBB")
        end
    end

    nothing
end

"""
### function rlocus_on_mouse_press(event)

Function called if the used click on a root locus plot. This is only available
for PyPlot

##### Args

* event: Click event.

"""

function rlocus_on_mouse_press(event)
    rlocus_data = rlocus_dict[event[:canvas]]

    root_locus   = rlocus_data[1]
    root_locus_k = rlocus_data[2]

    # The click point inside the plot.
    x = event[:xdata]
    y = event[:ydata]

    # If `x`  or `y` is `nothing`, then the user clicked outside the plot area.
    ( (x == nothing) || (y == nothing) ) && return

    # If there is no point ploted in the root locus, then nothing should be
    # done.
    (isempty(root_locus))    && return

    # If the system does not have any closed-loop poles, then nothing should be
    # done.
    (isempty(root_locus[1])) && return

    # Get the closest point in the graphic to the click point.
    ind1  = indmin(map(i->minimum(abs.(i - Complex(x,y))), root_locus))
    ind2  = indmin(abs.(root_locus[ind1] - Complex(x,y)))
    poles = root_locus[ind1]
    gain  = root_locus_k[ind1]

    # Print the information about the poles.
    println("")
    print((@bold)*(@green))
    println("Closed-Loop Poles")
    print(@creset)
    println("=================")
    println("")
    print("Selected gain: ")
    print(@treset)
    println("$gain")
    print(@bold)
    println("Closed-loop poles (selected one is marked in yellow):")
    println(@treset)
    pole_info(poles; highlight_ind=ind2)

    nothing
end


