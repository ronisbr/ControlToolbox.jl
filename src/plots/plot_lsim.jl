"""
### function plot_lsim(y, t)

Plot the result computed by `lsim` function.

##### Args

* y: Output vector.
* t: Time vector.

##### Parameters

* title: Title of the figure.

"""

function plot_lsim(y, t; title="")
    # Get the number of outputs in y.
    num_outputs = size(y[1])[1]

    # Compute the grid of the subplots.
    nx = (num_outputs > 3) ? 3 : num_outputs
    ny = (num_outputs > 3) ? num_outputs%3 : 1

    if (plot_system == :pyplot)
        fig = plt.figure(title)

        ax = Array{Any}(num_outputs)

        for i = 1:num_outputs
            ax[i] = fig[:add_subplot](nx, ny, i)
        end

        for i = 1:num_outputs
            ax[i][:plot](t, map(x->x[i], y), linewidth=plw)
            ax[i][:set_title]("Y$i")
            ax[i][:set_xlabel]("Time")
            ax[i][:set_ylabel]("Output")

            # Change figure border width.
            ax[i][:spines]["top"][:set_linewidth](blw)
            ax[i][:spines]["bottom"][:set_linewidth](blw)
            ax[i][:spines]["left"][:set_linewidth](blw)
            ax[i][:spines]["right"][:set_linewidth](blw)
        end

        fig[:tight_layout]()
    end
end

"""
### function plot_lsim(y, t, u::Function)

Plot the result computed by `lsim` function together with the input function
`u`.

##### Args

* y: Output vector.
* t: Time vector.
* u: Input function.

##### Parameters

* title: Title of the figure.

##### Remarks

It is assumed that `u` function always returns a scalar.

"""

function plot_lsim(y, t, u::Function; title="")
    # Get the number of outputs in y.
    num_outputs = size(y[1])[1]

    # Compute the grid of the subplots.
    nx = (num_outputs > 3) ? 3 : num_outputs
    ny = (num_outputs > 3) ? num_outputs%3 : 1

    if (plot_system == :pyplot)
        fig = plt.figure(title)

        ax = Array{Any}(num_outputs)

        for i = 1:num_outputs
            ax[i] = fig[:add_subplot](nx, ny, i)
        end

        for i = 1:num_outputs
            # Plot the input function.
            ax[i][:plot](t, map(t->u(t), t),
                         color="#BBBBBB", linewidth=plw, linestyle="--")

            ax[i][:plot](t, map(x->x[i], y), linewidth=plw)
            ax[i][:set_title]("Y$i")
            ax[i][:set_xlabel]("Time")
            ax[i][:set_ylabel]("Output")

            # Change figure border width.
            ax[i][:spines]["top"][:set_linewidth](blw)
            ax[i][:spines]["bottom"][:set_linewidth](blw)
            ax[i][:spines]["left"][:set_linewidth](blw)
            ax[i][:spines]["right"][:set_linewidth](blw)
        end

        fig[:tight_layout]()
    end
end
