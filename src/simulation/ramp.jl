export ramp

################################################################################
#                                     ramp
################################################################################

"""
### function ramp(sys::StateSpace, t::Vector = []; solver=Tsit5(), plot=true)

Compute the system response to the unitary ramp input.

##### Args

* sys: Dynamic system in state space form.
* t: Time vector. The solver will add the values in this vector to those in
  which the solution will be evaluated. If it is empty, then the simulation will
  occur up to 10x the slowest pole.

##### Returns

* y: The output of the system.
* t: Array with the time instants in which the solver evaluated the system.
* x: State vector at each time instant.

##### Remarks

If the system has more than one input, then the step simulation will be applied
to each input. In this case, the output vectors `y`, `t`, and `x` will be arrays
of the solution data. Hence, in this case, the data from the first input can be
accessed by mean of:

    y[1]
    t[1]
    x[1]

Moreover, a new plot window will be created for each simulated input.


"""

function ramp(sys::StateSpace, t::Vector = []; solver=Tsit5(), plot=true)
    # Check if the time vector was specified.
    if (isempty(t))
        # If not, simulate the system up to 10x the slowest pole.
        poles = pole(sys)

        # Get the slowest frequency of the system.
        slowest_freq = minimum(abs.(real(poles)))

        t_final = (slowest_freq != 0) ? 10/slowest_freq : 40

        # Do not simulate the entire interval if t_final is too large.
        (t_final > 1000) && (t_final = 1000)

        t = [0,t_final]
    end

    # Simulate the system considering each input.
    num_inputs = size(sys.B,2)

    if num_inputs == 1
        lsim(sys, (t)->(t > 0) ? t : 0, t, [];
             solver=solver, plot=plot, title="Ramp Response")
    else
        # Outputs.
        yo = []
        to = []
        uo = []

        for i = 1:num_inputs
            yi, ti, ui = lsim(sys, (t)->begin
                                  u = zeros(num_inputs)
                                  (t > 0) && (u[i] = t)
                                  u
                              end, t, [];
                              solver=solver, plot=plot, title="Ramp Response - From Input #$i")
            push!(yo,yi)
            push!(to,ti)
            push!(uo,ui)
        end
        yo, to, uo
    end
end

"""
### function ramp(G::TransferFunction, t::Vector = []; solver=Tsit5(), plot=true)

Compute the system response to the unitary ramp input.

##### Args

* G: Transfer function of the dynamic system.
* t: Time vector. The solver will add the values in this vector to those in
  which the solution will be evaluated. If it is empty, then the simulation will
  occur up to 10x the slowest pole.

##### Returns

* y: The output of the system.
* t: Array with the time instants in which the solver evaluated the system.
* x: State vector at each time instant.

"""

function ramp(G::TransferFunction, t::Array = []; solver=Tsit5(), plot=true)
    # Convert the transfer function to state space and call step.
    ramp(tf2ss(G), t; plot=plot, solver=solver)
end
