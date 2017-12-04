export  change_lsim_tolerances, lsim

################################################################################
#                            Global Configurations
################################################################################

# Tolerances for the ODE solver.
global abstol = 1e-8
global reltol = 1e-8

"""
### function change_lsim_tolerances(new_abstol, new_reltol)

Change the tolerances for the ODE solver used in `lsim()`. A lower tolerance
yields a more accurate simulation, whereas a higher tolerance yields a faster
simulation. For more information, see OrdinaryDiffEq.

##### Args

* new_abstol: New absolute tolerance.
* new_reltol: New relative tolerance.

"""

function change_lsim_tolerances(new_abstol, new_reltol)
    global abstol
    global reltol

    abstol = new_abstol
    reltol = new_reltol

    nothing
end

################################################################################
#                                     lsim
################################################################################

"""
### function lsim(sys::StateSpace, u::Function, t::Vector, x0::Vector=[], solver=Tsit5(), plot=true)

Simulate a dynamic system using an user-defined function.

##### Args

* sys: Dynamic system in state space form.
* u: Input function.
* t: Time vector. The solver will add the values in this vector to those in
  which the solution will be evaluated.
* x0: Initial state. If it is not present, the 0 will be used.

##### Parameters

* solver: Algorithm to solve the differential equation (See OrdinaryDiffEq).
  Default is `Tsit5()`.
* plot: Plot the result (default = `true`).
* title: Title of the figure.

##### Returns

* y: The output of the system.
* t: Array with the time instants in which the solver evaluated the system.
* x: State vector at each time instant.

"""

function lsim(sys::StateSpace,
              u::Function,
              t::Vector,
              x0::Vector=[];
              solver=Tsit5(),
              plot=true,
              title="")
    # Create the dynamic function.
    f(t,x) = sys.A*x + sys.B*u(t)

    # Check if the initial value is available.
    if (isempty(x0))
        x0 = zeros(size(sys.A,1))
    else
        ( length(x0) != size(sys.A,1) ) &&
        error("The dimension of x0 is not correct.")

        x0 = convert(Array{Float64,1},x0)
    end

    # Simulate the system.
    prob = ODEProblem(f, x0, ( Float64(t[1]), Float64(t[end]) ) )
    sol  = solve(prob,
                 solver,
                 tstops = convert(Vector{Float64}, t),
                 abstol=abstol,
                 reltol=reltol)

    # Compute the measurement.
    y = map(x->sys.C*x, sol.u) +
        map(u->sys.D*u, [u(k) for k in sol.t])

    if plot
        # If the size of input function `u` is 1, then plot it.
        ( length(u(0)) == 1 ) ? plot_lsim(y, sol.t, u; title=title) :
                                plot_lsim(y, sol.t; title=title)
    end

    # Return.
    y, sol.t, sol.u
end

# Initial state `x0` is a number.
function lsim(sys::StateSpace,
              u::Function,
              t::Vector,
              x0::Number;
              solver=Tsit5(),
              plot=true,
              title="")
    lsim(sys, u, t, [x0]; solver=solver, plot=plot)
end

"""
### function lsim(G::TransferFunction, u::Function, t::Vector, x0::Vector=[], solver=Tsit5(), plot=true)

Simulate a dynamic system using an user-defined function.

##### Args

* G: Transfer function of the dynamic system.
* u: Input function.
* t: Time vector. The solver will add the values in this vector to those in
  which the solution will be evaluated.
* x0: Initial state. If it is not present, the 0 will be used.

##### Parameters

* solver: Algorithm to solve the differential equation (See OrdinaryDiffEq).
  Default is `Tsit5()`.
* plot: Plot the result (default = `true`).
* title: Title of the figure.

##### Returns

* y: The output of the system.
* t: Array with the time instants in which the solver evaluated the system.
* x: State vector at each time instant.

"""

function lsim(G::TransferFunction,
              u::Function,
              t::Vector,
              x0::Union{Number,Vector}=[];
              solver=Tsit5(),
              plot=true,
              title="")
    # Convert the transfer function to state space and call lsim.
    lsim(tf2ss(G), u, t, x0; solver=solver, plot=plot)
end
