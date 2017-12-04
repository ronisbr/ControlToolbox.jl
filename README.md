# ControlToolbox.jl

This is a Control Toolbox for design and test control systems using Julia
language. The package development started in 2014 and was created to solve two
problems that the other packages I found on Internet could not solve:

1. The continuous systems (transfer functions and state-space) must be simulated
   using an ODE solver instead of the usual approach to convert them to a
   discrete system. Hence, we could select the solver algorithm and tolerance.
2. The root locus plot must be interactive so that the design can be simplified.

The first point was easily achieved by using the package **OrdinaryDiffEq.jl**.
The second was only recently possible due to the following bug:

https://github.com/JuliaPy/PyPlot.jl/issues/21

Hence, I decided to make this package public! Notice that it is still
experimental and many improvements must be done (help is welcome :) ).

**NOTE**: This package was not initally created prioritizing code performance.
Hence, some functions can still be optimized to decrease the computational
workload.

## Characteristics

The following major characteristics is already implemented:

* Creation of Space-States and Transfer Functions.
* System simulation with different kinds of inputs:
    * Step;
    * Ramp;
    * Impulse;
    * User-defined.
* Graphical analysis tools:
    * Bode;
    * Nyquist;
    * Root locus.

## Compatibility

The plotting system, including the interactive root locus plot, was tested under
the following systems:

|                | Linux (openSUSE Tumbleweed) | macOS 10.13 |
|----------------|:---------------------------:|:-----------:|
| Python 2 + Qt4 | X                           |             |
| Python 2 + Qt5 | X                           |             |
| Python 2 + Tk  | X                           |             |
| Python 3 + Qt4 | X                           |             |
| Python 3 + Qt5 | X                           | X           |
| Python 3 + Tk  | X                           |             |

However, the recommended setup is **Python 3 with Qt5**.

P.S.: If you test this package under other environment, please, let me know so
that I can update this table!

## Installation and Initialization

To install the package, execute the following command:

    Pkg.clone("https://github.com/ronisbr/ControlToolbox.jl")

To use the package, execute the following command:

    using ControlToolbox

If you want to use the graphical analysis tools, then you **must** install
**PyPlot** (https://github.com/JuliaPy/PyPlot.jl) and **matplotlib**
(https://matplotlib.org/users/installing.html) with a backend (the Qt5 backend
is highly recommended). After this, the graphical system can be initialized by":

    use_pyplot(:qt5)

Change `:qt5` to the selected backend. For more information, see `pygui`
function in **PyCall** package (https://github.com/JuliaPy/PyCall.jl). 

## Examples

### Creating a state-space

```julia
A = [0 1; -1 -1]
B = [0; 1]
C = [1 1]
D = []

sys = ss(A,B,C,D)
```

### Creating a transfer function

```julia

# Creating a transfer function with a numerator and denominator.
G = tf([1;2], [1;1;2;3;4])


                   1.0*s + 2.0
    -----------------------------------------
    1.0*s^4 + 1.0*s^3 + 2.0*s^2 + 3.0*s + 4.0

# Creating a transfer function with a gain, a set of zeros, and a set of poles.
G = zpk([-1],[-1;-2;-3;-4],2)

                     2.0*s + 2.0
    ---------------------------------------------
    1.0*s^4 + 10.0*s^3 + 35.0*s^2 + 50.0*s + 24.0

# Computing the minimal realization of a transfer function.
Gm = minreal(G)

                   2.0
    ---------------------------------
    1.0*s^3 + 9.0*s^2 + 26.0*s + 24.0
```

### Converting between state-space and transfer function

```julia
G = tf([1;1],[1;2;3;4;5])

                   1.0*s + 1.0
    -----------------------------------------
    1.0*s^4 + 2.0*s^3 + 3.0*s^2 + 4.0*s + 5.0

sys = tf2ss(G)

    A = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; -5.0 -4.0 -3.0 -2.0]

    B = [0.0, 0.0, 0.0, 1.0]

    C = [1.0 1.0 0.0 0.0]

    D = [0.0]

G2 = ss2tf(sys)

                   1.0*s + 1.0
    -----------------------------------------
    1.0*s^4 + 2.0*s^3 + 3.0*s^2 + 4.0*s + 5.0
```

### Analysis

```julia
# Root locus
G = tf([1;2;2],[1;9;33;51;26])

                1.0*s^2 + 2.0*s + 2.0
    --------------------------------------------
    1.0*s^4 + 9.0*s^3 + 33.0*s^2 + 51.0*s + 26.0

rlocus(G)
```
![Root Locus](https://github.com/ronisbr/ControlToolbox.jl/raw/master/figs/rlocus.gif "Root Locus")

```julia
# Bode diagram
bode(G)

# Nyquist diagram
nyquist(G)

# System simulation
step(G)
ramp(G)
impulse(G)
lsim(G, (t)->sin(t), [0;10], []; title="Simulation #1")

# Compute poles and invariant zeros.
pole(G)
    4-element Array{Complex{Float64},1}:
     -3.0+2.0im
     -3.0-2.0im
     -2.0+0.0im
     -1.0+0.0im

izero(G)
    2-element Array{Complex{Float64},1}:
     -1.0+1.0im
     -1.0-1.0im

# Pole information
pole_info(pole(G))
                      Pole |         Damping |   Freq. [rad/s] | Time Const. [s]
    -----------------------+-----------------+-----------------+----------------
              -3.0 + 2.0im |          3.6056 |         0.83205 |         0.33333
              -3.0 - 2.0im |          3.6056 |         0.83205 |         0.33333
              -2.0 + 0.0im |               2 |               1 |             0.5
              -1.0 + 0.0im |               1 |               1 |               1

```

## TODO

Currently, the TODO list is huge! :) The idea is to create a toolbox that can be
used by the students and researchers at my Institution to design and simulate
Control Systems.
