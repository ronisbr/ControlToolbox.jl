export use_pyplot

global plot_system = :none

################################################################################
#                         Selection of plotting system
################################################################################

"""
### function use_pyplot(pyplot_gui)

Select the PyPlot as the plotting system.

##### Args

* pyplot_gui: PyPlot GUI (see PyPlot).

"""

function use_pyplot(pyplot_gui)
    # List obtained from PyPlot.jl package.
    gui2matplotlib = Dict(:wx=>"WXAgg",:gtk=>"GTKAgg",:gtk3=>"GTK3Agg",
                          :qt_pyqt4=>"Qt4Agg",
                          :qt_pyqt5=>"Qt5Agg",
                          :qt_pyside=>"Qt4Agg",
                          :qt4=>"Qt4Agg",
                          :qt5=>"Qt5Agg",
                          :qt=>"Qt4Agg",:tk=>"TkAgg")

    # Import matplotlib.
    const matplotlib = pyimport("matplotlib")

    # Select the GUI and set to interactive mode.
    pygui_start(pyplot_gui)
    matplotlib[:use](gui2matplotlib[pyplot_gui])
    matplotlib[:interactive](true)

    # Import pyplot.
    pltm = pyimport("matplotlib.pyplot") # raw Python module
    global plt = pywrap(pltm)

    # Set the plotting system as PyPlot.
    global plot_system = :pyplot

    nothing
end

################################################################################
#                                    Plots
################################################################################

# Definitions

const blw   = 2.0   # Border line width.
const plw   = 2.5   # Plot line width.
const glw   = 1.3   # Grid line width.
const sglw  = 1.0   # Subgrid line width.
const tifsz = 14    # Title font size.
const mksz  = 80    # Marker size.
const arrsz = 20    # Arrow height.
const alw   = 1.5   # Axes line width.

################################################################################
#                                   Includes
################################################################################

include("plot_bode.jl")
include("plot_lsim.jl")
include("plot_nyquist.jl")
include("plot_rlocus.jl")
