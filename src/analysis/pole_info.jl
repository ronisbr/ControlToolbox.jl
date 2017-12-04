export pole_info

"""
### function pole_info_header()

Print the header of the table that shows the information about the poles:

                      Pole |         Damping |   Freq. [rad/s] | Time Const. [s]
    -----------------------+-----------------+-----------------+----------------

"""

function pole_info_header()
    println(((@bold)) *
            lpad("Pole",            22, " ") * " |" *
            lpad("Damping",         16, " ") * " |" *
            lpad("Freq. [rad/s]",   16, " ") * " |" *
            lpad("Time Const. [s]", 16, " ") *
            ((@normal)))

    println(((@bold)) * lpad("", 22, "-") * "-+" *
            lpad("", 16, "-")           * "-+" *
            lpad("", 16, "-")           * "-+" *
            lpad("", 16, "-")           *
            ((@normal)))
end

"""
### function pole_info(pole::Complex128; header = true, highlight = false)

Print the pole information:

* Pole: The location of the pole in the complex plane.
* Damping: The damping of the pole.
* Frequency: The natural frequency of the pole [rad/s].
* Time constant: The time constant of the pole [s].

##### Args

* pole: The pole.

##### Parameters

* header: If `true`, then the header will be printed (see `pole_info_header()`).
* highlight: If `true`, then the text will be printed in yellow.

"""

function pole_info(pole::Complex128; header = true, highlight = false)
    # Compute the damping, frequency, and time constant of the pole.
    ω_n = norm(pole)
    ζ = -real(pole)/ω_n
    t = 1/(ω_n*ζ)

    # Convert the numbers into text.
    pole_str = lpad("$(round(pole,5))",      22, " ")
    ω_n_str  = lpad((@sprintf "%16.5g" ω_n), 16, " ")
    ζ_str    = lpad((@sprintf "%16.5g" ζ),   16, " ")
    t_str    = lpad((@sprintf "%16.5g" t),   16, " ")

    (header) && pole_info_header()

    println(( highlight ? (@yellow) : "") * pole_str * ((@boldn)) * " |" * ((@normal)) *
            ( highlight ? (@yellow) : "") * ω_n_str  * ((@boldn)) * " |" * ((@normal)) *
            ( highlight ? (@yellow) : "") * ζ_str    * ((@boldn)) * " |" * ((@normal)) *
            ( highlight ? (@yellow) : "") * t_str    * ((@treset)))

    nothing
end

"""
### function pole_info(poles::Array{Complex128,1}; highlight_ind = 0)

Print the information of an array of poles.

##### Args

* poles: The array of poles.

##### Parameters

* highlight_ind: The index of a pole to be highlighted.

"""

function pole_info(poles::Array{Complex128,1}; highlight_ind = 0)
    pole_info_header()

    for i = 1:length(poles)
        if (i == highlight_ind)
            pole_info(poles[i]; header = false, highlight = true)
        else
            pole_info(poles[i]; header = false)
        end
    end

    nothing
end
