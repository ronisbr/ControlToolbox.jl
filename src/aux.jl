export  issiso

################################################################################
#                             Auxiliary Functions
################################################################################

"""
### function issiso(sys::StateSpace)

Check if a system is SISO (single input, single output).

##### Args

* sys: System that will be checked in state space representation.

##### Returns

`true` if the system is SISO, `false` otherwise.

"""

function issiso(sys::StateSpace)
    num_inputs  = size(sys.B,2)
    num_outputs = size(sys.C,1)

    ( ( num_inputs == 1 ) && ( num_outputs == 1 ) )
end

################################################################################
#                          Text Formatting Functions
################################################################################

"""
### macro bold()

Bold text.

"""

macro bold()
    "\x1b[1m"
end

"""
### macro boldn()

Bold text with default foreground color.

"""

macro boldn()
    "\x1b[1m\x1b[39m"
end

"""
### macro normal()

Normal text (removes bold).

"""

macro normal()
   "\x1b[22m"
end

"""
### macro green()

Green text.

"""

macro green()
    "\x1b[32m"
end

"""
### macro yellow()

Yellow text.

"""

macro yellow()
    "\x1b[33m"
end

"""
### macro creset()

Reset text color, i.e. change color to the default foreground color.

"""

macro creset()
    "\x1b[39m"
end

"""
### macro treset()

Reset text color and text formatting.

"""

macro treset()
    "\x1b[0m"
end

