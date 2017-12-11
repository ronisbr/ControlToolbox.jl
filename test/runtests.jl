using ControlToolbox
using Polynomials
using Base.Test

################################################################################
#                              Transfer Function
################################################################################

# Basic Operations
# ================

G = tf([1],[1;1])
H = tf([1],[1;2])

G1 = G-2G

@test G1.num.a == [-1]
@test G1.den.a == [1;1]

G1 = minreal(G + H)

@test G1.num.a == [3;2]
@test G1.den.a == [2;3;1]

G1 = G + H
G2 = parallel(G,H)

@test G1.num == G2.num
@test G1.den == G2.den

G1 = G*H
G2 = series(G,H)

@test G1.num == G2.num
@test G1.den == G2.den

G1 = minreal(G/(1+G*H))
G2 = feedback(G,H)

@test G1.num == G2.num
@test G1.den == G2.den

G = zpk([-1;-1],[+1;+1],4)

@test G.num == 4*poly([-1;-1])
@test G.den == 1*poly([+1,+1])
