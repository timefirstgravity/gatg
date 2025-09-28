# validate_lorentz_from_phi.sage
# SageMath script to verify the Lorentz-factor derivation from the time-first lapse Φ
#
# Metric (zero shift, spherical, radial motion):
#   ds^2 = -A dt^2 + A^{-1} dr^2  with  A = e^{2Φ(x,t)}
#
# Static clock:  dτ_stat = sqrt(A) dt
# Local ruler:   dl      = dr / sqrt(A)
# Local speed:   v = (dl/dτ_stat) = (dr/dt) / A
# Worldline:     dτ^2 = A dt^2 - A^{-1} dr^2 = A dt^2 (1 - v^2)
# Hence:         dτ = dτ_stat * sqrt(1 - v^2),  γ = 1/sqrt(1 - v^2)

from sage.all import var, function, exp, sqrt, simplify, SR

# Coordinates and differentials (symbolic)
t, r = var('t r', domain='real')
dt, dr = var('dt dr', domain='real')

# Lapse exponent Φ(t,r) and A = e^{2Φ}
Phi = function('Phi')(t, r)
A   = exp(2*Phi)

# 1) Static observer proper time and local radial ruler
dtaus_stat = sqrt(A) * dt                # dτ_stat
dl         = dr / sqrt(A)                # local proper radial length element

# 2) Define coordinate speed u := dr/dt (symbolic placeholder)
u = var('u', domain='real', latex_name='u')

# Local speed measured by static observers: v = (dl/dτ_stat)
# Carefully perform the ratio so that dt, dr remain symbolic, then substitute dr = u*dt
v_expr = (dl / dtaus_stat).subs(dr == u*dt).simplify_full()
# Expected: v = u / A
v_expected = (u / A).simplify_full()

# 3) Proper time along moving worldline
# For pure radial motion (dΩ = 0):  dτ^2 = A dt^2 - A^{-1} dr^2
dtausq = (A * dt**2 - (1/A) * dr**2).simplify_full()

# Substitute dr = u dt to express dτ^2 in terms of u, then also in terms of v
dtausq_u = dtausq.subs(dr == u*dt).simplify_full()           # A dt^2 - (1/A) u^2 dt^2
dtausq_v = (A * dt**2 * (1 - v_expr**2)).simplify_full()     # A dt^2 (1 - v^2)

# 4) Check key identities

def check_zero(expr, label):
    expr_simplified = expr.simplify_full()
    ok = (expr_simplified == 0)
    print(f"[{label}] ->", "OK" if ok else f"FAILED (residual: {expr_simplified})")
    return ok

all_ok = True

# v = u / A
all_ok = all_ok and check_zero((v_expr - v_expected), "v == u/A (local speed match)")

# dτ^2 with u vs with v: A dt^2 - A^{-1} u^2 dt^2  ==  A dt^2 (1 - v^2)
all_ok = all_ok and check_zero((dtausq_u - dtausq_v), "dτ^2(u) == dτ^2(v)")

# Factorization: dτ^2 = (dτ_stat)^2 * (1 - v^2)
lhs = dtausq.subs(dr == u*dt)  # Make sure we substitute consistently
rhs = (dtaus_stat**2 * (1 - v_expr**2)).simplify_full()
all_ok = all_ok and check_zero((lhs - rhs), "dτ^2 == (dτ_stat)^2 * (1 - v^2)")

# 5) Report γ and the final relation dτ = dτ_stat * sqrt(1 - v^2)
# (Symbolic statement; the checks above already validated the squared relation.)
gamma = (1 / sqrt(1 - v_expr**2)).simplify_full()

print("\nSummary:")
print("  A(t,r) = e^{2Φ(t,r)}")
print(f"  v = {v_expr} = u/A  (local speed)")
print("  γ = 1/sqrt(1 - v^2)  (Lorentz factor)")
print("  dτ = dτ_stat * sqrt(1 - v^2)  (proper time relation)")
print("\nPhysical interpretation:")
print("  - Static observers measure proper time dτ_stat = √A dt")
print("  - Moving objects with coordinate speed u = dr/dt have local speed v = u/A")
print("  - The usual special relativistic time dilation applies with this local speed")

print("\nOverall:", "ALL CHECKS PASSED ✅" if all_ok else "Some checks FAILED ❌")
