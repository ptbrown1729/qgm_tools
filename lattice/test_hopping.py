from lattice_single_particle import *

theta = 91.6267
r = 0.6
# theta = 90
# r = 1
alpha = 0
depth = 4

hoppings_bs, hoppings_bw, hoppings_wannier, hoppings_wannier_dens_dep, wannier_integral, energies, wannier_fn, potential, xpts, ypts = \
        solve_d4_lattice(depth, theta, r, alpha, dx=0.05, n1_sites=13, n2_sites=13, n1_recp_max=10, n2_recp_max=10, normalize_er=True)

print("bandw tx = %.5e Er" % hoppings_bw[1])
print("bandstru tx = %.5e Er" % hoppings_bs[0])
print("wannier tx = %0.5e Er" % hoppings_wannier[0])
print("")
print("bandw ty = %.5e Er" % hoppings_bw[2])
print("bandstru ty = %.5e Er" % hoppings_bs[1])
print("wannier ty = %0.5e Er" % hoppings_wannier[1])
print("")
print("bandw td = %.5e Er" % hoppings_bw[3])
print("bandstru td = %.5e Er" % hoppings_bs[2])
print("wannier td = %0.5e Er" % hoppings_wannier[2])

