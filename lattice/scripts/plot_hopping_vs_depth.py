"""
Plot hopping versus depth of lattice using a variety of calculation techniques
"""
import numpy as np
from lattice_single_particle import *
import datetime
import matplotlib.pyplot as plt

# #################################
# general parameters for the calculation
# #################################
# depths = np.array([3, 7, 10, 20, 100])
depths = np.logspace(0, 2, 10)

theta = 91.6267
r = 0.4

hoppings_bs, hoppings_bw, hoppings_wannier, hoppings_wannier_dens_dep, wannier_integral, energies, wannier_fn, potential, xpts, ypts = \
        solve_d4_lattice(depths, theta, 1, 0, dx=0.1, n1_sites=12, n2_sites=12, n1_recp_max=10, n2_recp_max=10, normalize_er=True)

# plot hopping vs. depth
figh_log = plt.figure()
nrows = 2
ncols = 2

# potential
plt.subplot(nrows, ncols, 1)
plt.imshow(potential)
plt.title('Potential')

# wannier function
plt.subplot(nrows, ncols, 3)
plt.imshow(np.abs(wannier_fn[:, :, 0])**2)
plt.title('Wannier function |W(x,y)|^2')

# nearest-neighbor hopping
plt.subplot(nrows, ncols, 2)
# x
ph_x = plt.semilogy(depths, hoppings_bs[:, 0], '.')

h2 = plt.semilogy(depths, hoppings_bw[:, 0], 'x')
h2[0].set_color(ph_x[0].get_color())

ph_xw = plt.semilogy(depths, hoppings_wannier[:, 1], '+')
ph_xw[0].set_color(ph_x[0].get_color())

# y
ph_y = plt.semilogy(depths, hoppings_bs[:, 1], '.')

h2 = plt.semilogy(depths, hoppings_bw[:, 2], 'x')
h2[0].set_color(ph_y[0].get_color())

h3 = plt.semilogy(depths, np.abs(hoppings_wannier[:, 1]), '+')
h3[0].set_color(ph_y[0].get_color())

plt.title('NN hopping')
plt.xlabel('Depth (Er)')
plt.ylabel('Hopping (Er)')
leg = plt.legend([ph_x[0], ph_xw[0], ph_y[0], h3[0]], ['x bs', 'x wann', 'y bs', 'y wann'])
leg.draggable()
plt.grid()

# diagonal and second neighbor hopping
plt.subplot(nrows, ncols, 4)

# x
ph_xnnn = plt.semilogy(depths, np.abs(hoppings_bs[:, 3]), '.')

h3 = plt.semilogy(depths, np.abs(hoppings_wannier[:, 3]), '+')
h3[0].set_color(ph_xnnn[0].get_color())

# y
ph_ynnn = plt.semilogy(depths, np.abs(hoppings_bs[:, 4]), '.')

h3 = plt.semilogy(depths, np.abs(hoppings_wannier[:, 4]), '+')
h3[0].set_color(ph_ynnn[0].get_color())

# diag
ph_d = plt.semilogy(depths, hoppings_bs[:, 2], '.')

h2 = plt.semilogy(depths, np.abs(hoppings_bw[:, 3]), 'x')
h2[0].set_color(ph_d[0].get_color())

h3 = plt.semilogy(depths, np.abs(hoppings_wannier[:, 2]), '+')
h3[0].set_color(ph_d[0].get_color())

plt.grid()
plt.xlabel('Lattice Depth (Er)')
plt.ylabel('Hopping (Er)')
leg = plt.legend([ph_xnnn[0], ph_ynnn[0], ph_d[0], h3[0]],
                 ['x2n', 'y2n', 'diag bs', 'diag wannier'])
leg.draggable()
plt.show()

# save results
if 0:
    now = datetime.datetime.now()
    date_str = "%04d;%02d;%02d_%02dh_%02dm" % (now.year, now.month, now.day, now.hour, now.minute)
    fname = "%s_log_hopping_v_latt_depths.png" % date_str
    figh_log.savefig(fname)
