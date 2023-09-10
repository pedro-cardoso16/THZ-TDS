import thztds.spectrum as spec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
d = 600
n1 = 2 + 4j
n2 = 1
w = 1

material_ex = spec.Medium("example", d)
p_ex = spec.PropagationTerms(n1, n2, 1, d)


n1_set = np.mgrid[0:5:0.1,0:9:0.1]
n1_set = n1_set[0] + 1j*n1_set[1]

Tc_p = spec.PropagationTerms(n1_set, n2, 1, d)
Tc, Tc_arg = Tc_p.total_propagation_term(), Tc_p.argument()


Tm = p_ex.total_propagation_term()
Tm_arg = p_ex.argument()


y = spec.error(Tc, Tc_arg, Tm, Tm_arg)
x1,x2 = np.real(n1_set), np.imag(n1_set)

fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
surf = ax.plot_surface(x1,x2,y,cmap=cm.coolwarm)
ax.set_xlabel("n")
ax.set_ylabel("k")
ax.set_zlabel(r"$\delta$")
plt.show()