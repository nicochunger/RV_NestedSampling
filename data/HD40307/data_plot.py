import pandas as pd
import matplotlib.pyplot as plt


# ----------- Data from 2003-2015 -------------------
# filepath = "HD40307_vizier.txt"

# data = pd.read_csv(filepath, sep='\t', comment='#', skiprows=[1, ])

# time = data["rjd"]
# rv = data["vrad"]
# err = data["svrad"]
# rhk = data["rhk"]

# plt.figure(0)
# plt.errorbar(time, rv, yerr=err, fmt='r.', ecolor='k')
# plt.xlabel("Tiempo (BJD-2,400,000)")
# plt.ylabel("Velocidad radial (m/s)")


# plt.figure(1)
# plt.plot(time, rhk, 'k.')
# plt.xlabel("Tiempo (BJD-2,400,000)")
# plt.ylabel(r"$\log(R'_{HK})$")

# plt.show()
# ---------------------------------------------------


# ---------- Both Datasets (2003-2018) ---------------
H03_filepath = "HD40307_HARPS03(DRS-3-5).rdb"
H15_filepath = "HD40307_HARPS15(DRS-3-5).rdb"

H03 = pd.read_csv(H03_filepath, sep='\t', comment='#', skiprows=[1, ])
H15 = pd.read_csv(H15_filepath, sep='\t', comment='#', skiprows=[1, ])

time_H03 = H03["rjd"][6:]
rv_H03 = H03["vrad"][6:] - 31334.664
err_H03 = H03["svrad"][6:]
rhk_H03 = H03["rhk"][6:]

time_H15 = H15["rjd"]
rv_H15 = H15["vrad"] - 31349.31
err_H15 = H15["svrad"]
rhk_H15 = H15["rhk"]

plt.figure(0)
plt.errorbar(time_H03, rv_H03, yerr=err_H03,
             fmt='r.', ecolor='k', label="HARPS03")
# plt.errorbar(time_H15, rv_H15, yerr=err_H15,
#              fmt='y.', ecolor='k', label="HARPS15")
plt.xlabel("Tiempo (BJD-2,400,000)", fontsize=12)
plt.ylabel(r"$\Delta$ Velocidad radial (m/s)", fontsize=12)
# plt.legend(loc=2)


plt.figure(1)
plt.plot(time_H03, rhk_H03, 'k.', label="HARPS03")
# plt.plot(time_H15, rhk_H15, 'y.', label="HARPS15")
plt.xlabel("Tiempo (BJD-2,400,000)", fontsize=12)
plt.ylabel(r"$\log(R'_{HK})$", fontsize=12)
# plt.legend(loc=2)

plt.show()
