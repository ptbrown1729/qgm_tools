# run clusters and store results in MySql database
import os
import glob
import re
import numpy as np
import mysql.connector
import convert_qmc_file

# dqmc data path. Paths with in and out dirs.
dqmc_paths = ["Z://QuestQMC//dqmc_data//2019_05_14",
              "Z://QuestQMC//dqmc_data//2019_05_13",
              "Z://QuestQMC//dqmc_data//2019_05_12",
              "Z://QuestQMC//dqmc_data//repulsive_hubbard_dqmc_fitting_results//2018_08_21",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=0,1,9to11,T=0.15to0.35,mu=-5to-1,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=16,T=0.3to10,mu=-3to-1.25,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=2to8,T=0.15to0.35,mu=-3to1,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=2to8,T=0.15to0.35,mu=-5to3.25,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=0.3to10,mu=-3to-1.25,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=0.5to15,mu=-120to0,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=1to15,mu=-5to0,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=9to15,mu=-18.5to0,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=9to15,mu=-30to-19,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=9to11,T=0.2to0.35,mu=-1to1,npass=50000",
              "Z://QuestQMC//dqmc_data//VaryFilling//sign_test",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=12_test",
              "Z://QuestQMC//dqmc_data//VaryFilling//n=16_test",]

# # paths with more complicated issues...
# dqmc_paths = ["Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=0.25to0.45,mu=-4.5to0,npass=50000",
#               "Z://QuestQMC//dqmc_data//VaryFilling//n=8,U=8,T=0.3,mu=-3to1.5,npass=10x20000",]

# dqmc_paths = ["Z://QuestQMC//dqmc_data//FFLO//n=50x1,U=-1to-5,T=0.05to0.4,ty=0,dmu=0to8,mu=0"]

# database info
db_name = "dqmc_data"
table_name = "data"

# connect to mysql server
db = mysql.connector.connect(
    host = "phy-wbakr4.princeton.edu",
    user = "ptbrown",
    passwd = "b4krl4b",
    database = db_name
)

cursor = db.cursor()

# SHOW DATABASES
# USE dqmc_data

# drop table, for the moment as I am testing this repeatedly...
# query = "DROP TABLE %s" % table_name
# cursor.execute(query)


# create table if it doesn't exist
query = 'SHOW TABLES LIKE "%s"' % table_name
cursor.execute(query)
result = cursor.fetchall()
if result == []:

    query = 'CREATE TABLE data (' \
            'id INT(11) NOT NULL AUTO_INCREMENT PRIMARY KEY,' \
            'fname_out TEXT, dir_out TEXT, fname_in TEXT, dir_in TEXT, fname_geom TEXT, dir_geom TEXT,' \
            'year INT(4), month INT(2), day INT(2),'\
            'nsites INT(4), nx INT(4), ny INT(4), periodic_bc_x BOOLEAN, periodic_bc_y BOOLEAN,' \
            'u DOUBLE, T DOUBLE, tups DOUBLE, tdns DOUBLE, mu_up DOUBLE, mu_dn DOUBLE,' \
            'dtau DOUBLE, time_slices INT(4),' \
            'nwarmup INT(10), nmeasure INT(10),  nbin INT(10),' \
            'seed INT(10), frq_measurement INT(10), frq_recomputing_g INT(10),' \
            'naccept INT(10), nreject INT(10), gamma DOUBLE,' \
            'sign_avg DOUBLE, sign_avg_unc DOUBLE,' \
            'sign_up_avg DOUBLE, sign_up_avg_unc DOUBLE,' \
            'sign_dn_avg DOUBLE, sign_dn_avg_unc DOUBLE,' \
            'nup DOUBLE, nup_unc DOUBLE, ndn DOUBLE, ndn_unc DOUBLE' \
            ')'
    cursor.execute(query)


# print tables
cursor.execute('SHOW TABLES')
tables = cursor.fetchall()
print(tables)

# print description
cursor.execute('DESC %s' % table_name)
print(cursor.fetchall())

for jj, dqmc_path in enumerate(dqmc_paths):
    dqmc_out_path = os.path.join(dqmc_path, "out")
    dqmc_in_path = os.path.join(dqmc_path, "in")

    # load dqmc data
    # TODO: modify so can handle output files w/o input files

    dqmc_input_files = glob.glob(os.path.join(dqmc_in_path, "*.in"))
    for ii, fpath_in in enumerate(dqmc_input_files):
        print("%d/%d for path %d/%d" % (ii + 1, len(dqmc_input_files), jj + 1, len(dqmc_paths)) )

        fname_in = os.path.basename(fpath_in)

        # load input file data
        input_dict = convert_qmc_file.parse_input_file(fpath_in)

        # find output file
        fname_out = input_dict["ofile"]
        # in case has some directory also, strip it off
        fname_out = os.path.basename(input_dict["ofile"])
        fpath_out = os.path.join(dqmc_out_path, fname_out + ".out")

        # find geometry file
        if "gfile" in input_dict.keys():
            fname_geom = os.path.basename(input_dict["gfile"])
            fpath_geom = os.path.join(dqmc_in_path, fname_geom)
        else:
            fpath_geom = ''
            fname_geom = ''

        # load file data
        # TODO: loading input file twice right now...
        dqmc_dat = convert_qmc_file.qmc_results.from_dqmc_text_file(fpath_out, fnames_geom=fpath_geom, fnames_input=fpath_in)

        if not np.any(np.isnan(dqmc_dat.settings.super_cell_vectors)):
            nx = int(dqmc_dat.settings.super_cell_vectors[0, 0])
            ny = int(dqmc_dat.settings.super_cell_vectors[1, 1])
            if dqmc_dat.settings.super_cell_vectors[1, 0] != 0 or dqmc_dat.settings.super_cell_vectors[0, 1] != 0:
                raise Exception("DQMC geometry not square.")
        else:
            nx = dqmc_dat.settings.nx
            ny = dqmc_dat.settings.ny

        # store in SQL database
        query = "INSERT INTO" + " " + table_name + " " + \
                "(fname_out, dir_out, fname_in, dir_in, fname_geom, dir_geom," \
                " nsites, nx, ny," \
                " u, T, tups, tdns," \
                " mu_up, mu_dn, dtau, time_slices," \
                " nwarmup, nmeasure, nbin," \
                " seed, frq_measurement, frq_recomputing_g," \
                "naccept, nreject, gamma," \
                "sign_avg, sign_avg_unc," \
                " sign_up_avg, sign_up_avg_unc," \
                "sign_dn_avg, sign_dn_avg_unc," \
                "nup, nup_unc," \
                "ndn, ndn_unc)" \
                " VALUES (\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\"," \
                " \"%s\", \"%s\", \"%s\"," \
                " %s, %s, %s, %s," \
                " %s, %s, %s, %s," \
                " %s, %s, %s," \
                " %s, %s, %s," \
                " %s, %s, %s," \
                " %s, %s," \
                " %s, %s," \
                " %s, %s," \
                " %s, %s," \
                " %s, %s" \
                ")"

        # TODO: protect against nans
        # terminaed @ 27/1680 in folder 5/17 because one vaule was nan
        values = (fname_out, dqmc_out_path, fname_in, dqmc_in_path, fname_geom, dqmc_in_path,
                  dqmc_dat.settings.nsites, nx, ny,
                  dqmc_dat.settings.u, 1./dqmc_dat.settings.beta, dqmc_dat.settings.t_up, dqmc_dat.settings.t_dn,
                  dqmc_dat.settings.mu_up, dqmc_dat.settings.mu_dn, dqmc_dat.settings.dtau, dqmc_dat.settings.n_time_slices,
                  dqmc_dat.settings.nwarmup, dqmc_dat.settings.nmeasure, dqmc_dat.settings.nbin,
                  dqmc_dat.settings.seed, dqmc_dat.settings.frq_measurement, dqmc_dat.settings.frq_recomputing_g,
                  dqmc_dat.settings.naccept, dqmc_dat.settings.nreject, dqmc_dat.settings.gamma,
                  dqmc_dat.static_measurements.avg_sign.data, dqmc_dat.static_measurements.avg_sign.unc,
                  dqmc_dat.static_measurements.avg_up_sign.data, dqmc_dat.static_measurements.avg_up_sign.unc,
                  dqmc_dat.static_measurements.avg_dn_sign.data, dqmc_dat.static_measurements.avg_dn_sign.unc,
                  dqmc_dat.static_measurements.n_up.data, dqmc_dat.static_measurements.n_up.unc,
                  dqmc_dat.static_measurements.n_dn.data, dqmc_dat.static_measurements.n_dn.unc,
                  )

        cursor.execute(query, values)
        db.commit()

    # look at output
    # query = "SELECT fname_out, u, T, mu_up, mu_dn, sign_avg, sign_avg_unc FROM data WHERE u=-3 ORDER BY T,  mu_up, mu_dn, ASC"

