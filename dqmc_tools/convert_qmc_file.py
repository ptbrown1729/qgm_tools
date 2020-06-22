from __future__ import print_function
import numpy as np
import os
import glob
import re
import fnmatch
import copy
import warnings

# functions to parse different QUEST input/output files
def parse_dqmc_output_file(fname):
    """
    Extract data from DQMC output file (static, not time dependent), and return the results as a set of dictionaries
    which can be used to instantiate a qmc_results class object. Using that object, it is easy to average, etc,
    these results.

    Note that different QUEST-QMC fortran programs produce slightly different output files_out. Therefore, this function
    handles several different cases. It may be necessary to add more cases.

    :param fname:
    :return:
    settings: dictionary of floating point numbers representing parameters of the qmc_avg simulation, such as chemical
    potential, temperature, etc.
    meas: dictionary of measurement class instances or equal_time_fn instances, representing qmc_avg measurements.
     Measurement instances are simple quantities, such as density, kinetic energy, etc. equal_time_fn instances
     represent measurements that have either spatial or k-space resolution, such as SzSz correlations, green's
     functions, etc.
    k_indices: NumPy array of integer indices specifying the unique k-vectors in the system. For example, an 8x8 system
     has only 15 unique k-vectors
    kvects_reps: a list of NumPy arrays which give the representations of each unique k-vector in the Brillouin zone.
     This accounts for symmetry, so for a square system (0, pi/2) and (pi/2, 0) (e.g.) are equivalent

    """
    settings = {}
    meas = {}
    kvects_reps = []
    k_indices = []

    # regular expressions to match various lines in the DQMC output text files_out
    # add a new expression here (or extend an old one) to handle the new case.

    # matches a decimal number with or without +/-, witht or without decimal point.
    number_exp = "([+-]?\d+[.]?\d*)"

    # matches a number in scientific notation
    sci_not_exp = number_exp + "E" + number_exp
    settings_exp = " *(.*) : +%s$" % number_exp
    # number with uncertainty
    meas_exp = " *(.*) : +%s *\+- *%s$" % (sci_not_exp, sci_not_exp)

    # different possible real-space expressions depending on what Quest-QMC
    # Fortran program we are using. Possibly can come up with a nice general expression that encompasses all
    # (?:...) is a non-capturing group
    # This should match files_out using the "ggeom" program
    gfun_real_space_exp = " *%s {1,2}\d+[*]{0,2}(?: \d+)? {2,3}%s {3}%s {3}%s *%s *\+- *%s *$" % (number_exp,
                                                                                number_exp, number_exp, number_exp,
                                                                                sci_not_exp, sci_not_exp)
    # e.g.  0 064   0.00   0.00   0.00     0.50000000E+00 +-  0.00000000E+00
    # e.g.  0 0 64   0.00   0.00   0.00  0.73011108E+00 +-  0.15783433E-03
    # e.g.  0 0**   1.00   0.00   0.00    -0.17398298E+00 +-  0.26398018E-04

    # this should match files_out using the "test" program
    gfun_real_space_exp2 = " *[a-zA-z]* ; dx= {2}%s ; dy= {2}%s : *%s *\+- *%s *" % (number_exp, number_exp, sci_not_exp, sci_not_exp)
    # e.g. "point;    dx = 3;    dy = 2: -0.47704462E-14 + -  0.62600131E-16

    # match fourier space functions
    gfun_fourier_space_exp = " *%s *%s *%s *%s *\+- *%s *$" % (number_exp, number_exp, number_exp,
                                                               sci_not_exp, sci_not_exp)
    # match k-vector expression
    kvect_first_exp = " *%s +%s +%s *$" % (number_exp, number_exp, number_exp)
    # e.g. #9        -1.57080  -2.35619
    #
    kvect_next_exp = " *%s +%s *$" % (number_exp, number_exp)

    # parse file
    with open(fname, "r") as fid:

        # read file one section at a time
        section = ["not-empty"]
        section_title = ''
        while section != [] or section_title != '':
            section, section_title = get_file_section(fid)
            section_title = reduce_str(section_title)
            if section == []:
                continue

            #first, determine what kind of section we have

            # now parse section

            # real space functions
            if re.match(gfun_real_space_exp, section[0]):
                data = np.zeros([len(section)])
                unc = np.zeros([len(section)])
                dx = np.zeros([len(section)])
                dy = np.zeros([len(section)])
                for ii, line in enumerate(section):
                    m = re.match(gfun_real_space_exp, line)
                    # dx[ii] = int(float(m.group(4)))
                    # dy[ii] = int(float(m.group(5)))
                    # data[ii] = float(m.group(7)) * 10 ** float(m.group(8))
                    # unc[ii] = float(m.group(9)) * 10 ** float(m.group(10))
                    dx[ii] = int(float(m.group(2)))
                    dy[ii] = int(float(m.group(3)))
                    data[ii] = float(m.group(5)) * 10 ** float(m.group(6))
                    unc[ii] = float(m.group(7)) * 10 ** float(m.group(8))
                meas.update({section_title : equal_time_fn(section_title, data, unc, 1, dx=dx, dy=dy)})
                continue

            # another way real space functions can be defined
            if re.match(gfun_real_space_exp2, section[0]):
                data = np.zeros([len(section)])
                unc = np.zeros([len(section)])
                dx = np.zeros([len(section)])
                dy = np.zeros([len(section)])
                for ii, line in enumerate(section):
                    m = re.match(gfun_real_space_exp2, line)
                    dx[ii] = int(float(m.group(1)))
                    dy[ii] = int(float(m.group(2)))
                    data[ii] = float(m.group(3)) * 10 ** float(m.group(4))
                    unc[ii] = float(m.group(5)) * 10 ** float(m.group(6))
                meas.update({section_title: equal_time_fn(section_title, data, unc, 1, dx=dx, dy=dy)})
                continue

            # fourier space functions
            if re.match(gfun_fourier_space_exp, section[0]):
                data = np.zeros([len(section)])
                unc = np.zeros([len(section)])
                for ii, line in enumerate(section):
                    m = re.match(gfun_fourier_space_exp, line)
                    k_index = int(float(m.group(1)))
                    data[ii] = float(m.group(4)) * 10 ** float(m.group(5))
                    unc[ii] = float(m.group(6)) * 10 ** float(m.group(7))
                meas.update({section_title :equal_time_fn(section_title, data, unc, 1, k_index=k_index)})
                continue

            if re.match("Grid_for_Greens_function", section_title):
                kvects_reps = []
                k_indices = []
                for line in section:
                    m = re.match(kvect_first_exp, line)
                    if m:
                        if k_indices:
                            kvects_reps.append(rep_current)
                        k_indices.append(int(float(m.group(1))))
                        kx = float(m.group(2))
                        ky = float(m.group(3))
                        rep_current = np.array([kx, ky])[None, :]
                        continue
                    m = re.match(kvect_next_exp, line)
                    if m:
                        kx = float(m.group(1))
                        ky = float(m.group(2))
                        rep_current = np.concatenate((rep_current, np.array([kx, ky])[None, :]), 0)
                kvects_reps.append(rep_current) # this otherwise not done at the end
                kvects_reps = [kr for _,kr in sorted(zip(k_indices, kvects_reps))]
                k_indices = sorted(k_indices)
                continue

            if re.match("Grid_for_spin_charge_correlations", section_title):
                continue

            # measurements/settings sections
            #if re.match(settings_exp, section[0]) or re.match(meas_exp, section[0]):
            for line in section:
                # settings parameter
                m = re.match(settings_exp, line)
                if m is not None:
                    var_name = re.sub(r"[:<>*;=()-]", "", m.group(1))
                    var_name = re.sub(r" *$", "", var_name)
                    var_name = re.sub(r" +", "_", var_name)
                    value = float(m.group(2))
                    settings.update({var_name : value})
                    continue

                # measurement result
                m = re.match(meas_exp, line)
                if m is not None:
                    var_name = re.sub(r"[:<>*;=()-]", "", m.group(1))
                    var_name = re.sub(r" *$", "", var_name)
                    var_name = re.sub(r" +", "_", var_name)
                    value = float(m.group(2)) * 10 ** float(m.group(3))
                    err = float(m.group(4)) * 10 ** float(m.group(5))
                    meas.update({var_name : measurement(value, err)})
                    continue
            continue


        return settings, meas, k_indices, kvects_reps

def parse_dqmc_tdm_output_file(fname, match_str="*"):
    """
    Extract data from DQMC time-dependant output file , and return the results as a set of dictionaries which can
     be used to instantiate a qmc_results class object. Using that object, it is easy to average, etc, these
     results. The dictionaries values are instances of the unequal_time_fn class.

    :param fname: path of file to parse
    :param match_str: only time-dependent measurements with names that match this string (using unix wildcard matching
         conventions) will be returned. For example, "*" will return all possible measurements, and "Gfun*" will return
         all Green's function measurements. If you want to e.g. exclude a certain character or group from the match
         in unix style matching, you can write "Gfun*[!k]", which would match any name starting with "Gfun" and not
         ending with "k"
    :return:
    measurements_dict, dictionary of qmc_tdm_measurement class instances
    unparsed_section_list, list of sections that were not parsed
    """

    measurements_dict = {}
    # these just for reference...do not use
    # real space
    #["Gfun", "Gfun_up", "Gfun_dn", "SxSx", "SzSz", "Den-Den", "S-wave", "Conductivity"]
    # k-space
    # ["Gfun", "Gfun_up", "Gfun_dn", "SxSx", "SzSz", "Den-Den", "S-wave", "Gfun_SelfEn", "Gfun_up_SelfEn", "Gfun_dn_SelfEn"]

    # regular expression strings to match various pieces of the output text file
    # some limited testing of these

    # matches a single number
    number_exp = "([+-]?\d+[.]?\d*)"
    # TODO: incorporate possibility of nans
    number_or_nan = number_exp + "|(?i)(nan)"

    # matches a number in scientific notation
    sci_not_exp = number_exp + "E" + number_exp
    # TODO: incorporate possibility of nans
    sci_not_exp_or_nan = sci_not_exp + "|(?i)(nan)"

    # matches the title of a section
    global_title_exp = " *([a-zA-z -]+) *$"

    # matches the title of a section containing a function of real space
    real_space_title_exp = " *([a-zA-z -]+) *%s *(\d) *%s *%s *$" % (number_exp, number_exp, number_exp)
    # e.g. #Gfun         0  0256   1.00   \n
    # matches data line for a function of real space
    real_space_data_exp = " *%s *%s \+ *- *%s *$" % (number_exp, sci_not_exp, sci_not_exp)
    # e.g. #   0.75000  0.20405353E+00 +-  0.15527078E-03\n

    # matches the title of a seciton containing a function of k-space
    k_space_title_exp = " *([a-zA-z -]+) *k *= *%s *pair *= *%s, *%s *$" % (number_exp, number_exp, number_exp)
    # e.g. #Gfun    k = 4    pair = 1, 1
    # matches data line for a function of k-space
    k_space_data_exp = " *%s *\( *%s *\+- *%s *\) *\+i *\( *%s *\+- *%s *\) *$" % (number_exp, sci_not_exp, sci_not_exp, sci_not_exp, sci_not_exp)
    # e.g. #0.00000(0.11134211E+00 + - 0.42901441E-01) + i(0.34694470E-17 + - 0.62530599E-16)

    with open(fname, "r") as fid:

        # for any sections we fail to parse
        unparsed_section_list = []

        # read file one section at a time
        section = ["not-empty"]
        section_title = ''
        while section != [] or section_title != '':
            section, section_title = get_file_section(fid)
            section_title = reduce_str(section_title)
            if section == []:
                continue

            # now parse different types of sections

            # real-space functions which are not spatially resolved (conductivity is the only one)
            if re.match(global_title_exp, section_title):
                m = re.match(global_title_exp, section_title)
                name = reduce_str(m.group(1))

                taus = np.zeros([len(section)])
                data = np.zeros([len(section)])
                unc = np.zeros([len(section)])
                for ii, line in enumerate(section):
                    m = re.match(real_space_data_exp, line)
                    taus[ii] = float(m.group(1))
                    data[ii] = float(m.group(2)) * 10 ** float(m.group(3))
                    unc[ii] = float(m.group(4)) * 10 ** float(m.group(5))
                measurements_dict.update({name: unequal_time_fn(name, taus, data, unc, 1)})
                continue

            # real space functions
            if re.match(real_space_title_exp, section[0]):
                # get name and etc
                m = re.match(real_space_title_exp, section[0])
                # TODO: what are these parameters
                a = int(float(m.group(2)))
                b = int(float(m.group(3)))
                c = int(float(m.group(4)))
                d = int(float(m.group(5)))
                name = reduce_str(m.group(1))
                if not fnmatch.fnmatch(name, match_str):
                    continue
                #print name

                taus = np.zeros([len(section) - 1])
                data = np.zeros([len(section) - 1])
                unc = np.zeros([len(section) - 1])
                for ii, line in enumerate(section[1:]):
                    m = re.match(real_space_data_exp, line)
                    taus[ii] = float(m.group(1))
                    data[ii] = float(m.group(2)) * 10 ** float(m.group(3))
                    unc[ii] = float(m.group(4)) * 10 ** float(m.group(5))
                if name not in measurements_dict.keys():
                    measurements_dict.update({name : unequal_time_fn(name, taus, data, unc, 1)})
                else:
                    measurements_dict[name] = measurements_dict[name] + unequal_time_fn(name, taus, data, unc, 1, dx = d)
                continue

            # fourier space functions
            if re.match(k_space_title_exp, section[0]):
                m = re.match(k_space_title_exp,section[0])
                k_index = int(float(m.group(2)))
                pair_a = int(float(m.group(3)))
                pair_b = int(float(m.group(4)))
                name = "%s_k" % reduce_str(m.group(1))
                if not fnmatch.fnmatch(name, match_str):
                    continue
                #print name, k_index

                taus = np.zeros([len(section) - 1])
                data = np.zeros([len(section) - 1], dtype=complex)
                unc = np.zeros([len(section) - 1], dtype=complex)
                for ii, line in enumerate(section[1:]):
                    m = re.match(k_space_data_exp, line)
                    taus[ii] = float(m.group(1))
                    data[ii] = float(m.group(2)) * 10 ** float(m.group(3)) + 1j * float(m.group(6)) * 10 ** float(m.group(7))
                    unc[ii] = float(m.group(4)) * 10 ** float(m.group(5)) + 1j * float(m.group(8)) * 10 ** float(m.group(9))
                if not name in measurements_dict.keys():
                    measurements_dict.update({name :unequal_time_fn(name, taus, data, unc, 0, k_index=k_index)})
                else:
                    measurements_dict[name] = measurements_dict[name] + unequal_time_fn(name, taus, data, unc, 0, k_index=k_index)
                continue

            #TO DO: missing conductivity because global, not spatially resolved

            unparsed_section_list.append(section)

    return measurements_dict

def parse_geometry_file(fname):
    """
    Parse input geometry file. This contains information about the cluster dimensions and etc.

    :param fname: file name
    :return:
    """

    number_exp = "([+-]?\d+[.]?\d*)"

    geom_dict = {}

    # read file
    with open(fname, "r") as fid:

        # one section at a time
        section = ["not-empty"]
        section_title = ''
        while section != [] or section_title != '':
            section, section_title = get_file_section(fid, section_title_exp='(#.*)')
            section_title = reduce_str(section_title)
            if section == []:
                continue


            if re.match("NDIM.*", section_title):
                geom_dict["NDIM"] = int(section[0])
                continue

            elif re.match("PRIM.*", section_title):
                # this can be a 2x2 or 3x3 matrix. Handle both cases
                ms = []
                for line in section:
                    ms.append(re.match("%s %s ?%s?$" % (number_exp, number_exp, number_exp), line))

                # assuming must be square
                prim = np.zeros([len(ms), len(ms)])
                for ii, m in enumerate(ms):
                    prim[ii, :] = np.array(m.groups()[0:len(ms)])

                geom_dict["PRIM"] = prim
                continue

            elif re.match('SUPER.*', section_title):
                # have only seen these as 2x2 matrices. Probably could be 3x3 in principle?
                ms = []
                for line in section:
                    ms.append(re.match("%s %s ?%s?$" % (number_exp, number_exp, number_exp), line))

                super = np.zeros([len(ms), len(ms)])
                for ii, m in enumerate(ms):
                    super[ii, :] = np.array(m.groups()[0:len(ms)])

                geom_dict["SUPER"] = super
                continue

            elif re.match("ORB.*", section_title):
                # TODO: what to do?
                geom_dict["ORB"] = section
                continue

            elif re.match("HAMILT.*", section_title):
                # last three lines are tup, tdn, and U
                exp = "%s +%s +%s +%s +%s +%s +%s +%s *$" % (number_exp, number_exp, number_exp, number_exp, number_exp, number_exp, number_exp, number_exp)
                m1 = re.match(exp, section[0])
                m2 = re.match(exp, section[1])
                m3 = re.match(exp, section[2])

                hamilt = np.zeros([3, 8])
                hamilt[0, :] = np.array(m1.groups())
                hamilt[1, :] = np.array(m2.groups())
                hamilt[2, :] = np.array(m3.groups())

                geom_dict["HAMILT"] = hamilt
                continue

            elif re.match("SYMM.*", section_title):
                # TODO: what is format and how to handle?
                geom_dict["SYMM"] = section
                continue

            elif re.match("PHASE.*", section_title):
                # TODO: what is format and how to handle?
                geom_dict["PHASE"] = section
                continue

            elif re.match("BONDS.*", section_title):
                # TODO: what is format and how to handle?
                geom_dict["BONDS"] = section
                continue

            elif re.match("PAIR.*", section_title):
                # TODO: what is format and how to handle?
                geom_dict["PAIR"] = section
                continue

            else:
                # only print warning if section non-empty. Otherwise probably picking up on a comment
                if section == []:
                    print("Warning: received unhandled or unexpected section in geometry file: `%s`.", section_title)

    return geom_dict

def parse_input_file(fname):
    """
    Parse QuestQMC input file. This contains settings related to both physical quantities (e.g. chemical potential,
    interaction) and QMC measurement settings (e.g. nbins).

    :param fname: path to the text file to be read
    :return:
    settings_dict: a dictionary of the settings loaded from the input file
    """

    # dictionary to hold results
    settings_dict = {}

    # regular expressions for matching lines
    number_exp = "([+-]?\d+[.]?\d*)"
    numerical_setting = "^ *(.*) *= *%s *$" % (number_exp)
    # don't want to match spaces in settings...
    #other_setting = "^ *(.*) *= *(.*) *$"
    other_setting = "^ *(.*) *= *([^ $]*) *$"

    # read file
    with open(fname, "r") as fid:

        line = "not empty"
        while line != "":

            line = fid.readline()
            # test if line looks like a numerical setting or another type of setting
            if re.match(numerical_setting, line):
                m = re.match(numerical_setting, line)
                var_name = reduce_str(m.group(1))
                value = float(m.group(2))

                settings_dict.update({var_name : value})

            elif re.match(other_setting, line):
                m = re.match(other_setting, line)
                var_name = reduce_str(m.group(1))
                value = m.group(2)

                if var_name != '' and value != '':
                    settings_dict.update({var_name : value})

    return settings_dict

# helper functions for parsing files
def get_file_section(fid, section_title_exp=" *([a-zA-z /'\-\(\)]*):? *([a-zA-z /'\-\(\)]*)?$", section_end_exp="=+"):
    """
    Return section of DQMC file and title. Sections are delimited by strings of equal signs.
    :param fid:
    :param section_title_exp: regular expression that section titles should match. Any line matching this
    expression is interpreted as a section title.
    :param section_end_exp: regular expression that section delimiters should match. Any line matching this
    expression is interpreted as a section delimiter.
    :return:
    section: list of strings, where each is one line of the section.
    section_title: string, title of section
    """

    # section_title_exp = " *([a-zA-z /'\-]*):? *([\(\)a-zA-z /'\-]*)?$"
    # section_end_exp = "=+"

    # maybe want helper function get_section
    section = []
    section_title = ""
    line = fid.readline()

    # keep reading section until encounter section delimiter or next section start
    while not re.match(section_end_exp, line) and line != "":
        m = re.match(section_title_exp, line)

        if not m:
            section.append(line)
        else:
            if section_title == "":
                # TODO: probably should not reduce the string here. Return titles and let caller handle them however.
                #section_title = reduce_str(m.group())
                section_title = m.group()

        file_pointer_pos = fid.tell()
        line = fid.readline()

        # if next line is start of next section, break
        if re.match(section_title_exp, line):
            # rewind file so the next call will pick up the section pointer
            fid.seek(file_pointer_pos)
            break

    return section, section_title

def reduce_str(str):
    """
    Remove extra white space and other character from a string, and replace spaces with underscores. Useful for
    creating a variable name from a string.
    :param str: sring to be reduced
    :return: 
    str_reduced: string after reduction
    """
    # get rid of special characters
    str_reduced = re.sub(r"[:<>*;=()'#]", "", str)
    # trim leading or trailing whitespace
    str_reduced = re.sub(r" *$", "", str_reduced)
    str_reduced = re.sub(r"^ *", "", str_reduced)
    # replace slashes or hyphens with spaces
    str_reduced = re.sub(r"[\-/]", " ", str_reduced)
    # consecutive spaces reduced to single hyphen
    str_reduced = re.sub(r" +", "_", str_reduced)

    return str_reduced

def weighted_avg(values, sigmas):
    """
    based on pschauss fn weighted_average.m
    see https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Standard_error
    https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
    TODO: I think this only handles 1D arrays? Need to check
    :param values: NumPy array of values to be averaged
    :param sigmas: NumPy array of standard deviations associated with values
    :return:
    weighted_mean: the weighted mean of values based on the uncertainties sigma
    weighted_sigma: the uncertainty in the weighted mean
    """
    if np.all(values == values[0]):
        #print("All values equal. Sigma set to mean")
        weighted_mean = values[0]
        weighted_sigma = np.mean(sigmas)
    elif np.any(sigmas == 0):
        print("one or more sigmas was 0. Ignoring sigmas")
        weighted_mean = np.mean(values)
        weighted_sigma = np.std(values)
    else:
        weights = np.divide(1, np.square(sigmas))
        weight_sum = np.sum(weights) #V1
        weight_sqr_sum = np.sum(np.square(weights)) #V2

        weighted_mean = np.sum(values * weights) / weight_sum
        biased_var = np.sum(np.square(values - weighted_mean) * weights) / weight_sum
        # corrector on bottom is analog of (N-1)/N to make this an unbiased estimator
        unbiased_var = biased_var / (1 - (weight_sqr_sum / weight_sum**2))
        weighted_sigma = np.sqrt(unbiased_var)

    return weighted_mean, weighted_sigma

def weighted_avg_array(array_list, sigma_list):
    """
    Weighted average of higher dimensional arrays.
    :param array_list:
    :param sigma_list:
    :return: weighted_mean, weighted_sigma
    """
    ndims = array_list[0].ndim
    shape_single = np.asarray(list(array_list[0].shape))
    stack_size = np.asarray([len(array_list)])

    shape = np.concatenate((shape_single, stack_size), 0)
    data_arr = np.zeros(shape)
    unc_arr = np.zeros(shape)
    for ii,(arr, sigma) in enumerate(zip(array_list, sigma_list)):
        data_arr[..., ii] = arr
        unc_arr[..., ii] = sigma

    last_dim = ndims
    weights = np.divide(1, np.square(unc_arr))
    weight_sum = np.sum(weights, last_dim)
    weights_sqr_sum = np.sum(np.square(weights), last_dim)

    weighted_mean = np.divide(np.sum(np.multiply(data_arr, weights), last_dim), weight_sum)
    weighted_mean_expanded = np.repeat(np.expand_dims(weighted_mean, last_dim), stack_size, last_dim)
    biased_var = np.divide(np.sum(np.multiply(np.square(data_arr - weighted_mean_expanded), weights), last_dim), weight_sum)
    unbiased_var = np.divide(biased_var, 1 - np.divide(weights_sqr_sum, np.square(weight_sum)))
    weighted_sigma = np.sqrt(unbiased_var)

    # TODO: maybe there is a better way of dealing with nans?
    weighted_mean[np.isnan(weighted_mean)] = 0
    weighted_sigma[np.isnan(weighted_sigma)] = 0

    return weighted_mean, weighted_sigma

# classes storing single measurements
class measurement():
    """
    General measurement class, making it easy to work with a value and its uncertainty at the same time.

    """
    def __init__(self, data, unc):
        """

        :param data:
        :param unc:
        """
        self.data = data
        self.unc = unc

    def average(self, inst_stack):
        """
        Average data for a list of measuremnt class instances
        :param inst_stack: a list of measurement instances.
        :return:
        """
        data = [inst.data for inst in inst_stack]
        data.append(self.data)
        unc = [inst.unc for inst in inst_stack]
        unc.append(self.unc)
        data_avg, unc_avg = weighted_avg(np.asarray(data), np.asarray(unc))

        out = measurement(data_avg, unc_avg)
        return out

    def __repr__(self):
        return "%0.6f +/- % 0.6f" % (self.data, self.unc)

    def __add__(self, other):
        if 0:
            #check that all fields except for val and unc are the same between the two...
            raise Exception()
        else:
            out = copy.deepcopy(self)
            # weighted average
            sum_weights = np.divide(1, np.square(self.unc)) + np.divide(1, np.square(other.unc))
            out.data = np.divide(np.divide(self.data, np.square(self.unc)) + np.divide(other.data, np.square(other.unc)), sum_weights)
            # TODO: calculation of uncertainty is not quite right
            out.unc = np.sqrt(np.divide(1, sum_weights))
            return out

    def __radd__(self, other):
        return self.__add__(other)

class equal_time_fn(measurement):
    """
    Class for storing a measurement which has an extra dimension, such as position or k. This is used
    for storing equal-time Green's functions e.g. Inhereits from measurement class.
    """

    def __init__(self, name, data, unc, is_real_space, k_index=0, dx=0, dy=0):
        """

        :param name: name of the equal time function
        :param data: numpy array containing the data of the function
        :param unc: numpy array containing uncertainty for the function
        :param is_real_space: boolean value. 1 if class represents a function of real space. 0 if it represents a function
        of momentum space
        :param k_index: The k-index for a function of momentum space
        :param dx: x-site runner_offset for a function of real space
        :param dy:
        """

        measurement.__init__(self, data, unc)
        self.name = name
        # self.data = data
        # self.unc = unc
        self.is_real_space = is_real_space
        self.k_index = k_index
        self.dx = dx
        self.dy = dy

    def average(self, instance_stack):
        """
        average a list of equal_time_fn instances
        :param instance_stack:
        :return:
        """
        if not instance_stack:
            return copy.deepcopy(self)

        for inst in instance_stack:
            if self != inst:
                raise Exception("At least one element in results stack used different settings. These sets may not be combined.")

        dat_list = [inst.data for inst in instance_stack]
        dat_list.append(self.data)
        unc_list = [inst.data for inst in instance_stack]
        unc_list.append(self.unc)
        # weighted average
        data_avg, unc_avg = weighted_avg_array(dat_list, unc_list)

        out = equal_time_fn(self.name, data_avg, unc_avg, self.is_real_space, self.k_index, self.dx, self.dy)
        return out

    def __eq__(self, other):
        """
        test equality for two equal_time_fn instances
        :param other:
        :return:
        """
        if self.name == other.name and self.is_real_space == other.is_real_space and np.array_equal(self.k_index, other.k_index) and np.array_equal(self.dx, other.dx) and np.array_equal(self.dy, other.dy):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.name

class unequal_time_fn(measurement):
    """
    Class for storing (imaginary) time dependent correlators or Green's functions from QUEST qmc_avg.
    """
    # unequal time function, inherites from measurement class
    def __init__(self, name, tau, data, unc, is_real_space, k_index=np.array([0]), dx=np.array([0]), dy=np.array([0])):
        """
        Class for storing imaginary time dependent Green's functions from QUEST qmc_avg
        :param name: name of the equal time function
        :param tau:
        :param data: numpy array containing the data of the function
        :param unc: numpy array containing uncertainty for the function
        :param is_real_space: boolean value. 1 if class represents a function of real space. 0 if it represents a function
        of momentum space
        :param k_index: The k-index for a function of momentum space
        :param dx: x-site runner_offset for a function of real space
        :param dy:
        """
        # TODO: should I reshape the data here to ensure the same shape?
        if tau.ndim == 1:
            tau = tau[:, None]

        if data.ndim == 1:
            data = data[:, None]

        # uncertainty for complex data types is also treated as a compelx number, where the real part of the uncertainty
        # is the uncertainty of the real part and the imaginary part of the uncertainty is the uncertainty in the
        # imaginary part
        if unc.ndim == 1:
            unc = unc[:, None]

        measurement.__init__(self, data, unc)
        # self.data = data
        # self.unc = unc
        self.name = name
        self.tau = tau
        self.is_real_space = is_real_space
        self.k_index = np.atleast_1d(k_index)
        self.dx = np.atleast_1d(dx)
        self.dy = np.atleast_1d(dy)

    def write_text_file(self, index, fname, delimiter='\t', include_error=0, include_titles=1):
        """
        Write unequal time measurement to text file. First column is imaginary time, second column is real part of
        green's function, and third column is imaginary part.
        :param index: index of the unequal time measurement function to be written to the text file. For k-space functions
        this is the k-vector index.
        :param fname: file name to save text file
        :param delimiter: text file delimiter, defaults to '\t'
        :param include_titles: if 1, include title heading at the top of each data column
        :return:
        """
        with open(fname, 'wb') as f:
            if not include_error:
                if include_titles:
                    f.write(delimiter.join(['tau', 'Re Gfn', 'Im Gfn']) + '\n')
                for ii in range(0, self.data.shape[0]):
                    f.write(delimiter.join(['%0.13f' % self.tau[ii], '%0.13f' % self.data[ii, index].real, '%0.13f' % self.data[ii, index].imag]) + '\n')
            else:
                if include_titles:
                    f.write(delimiter.join(['tau', 'Re Gfn', 'Re Gfn Err', 'Im Gfn', 'Im Gfn Err']) + '\n')
                for ii in range(0, self.data.shape[0]):
                    f.write(delimiter.join(['%0.13f' % self.tau[ii], '%0.13f' % self.data[ii, index].real, '%0.13f' % self.unc[ii, index].real,
                                            '%0.13f' % self.data[ii, index].imag, '%0.13f' % self.unc[ii, index].imag]) + '\n')
        # print("wrote data file %s" % fname)

    def write_err_file(self, index, fname, delimiter='\t', include_titles=1):
        """
        Write unequal time measurement uncertainty to text file. First column is imaginary time, second column is
        uncertainty in the real part of the green's function, and third column is imaginary part.
        :param index: index of the unequal time measurement function to be written to the text file. For k-space functions
        this is the k-vector index.
        :param fname: file name to save text file
        :param delimiter: text file delimiter, defaults to '\t'
        :param include_titles: if 1, include title heading at the top of each data column
        """

        with open(fname, 'wb') as f:
            if include_titles:
                f.write(delimiter.join(['tau', 'Err Re Gfn', 'Err Im Gfn']) + '\n')
            for ii in range(0, self.data.shape[0]):
                f.write(delimiter.join(['%0.13f' % self.tau[ii], '%0.13f' % self.unc[ii, index].real,
                                        '%0.13f' % self.unc[ii, index].imag]) + '\n')
            # print("wrote error file %s" % fname)

    def average(self, instance_stack):
        """
        Average a list of unequal_time_fn measurements
        :param instance_stack:
        :return:
        """
        if not instance_stack:
            return copy.deepcopy(self)

        for inst in instance_stack:
            if self != inst:
                raise Exception(
                    "At least one element in esults stack used different settings. These sets may not be combined.")
        instance_stack.append(self)

        # deal with real and imaginary inputs
        dat_list = [(inst.data.real, inst.data.imag) for inst in instance_stack]
        dat_list_real, dat_list_imag = zip(*dat_list)
        #dat_list.append(self.data)
        unc_list = [(inst.data.real, inst.data.imag) for inst in instance_stack]
        unc_list_real, unc_list_imag = zip(*unc_list)
        #unc_list.append(self.unc)
        # weighted average
        data_avg_real, unc_avg_real = weighted_avg_array(dat_list_real, unc_list_real)
        data_avg_imag, unc_avg_imag = weighted_avg_array(dat_list_imag, unc_list_imag)

        data_avg = data_avg_real + 1j * data_avg_imag
        unc_avg = unc_avg_real + 1j * unc_avg_imag

        out = unequal_time_fn(self.name, self.tau, data_avg, unc_avg, self.is_real_space, self.k_index, self.dx, self.dy)
        return out

    def __add__(self, other):
        """
        Adding equivalent to concatenating...
        :param other:
        :return:
        """
        if not (self.name == other.name and self.is_real_space == other.is_real_space and np.array_equal(self.tau, other.tau)):
            raise Exception("Cannot add unequal_time_fn instances, as either the names, real/k-space-ness, or taus differ")
        data = np.concatenate((self.data, other.data), 1)
        unc = np.concatenate((self.unc, other.unc), 1)
        if self.is_real_space:
            dx = np.concatenate((self.dx, other.dx))
            dy = np.concatenate((self.dy, other.dy))
            k_index = self.k_index
        else:
            k_index = np.concatenate((self.k_index, other.k_index))
            dx = self.dx
            dy = self.dy

        out = unequal_time_fn(self.name, self.tau, data, unc, self.is_real_space, k_index, dx, dy)
        return out

    def __radd__(self, other):
        return self.__add__(self, other)

    def __eq__(self, other):
        """
        Two unequal_time_fn instances evaluate as equal if their data can be averaged together.
        Test two instances of unequal_time_fn for equality by comparing their labels (e.g. name, k-index, dx, dy, etc.),
        but not their data.
        :param other:
        :return:
        """
        if self.name == other.name and self.is_real_space == other.is_real_space and \
                np.array_equal(self.k_index, other.k_index) and np.array_equal(self.dx, other.dx) and\
                np.array_equal(self.dy, other.dy) and np.array_equal(self.tau, other.tau):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.name

# classes storing many measurements		
class qmc_settings():
    """
    class for storing qmc_avg settings. This class should only include values such as U and t, which determine whether two dqmc sets can be averaged together or not. Other quantities, such as random seeds, etc, which might change between QMC runs with the same parameters should not be stored here.
    """
    # TODO: need to add lots of other parameters, like nwrap, etc. etc.

    # Dictionary keys are the field names of this class, and the dictionary values are the names that come from
    # parsing the DQMC text file. These are basically the names that appear in the text file, but with white space,
    # punctuation, and special characters removed.
    # NOTE: some fields may have multiple strings that map to them. This is due to the fact different types of QUST
    # files may use different names...
    measurement_names = {"mu_up" : "mu_up",
                         "mu_dn" : "mu_dn",
                         "t_up" : "t_up",
                         "t_dn" : "t_dn",
                         "U" : "u",
                         "beta" : "beta",
                         "dtau" : "dtau",
                         "Time_slice_L" : "n_time_slices",
                         "L" : "n_time_slices",
                         "Number_of_sites" : "nsites",
                         "nx" : "nx",
                         "ny" : "ny",
                         "Number_of_warmup_sweep" : "nwarmup",
                         "nwarm" : "nwarmup",
                         "Number_of_measurement_sweep" : "nmeasure",
                         "npass" : "nmeasure",
                         "Frequency_of_measurement" : "frq_measurement",
                         "Random_seed" : "seed",
                         "seed" : "seed",
                         "Frequency_of_recomputing_G" : "frq_recomputing_g",
                         "Global_move_number_of_sites" : "global_move_nsites",
                         "Accept_count" : "naccept",
                         "Reject_count" : "nreject",
                         "Approximate_accept_rate" : "accept_rate",
                         "gamma" : "gamma",
                         "NDIM" : "ndim",
                         "PRIM" : "lattice_vectors",
                         "SUPER" : "super_cell_vectors",
                         "ORB" : "orbitals",
                         "HAMILT" : "hamilt",
                         "SYMM" : "symmetry",
                         "PHASE" : "phase",
                         "BONDS" : "bonds",
                         "PAIR" : "pair",
                         "ofile" : "output_file",
                         "gfile" : "geometry_file",
                         "HSF" : "hsf",
                         "bcond" : "bcond",
                         "nbin" : "nbin",
                         "ntry" : "ntry",
                         "tausk" : "tausk",
                         "nhist" : "nhist",
                         "north" : "north",
                         "nwrap" : "nwrap",
                         "fixwrap" : "fixwrap",
                         "errrate" : "err_rate",
                         "difflim" : "diff_lim"

                         }
    # TODO: understand and possibly change ORB, PAIR, HAMILT, etc.

    # fields which should be used to test of sest are taken at the "same" parameters
    comparability_fields = ["mu_up", "mu_dn", "t_up", "t_dn", "u",
                            "beta", "dtau", "n_time_slices", "nsites"]
    # for storing qmc_avg settings -- quantities without uncertainties
    def __init__(self, results_dict, geometry_dict={}, input_dict={}):
        """
        Instantiate class using dictionary of values  obtained from parsing the QMC output file, input file, and input
        geometry file.

        :param results_dict: A dictionary of settings in the format generated by parse_dqmc_output_file().
        The dictionary has the following form: the keys are the names given by parsing the
        QMC output files_out, and the vals are the values associated with these. The dictionary values are mapped
        onto fields of the class, which stores an internal dictionary resolving the names in the QMC file to
        fields of the class.
        """

        # first initialize all fields to NaN's
        for _, field_name in self.measurement_names.items():
            setattr(self, field_name, np.NaN)

        # load settings data from Quest output file dictionary
        for key, val in results_dict.items():

            # check if key maps to any field name of the class
            if key not in self.measurement_names.keys():
                warnings.warn("%s was an input name, but this is not mapped to any field in qmc_static_meas class" % key)
                continue

            else:
                field_name = self.measurement_names[key]

                # check that field_name is already a field name of the class. (it should be because we initialize all
                # the fields we want to NaN's in the beginning
                # due to the way the class is constructed, this check should be redundant
                if field_name in self.__dict__:
                    setattr(self, field_name, val)
                else:
                    warnings.warn("%s is not a field name for qmc_static_meas class." % field_name)

        # load data from geometry input file dictionary
        for key, val in geometry_dict.items():

            # check if key maps to any field name of the class
            if key not in self.measurement_names.keys():
                warnings.warn("%s was an input name, but this is not mapped to any field in qmc_static_meas class" % key)
                continue

            else:
                field_name = self.measurement_names[key]

                # check that field_name is already a field name of the class. (it should be because we initialize all
                # the fields we want to NaN's in the beginning
                # due to the way the class is constructed, this check should be redundant
                if field_name in self.__dict__:
                    attr = getattr(self, field_name)

                    # if field was not initialized before, initialize now
                    if attr == '' or (isinstance(attr, float) and np.isnan(attr)):
                        setattr(self, field_name, val)
                    # if field was initialized, check value is compatible
                    else:
                        if np.any(attr != val):
                            print("WARNING: field %s values loaded from settings file and geometry file are inconsistant." % field_name)

                else:
                    warnings.warn("%s is not a field name for qmc_static_meas class." % field_name)

        # load data from input file dictionary
        # this file has least precedence, as any values obtained from the output or geometry file will be applied
        # instead of these values. Usually this should not be important, but occasionally the input file can have
        # "wrong" parameters, which are ignored by the code on account of values given in the geometry file.
        # As far as I know, the output file should always be authoritative.
        for key, val in input_dict.items():

            # check if key maps to any field name of the class
            if key not in self.measurement_names.keys():
                warnings.warn("%s was an input name, but this is not mapped to any field in qmc_static_meas class" % key)
                continue

            else:
                field_name = self.measurement_names[key]

                # check that field_name is already a field name of the class. (it should be because we initialize all
                # the fields we want to NaN's in the beginning
                # due to the way the class is constructed, this check should be redundant
                if field_name in self.__dict__:
                    attr = getattr(self, field_name)

                    # if field was not initialized before, initialize now
                    if attr == '' or (isinstance(attr, float) and np.isnan(attr)):
                        setattr(self, field_name, val)
                    # if field was initialized, check value is compatible
                    else:
                        if attr != val:
                            print("WARNING: field %s values loaded from settings file and/or geometry file are"
                                  " inconsistant with input settings file. If using geometry file, geometry file"
                                  "will generally take precedence over input file for e.g. chemical potentials, etc." % field_name)

                else:
                    warnings.warn("%s is not a field name for qmc_static_meas class." % field_name)

    @classmethod
    def construct_from_fieldname_dict(cls, fieldname_dict):
        """
        Alternate constructor for creating an instance of this class from a dictionary where the keys are the class
        field names (instead of names derived directly from the DQMC text files_out).
        :param fieldname_dict:
        :return:
        """
        # get reversed dictionary, which has field names as keys and qmc_avg text-file names as values
        rev_dict = {}
        #for key, val in cls.measurement_names.iteritems(): # python2 syntax
        for key, val in cls.measurement_names.items():
            print(key, val)
            rev_dict.update({val: key})

        # get dictionary with qmc_avg text-file names as keys and initial values as values
        qmc_text_file_style_dict = {}
        #for key, val in fieldname_dict.iteritems(): #python2 syntax
        for key, val in fieldname_dict.items():
            qmc_text_file_style_dict.update({rev_dict[key]: val})
        return cls(qmc_text_file_style_dict)

    def is_comparable(self, other):
        """
        Returns 1 if sets can be averaged together, i.e. they have the same U,T, number of sites
        etc. But we will ignore other fields, like naccept, nreject, etc.
        :param other:
        :return:
        """
        for field in self.comparability_fields:
            if self.__getattribute__(field) != other.__getattribute__(field):
                return 0
        return 1

    def __repr__(self):
        # this displayed when you print an instance of this class
        print_str = ""
        for key, val in self:

            print_str = print_str + "%s = %s\n" % (key, repr(val))
        return print_str

    def __eq__(self, other):
        """
        Two qmc_settings instances evaluate as equal if all their attributes are identical
        :param other:
        :return:
        """
        # so we can test different instances for equality
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        # not equal operator. must be defined separately from equality operator
        return not self.__eq__(other)

    def __iter__(self):
        #for key, val in self.__dict__.iteritems():
        for key, val in self.__dict__.items():
            yield key, val

class qmc_meas_base():
    """
    Base class for storing and combining QUESTQMC measurements.
    """

    # Dictionary keys are the field names of this class, and the dictionary values are the names that come from
    # parsing the DQMC text file. These are basically the names that appear in the text file, but with white space,
    # punctuation, and special characters removed.
    measurement_names = { }

    def __init__(self, results_dict):
        """
        Instantiate class from dictionary, where entries in dictionary match the field names
        Each field of this class must be an instance of a class which possess an average method. One point of having
        this class is to average many qmc_avg measurements together. But how to do averaging is not defined in this class
        because it holds different types of measurements. e.g., some are global measurements (e.g. density) and some
        are site resolved (e.g. green's functions). We let each field decide how it is going to be averaged.
        :param results_dict: A dictionary where the keys are strings extracted from a Quest-QMC output file and the values
         are instances of the measurement class

        """

        # first, initialize all fields as NaNs
        for _, field_name in self.measurement_names.items():
            setattr(self, field_name, measurement(np.NaN, np.NaN))

        # set fields based on input dictionary, assuming the input dictionary keys are derived from the dqmc output
        # text file.
        #for key, val in results_dict.iteritems():
        for key, val in results_dict.items():
            if key not in self.measurement_names.keys():
                warnings.warn("%s was an input name, but this is not mapped to any field in qmc_static_meas class" % key)
                continue
            else:
                field_name = self.measurement_names[key]

            if field_name in self.__dict__:
                setattr(self, field_name, val)
            else:
                warnings.warn("%s is not a field name for qmc_static_meas class." % field_name)

    @classmethod
    def construct_from_fieldname_dict(cls, fieldname_dict):
        """
        Alternate constructor for creating an instance of this class from a dictionary where the keys are the class
        field names (instead of names derived directly from the DQMC text files_out).
        :param fieldname_dict:
        :return:
        """
        # get reversed dictionary, which has field names as keys and qmc_avg text-file names as values
        rev_dict = {}
        #for key, val in cls.measurement_names.iteritems(): # python2 syntax
        for key, val in cls.measurement_names.items():
            print(key, val)
            rev_dict.update({val: key})

        # get dictionary with qmc_avg text-file names as keys and initial values as values
        qmc_text_file_style_dict = {}
        #for key, val in fieldname_dict.iteritems(): #python2 syntax
        for key, val in fieldname_dict.items():
            qmc_text_file_style_dict.update({rev_dict[key]: val})
        return cls(qmc_text_file_style_dict)

    @classmethod
    def average(cls, qmc_meas_instance_list):
        dict_combine = {}
        first_inst = qmc_meas_instance_list[0]
        for key, val in first_inst:
            meas_list = [getattr(first_inst, key)]
            # print(key)
            for inst in qmc_meas_instance_list[1:]:
                meas_list.append(getattr(inst, key))
            # let each field do averaging appropriate to its class. This means each field must be a class with an
            # "average" method
            meas = meas_list[0].average(meas_list[1:])
            dict_combine.update({key: meas})
        return cls.construct_from_fieldname_dict(dict_combine)

    def __iter__(self):
        # behavior when looping over instance, returns each measurement object in this class
        #for key, val in self.__dict__.iteritems():
        for key, val in self.__dict__.items():
            yield key, val

    def __repr__(self):
        # this displayed when you print an instance of this class
        # let each field decide how to do this and print that result
        print_str = ""
        for key, val in self:
            print_str = print_str + "%s = %s\n" % (key, repr(val))
        return print_str

class qmc_static_meas(qmc_meas_base):
    """
    Class for storing and combining QUESTQMC measurements of static quantities.
    """
    # TODO: organize this...
    # Dictionary keys are the field names of this class, and the dictionary values are the names that come from
    # parsing the DQMC text file. These are basically the names that appea in the text file, but with white space,
    # punctuation, and special characters removed.
    measurement_names = {"Avg_sign" : "avg_sign",
                         "Avg_up_sign": "avg_up_sign",
                         "Avg_dn_sign" : "avg_dn_sign",
                         "Density" : "n",
                         "Up_spin_occupancy" : "n_up",
                         "Down_spin_occupancy" : "n_dn",
                         "Double_occupancy" : "d",
                         "Kinetic_energy" : "e_kin",
                         "Hopping_energy" : "e_hop",
                         "Potential_energy" : "e_pot",
                         "UN_upN_dn" : "e_int",
                         "Total_energy" : "e",
                         "ZZ_AF_structure_factor" : "zz_af_sf",
                         "Root_Mean_Square_of_ZZ_AF" : "zz_af_sf_rms",
                         "ZZ_Ferro_structure_factor" : "zz_ferro_sf",
                         "XX_AF_structure_factor" : "xx_af_sf",
                         "Root_Mean_Square_of_XX_AF" : "xx_af_sf_rms",
                         "XX_Ferro_structure_factor" : "xx_ferro_sf",
                         "Magnetisation_squared" : "msqr",
                         "Specific_heat" : "c_mu",
                         "Chi_thermal" : "chi_thermal",
                         'Double_Double_correlation_function': "dd",
                         'FT_of_Average_spin_correlation_fn': "ss",
                         'Average_Spin_correlation_function': "ss_k",
                         'FT_of_XX_spin_correlation_fn': "sxsx_k",
                         'FT_of_Double_Double_correlation_fn': "dd_k",
                         'FT_of_ZZ_spin_correlation_fn': "szsz_k",
                         'ZZ_Spin_correlation_function': "szsz",
                         'XX_Spin_correlation_function': "sxsx",
                         'Pairing_correlation_function': "pp",
                         'Single_Single_correlation_function': "nsns",
                         'FT_of_Pairing_correlation_fn': "pp_k",
                         'FT_of_Single_Single_correlation_fn': "nsns_k",
                         "FT_of_Up_Equal_t_Greens_function" : "gfn_up_k",
                         "Down_Equal_time_Greens_function" : "gfn_dn",
                         "Up_Equal_time_Greens_function" : "gfn_up",
                         "Mean_Equal_time_Greens_function" : "gfn",
                         "FT_of_Ave_Equal_t_Greens_function" : "gfn_k",
                         "FT_of_Dn_Equal_t_Greens_function" : "gfn_dn_k",
                         "FT_of_Density_density_correlation_fn_up_up" : "nup_nup_k",
                         "Density_density_correlation_fn_up_dn" : "nup_ndn",
                         "Density_density_correlation_fn_up_up" : "nup_nup",
                         "FT_of_Density_density_correlation_fn_up_dn" : "nup_ndn_k"
                         }

    def __init__(self, results_dict):
        """
        Each field of this class must be an instance of a class which possess an average method. One point of having
        this class is to average many qmc_avg measurements together. But how to do averaging is not defined in this class
        because it holds different types of measurements. e.g., some are global measurements (e.g. density) and some
        are site resolved (e.g. green's functions). We let each field decide how it is going to be averaged.
        """

        # first, instantiate all class members as zeros
        # self.avg_sign = measurement(0, 0)
        # self.avg_up_sign = measurement(0, 0)
        # self.avg_dn_sign = measurement(0, 0)
        # # occupancies
        # self.n = measurement(0, 0)
        # self.n_up = measurement(0, 0)
        # self.n_dn = measurement(0, 0)
        # self.d = measurement(0, 0)
        # # energies
        # self.e_kin = measurement(0, 0) # hopping_e + (mu + U/2)*n
        # self.e_hop = measurement(0, 0)  # <-t \sum_{ij}...>
        # self.e_pot = measurement(0, 0)
        # self.e_int = measurement(0, 0)  # d * U
        # self.e = measurement(0, 0)# e_hop + e_pot
        # # structure factors
        # self.zz_af_sf = measurement(0, 0)
        # self.zz_af_sf_rms = measurement(0, 0)
        # self.zz_ferro_sf = measurement(0, 0)
        # self.xx_af_sf = measurement(0, 0)
        # self.xx_af_sf_rms = measurement(0, 0)
        # self.xx_ferro_sf = measurement(0, 0)
        # # other
        # self.msqr = measurement(0, 0)
        # self.c_mu = measurement(0, 0)
        # self.chi_thermal = measurement(0, 0)
        #
        # # functions of space
        # self.gfn_up_k = measurement(0, 0)
        # self.gfn_up = measurement(0, 0)
        # self.gfn_dn_k = measurement(0, 0)
        # self.gfn_dn = measurement(0, 0)
        # self.gfn = measurement(0, 0)
        # self.gfn_k = measurement(0, 0)
        #
        # self.dd = measurement(0, 0) # <d_i d_j>
        # self.dd_k = measurement(0, 0)
        # self.nsns = measurement(0, 0)
        # self.nsns_k = measurement(0, 0)
        # self.nup_ndn = measurement(0, 0)
        # self.nup_ndn_k = measurement(0, 0)
        # self.nup_nup = measurement(0, 0)
        # self.nup_nup_k = measurement(0, 0)
        #
        # self.sxsx = measurement(0, 0)
        # self.sxsx_k = measurement(0, 0)
        # self.szsz = measurement(0, 0)
        # self.szsz_k = measurement(0, 0)
        # self.ss = measurement(0, 0)  # average spin correlation function
        # self.ss_k = measurement(0, 0)
        #
        # self.pp = measurement(0, 0)
        # self.pp_k = measurement(0, 0)

        # for _, field_name in self.measurement_names:
        #     setattr(self, field_name, measurement(np.NaN, np.NaN))

        # call constructor of base class
        qmc_meas_base.__init__(self, results_dict)

class qmc_tdm_meas(qmc_meas_base):
    """
    # for storing time dependent qmc_avg measurements -- quantities with uncertainties
    """
    measurement_names = {"Gfun" : "gfn",
                         "Gfun_k" : "gfn_k",
                         "Gfun_up" : "gfn_up",
                         "Gfun_up_k" : "gfn_up_k",
                         "Gfun_dn" : "gfn_dn",
                         "Gfun_dn_k" : "gfn_dn_k",
                         "Den_Den": "nn",
                         "Den_Den_k": "nn_k",
                         "SzSz" : "szsz",
                         "SzSz_k" : "szsz_k",
                         "SxSx" : "sxsx",
                         "SxSx_k" : "sxsx_k",
                         "Conductivity" : "conductivity",
                         "S_wave" : "swave",
                         "S_wave_k" : "swave_k",
                         "Gfun_SelfEn_k" : "gfn_self_en_k",
                         "Gfun_up_SelfEn_k" : "gfn_up_self_en_k",
                         "Gfun_dn_SelfEn_k" : "gfn_dn_self_en_k"}

    def __init__(self, tdm_meas_dict):

        # first, instantiate all measurements as zeros
        # self.gfn = measurement(0, 0)
        # self.gfn_k = measurement(0, 0)
        # self.gfn_up = measurement(0, 0)
        # self.gfn_up_k = measurement(0, 0)
        # self.gfn_dn = measurement(0, 0)
        # self.gfn_dn_k = measurement(0, 0)
        #
        # self.nn = measurement(0, 0)
        # self.nn_k = measurement(0, 0)
        #
        # self.szsz = measurement(0, 0)
        # self.szsz_k = measurement(0, 0)
        # self.sxsx = measurement(0, 0)
        # self.sxsx_k = measurement(0, 0)
        #
        # self.conductivity = measurement(0, 0)
        #
        # self.swave = measurement(0, 0)
        # self.swave_k = measurement(0, 0)
        #
        # self.gfn_self_en_k = measurement(0, 0)
        # self.gfn_up_self_en_k = measurement(0, 0)
        # self.gfn_dn_self_en_k = measurement(0, 0)

        qmc_meas_base.__init__(self, tdm_meas_dict)

# class storing all qmc_avg results
class qmc_results():
    """
    store all information from QUESTQMC. This class is a wrapper for holding the various types of possible measurements
    """
    def __init__(self, settings, static_measurements, tdm_measurements, k_indices, kvects_reps):
        """

        :param settings: a qmc_settings class instance
        :param static_measurements: a qmc_static_meas class instance
        :param tdm_measurements: a qmc_tdm_meas class instance
        :param k_indices:
        :param kvects_reps:
        """
        self.settings = settings
        self.static_measurements = static_measurements
        self.tdm_measurements = tdm_measurements
        self.k_indices = k_indices
        self.kvects_reps = kvects_reps # all representations of k_vectors

        self.kvects = np.zeros([len(self.k_indices), 2]) # single convenient rep of each k_vector
        for ii, rep in enumerate(kvects_reps):
            self.kvects[ii, :] = rep[0, :]

        self.navgs = 1

    @classmethod
    def from_dqmc_text_file(cls, fnames, fnames_tdm='', fnames_geom='', fnames_input='', param_match_exp="*"):
        """
        instantiate class directly from DQMC output text files_out
        :param fnames: File name of static measurement output text file. Can also be a list of files.
        :param fnames_tdm: File name of time-dependent of output text file. Can also be a list of files.
        :param param_match_exp: wildcard pattern to filter green's functions contained in the file. This allows you to select
        only a certain type of green's function.
        :return:
        """

        # if fnames is a list, then call this constructor for each element of the list, and return a list of qmc files
        if isinstance(fnames, list):

            # check fnames_tdm input is appropriate
            if fnames_tdm == '':
                fnames_tdm = [''] * len(fnames)

            if not isinstance(fnames_tdm, list) or len(fnames_tdm) != len(fnames):
                raise Exception('fnames_tdm type of length not consistent with fnames.')

            # check fnames_geom is appropriate
            if fnames_geom == '':
                fnames_geom = [''] * len(fnames)

            if not isinstance(fnames_geom, list) or len(fnames_geom) != len(fnames):
                raise Exception('fnames_geom type of length not consistent with fnames.')

            # check fnames_input is appropriate
            if fnames_input == '':
                fnames_input = [''] * len(fnames)

            if not isinstance(fnames_input, list) or len(fnames_input) != len(fnames):
                raise Exception('fnames_geom type of length not consistent with fnames.')

            qmc_list = []
            for f, ftdm, fgeom in zip(fnames, fnames_tdm, fnames_geom):
                qmc_list.append(qmc_results.from_dqmc_text_file(f, fnames_tdm=ftdm, fnames_geom=fnames_geom,
                                                                fnames_input=fnames_input, param_match_exp=param_match_exp))

            return qmc_list

        else:

            # load input file
            if os.path.isfile(fnames_input):
                input_dict = parse_input_file(fnames_input)
            else:
                input_dict = {}
                if fnames_input != '':
                    print("Warning: the input file\n%s\nassociated with the the file\n%s\ndid not exist." % (fnames_input, fnames))

            # load geometry input file
            if os.path.isfile(fnames_geom):
                geom_dict = parse_geometry_file(fnames_geom)
            else:
                geom_dict = {}
                if fnames_geom != '':
                    print("Warning: the geometry input file\n%s\nassociated with the the file\n%s\ndid not exist." % (fnames_geom, fnames))


            # load static results
            if os.path.isfile(fnames):
                settings_dict, meas_dict, k_indices, kvects_reps = parse_dqmc_output_file(fnames)

            else:
                settings_dict = {}
                meas_dict = {}
                k_indices = []
                kvects_reps = []
                # if fnames != '':
                #     print("Warning: the input file %s did not exist." % (fnames))

                # maybe some cases could want to not raise exception? But that will require rewriting parts of the input
                # variable handling that assume fnames always exists and is the reference for other inputs...
                raise Exception("The file %s does not exist." % fnames)

            static_meas = qmc_static_meas(meas_dict)
            settings = qmc_settings(settings_dict, geom_dict, input_dict)

            # load time dependent results
            if os.path.isfile(fnames_tdm):
                meas_tdm_dict = parse_dqmc_tdm_output_file(fnames_tdm, match_str=param_match_exp)
                tdm_meas = qmc_tdm_meas(meas_tdm_dict)
            else:
                tdm_meas = qmc_tdm_meas({})
                if fnames_tdm != '':
                    print("Warning: the time dependent output file\n%s\nassociated with the the file\n%s\ndid not exist." % (fnames_tdm, fnames))



            # instantiate class from results
            qmc = qmc_results(settings, static_meas, tdm_meas, k_indices, kvects_reps)

            return qmc

    @classmethod
    def avg_folder(cls, folder_path, load_tdm_results=True, param_match_exp="*", static_extension=".out", tdm_extension=".tdm.out"):
        """
        instantiate class by averaging a folder containing Quest-QMC output text files_out. Assume that the text files_out
        holding static data and time-dependent data have the same names except for the extension portion

        :param folder_path: path to a folder containing text files_out
        :param load_tdm_results: Whether or not to load time dependent data
        :param param_match_exp: a wildcard string specifying a pattern which determines what Green's functions will be loaded.
        The name of the Green's function in the file must match this string.
        :param static_extension: file extension for static data text files_out
        :param tdm_extension: file extension for time-dependent data text files_out
        :return: qmc_avg
        """

        # get text files_out for static measurements
        dir_expr = os.path.join(folder_path, "*[!%s]%s" % (tdm_extension, static_extension))
        # note that these paths are not sorted.
        files_static_meas = glob.glob(dir_expr)

        # instantiate a stack of qmc_results instances
        qmc_stack = []
        for file_static_meas in files_static_meas:
            folder_path, fname = os.path.split(file_static_meas)
            fname_root,_ = os.path.splitext(fname)
            if load_tdm_results:
                file_tdm_meas = os.path.join(folder_path, fname_root + tdm_extension)
            else:
                file_tdm_meas = ''

            qmc_stack.append(qmc_results.from_dqmc_text_file(file_static_meas, file_tdm_meas, param_match_exp=param_match_exp))

        qmc_avg = qmc_stack[0].average(qmc_stack[1:])
        return qmc_avg

    @classmethod
    def from_folder_path(cls, folder_path, load_tdm_results=True, param_match_exp="*", static_extension=".out", tdm_extension=".tdm.out"):
        """
        instantiate a list of class instances from a folder containing Quest-QMC output text files_out.

        :param folder_path: path to a folder containing text files_out
        :param load_tdm_results: Whether or not to load time dependent data
        :param param_match_exp: a wildcard string specifying a pattern which determines what Green's functions will be loaded.
        The name of the Green's function in the file must match this string.
        :param static_extension: file extension for static data text files_out
        :param tdm_extension: file extension for time-dependent data text files_out
        :return: qmc_list
        """

        # get text files_out for static measurements
        # may there is a smarter way than this?
        dir_expr = os.path.join(folder_path, "*[!%s]%s" % (tdm_extension, static_extension))

        fnames_static_meas = sorted(glob.glob(dir_expr))
        # list for holding file names for tdm measurements
        fnames_tdm_meas = [''] * len(fnames_static_meas)

        if len(fnames_static_meas) > 0:
            # instantiate a stack of qmc_results instances
            qmc_list = [''] * len(fnames_static_meas)
            for ii, file_static_meas in enumerate(fnames_static_meas):
                folder_path, fname = os.path.split(file_static_meas)
                fname_root, _ = os.path.splitext(fname)

                # load tdm results
                if load_tdm_results:
                    fnames_tdm_meas[ii] = os.path.join(folder_path, fname_root + tdm_extension)

                qmc_list[ii] = qmc_results.from_dqmc_text_file(file_static_meas, fnames_tdm_meas[ii], param_match_exp=param_match_exp)
        else:
            qmc_list = ['']
            fnames_static_meas = ['']
            fnames_tdm_meas = ['']

        return qmc_list, fnames_static_meas, fnames_tdm_meas

    def average(self, dqmc_instance_list):
        """
        average results from different dqmc instances
        """
        if dqmc_instance_list == "":
            return copy.deepcopy(self)

        dqmc_instance_list.append(self)
        for ii, inst in enumerate(dqmc_instance_list):

            # compare qmc_avg settings for different sets. Only average qmc_instances if these match
            kvects_eq = np.asarray([np.array_equal(a, b) for a, b in zip(self.kvects_reps, inst.kvects_reps)]).all()
            k_indices_eq = np.array_equal(self.k_indices, inst.k_indices)

            # TODO: maybe need to distinguish between settins like U,T, etc. that really need to be
            # TODO: the same, and others like nwrap that don't
            #if self.settings != inst.settings or not k_indices_eq or not kvects_eq:
        if not self.settings.is_comparable(inst.settings) or not k_indices_eq or not kvects_eq:
                raise Exception("WARNING: Element in qmc_results stack with index % d used different settings. "
                                "Are these sets ok to average together?." % ii)

        # average static measurements using the average method for that class
        static_meas = qmc_static_meas.average([inst.static_measurements for inst in dqmc_instance_list])
        # average dynamic measurements using the average method for that class
        tdm_meas = qmc_tdm_meas.average([inst.tdm_measurements for inst in dqmc_instance_list])
        # static_meas = self.static_measurements.average([inst.static_measurements for inst in dqmc_instance_list])
        # tdm_meas = self.tdm_measurements.average([inst.tdm_measurements for inst in dqmc_instance_list])
        out = qmc_results(self.settings, static_meas, tdm_meas, self.k_indices, self.kvects_reps)
        out.navgs = len(dqmc_instance_list) + 1
        return out

    def __repr__(self):
        print_str = repr(self.settings) + "\n" + repr(self.static_measurements) + "\n" + repr(self.tdm_measurements)
        return print_str

# Misc helper functions
def get_spectral_fn_moments(qmc, kx, ky, spin_state=0):
    """
    Get firs two moments of spectral function for each spin state A_up(k,w) \propto -*Im(G^R_up(k,w)).

    Note that these moments correspond to a spectral function normalized so that \int dw A(k,w) = 1. Whereas, this is
    more usually taken to be 2*pi.
    :param qmc: a qmc_results class instance
    :param kx: value of kx in units of the inverse lattice constant
    :param ky: value of ky in units of the inverse lattice constant
    :param spin_state:
    :return: moments, moments_unc
    """
    qmeas = qmc.static_measurements
    qset = qmc.settings

    #TODO: update for different spin states, for imbalance

    ek = -2 * qset.t_up * (np.cos(kx) + np.cos(ky))
    # this is the same as e_hop I think?
    ke_mod = qmeas.e_kin.data + (qset.mu_dn + qset.u/2) * qmeas.n.data
    ke_mod_unc = np.sqrt(qmeas.e_kin.unc ** 2 + (qset.mu_dn + qset.u/2) ** 2 * qmeas.n.unc ** 2)

    mu_shift = qset.mu_dn + 0.5 * qset.u

    # moment_1 = ek - mu + (n/2)*u
    # mu = U/2 at half-filling in this convention
    # so at half-filling moment_1 = ek
    moment1 = (ek - mu_shift) + qmeas.n.data * qset.u / 2
    moment1_unc = np.abs(qset.u /2 ) * qmeas.n.unc

    # moment_2 = (ek - mu)^2 + U*(ek - mu)*n + 0.5 * u^2 * n
    # at half_filling moment_2 = ek^2 + U^2/4
    moment2 = (ek - mu_shift) ** 2 + qset.u * (ek - mu_shift) * qmeas.n.data + 0.5 * qset.u**2 * qmeas.n.data
    moment2_unc = np.sqrt((qset.u * (ek - mu_shift) * qmeas.n.unc ) **2 ** 2 + (qset.u * mu_shift * qmeas.n.unc) ** 2 + (0.5 * qset.u ** 2 * qmeas.n.unc) ** 2)

    moments = np.array([moment1, moment2])
    moments_unc = np.array([moment1_unc, moment2_unc])

    return moments, moments_unc


if __name__=="__main__":

    #################################################3
    #  load files
    #################################################3
    from convert_qmc_file import *
    # path = "sample_quest_files/ggeom/out"
    path = "C://Users//Peter//Documents//MATLAB//PeterB-analysis//dqmc//sample_quest_files//ggeom"
    path_in = os.path.join(path, "in")
    path_out = os.path.join(path, "out")

    files_out = glob.glob(os.path.join(path_out, '*[!.tdm].out'))
    files_tdm = glob.glob(os.path.join(path_out, '*tdm.out'))
    files_in = glob.glob(os.path.join(path_in, '*.in'))
    files_geom = glob.glob(os.path.join(path_in, '*.geom'))

    #################################################3
    #  static output files test
    #################################################3
    fname = files_out[0]

    # test splitting file into sections
    sections = []
    ttls = []
    with open(fname, "r") as fid:
        section = ["not empty"]
        ttl = ''
        while section != [] or ttl != '':
            section, ttl = get_file_section(fid)
            ttl = reduce_str(ttl)

            sections.append(section)
            ttls.append(ttl)

    # test full file parsing
    settings_dict, meas_dict, k_indices, kvects_reps = parse_dqmc_output_file(fname)
    settings_inst = qmc_settings(settings_dict)
    static_meas_inst = qmc_static_meas(meas_dict)
    #set, meas, ks = parse_dqmc_output_file(fname)

    # qmc_avg loading
    qmc = qmc_results.from_dqmc_text_file(fname)

    #################################################3
    # time dependent ouput files parsing test
    #################################################3
    if len(files_tdm) > 0:
        fname_tdm = files_tdm[0]
        meas_tdm_dict = parse_dqmc_tdm_output_file(fname_tdm, "Gfun*[!k]")
        tdm_meas_inst = qmc_tdm_meas(meas_tdm_dict)

        # this is the best way to average many qmc_avg files_out
        qmc_avg = qmc_results.avg_folder(path_out, "Gfun*[!k]")
        qmc_stack = qmc_results.from_folder_path(path_out, load_tdm_results=1)

        # get moments
        moments = np.zeros([len(qmc_avg.k_indices), 2])
        moments_unc = np.zeros([len(qmc_avg.k_indices), 2])
        for ii, _ in enumerate(qmc_avg.k_indices):
            kx = qmc_avg.kvects_reps[ii][0, 0]
            ky = qmc_avg.kvects_reps[ii][0, 1]
            moments[ii, :], moments_unc[ii, :] = get_spectral_fn_moments(qmc_avg, kx, ky)

    #################################################3
    # geometry file parsing test
    #################################################3
    if len(files_geom) > 0:
        fname_geom = files_geom[0]
        geom_dict = parse_geometry_file(fname_geom)


    #################################################3
    # qmc from all output files
    #################################################3
    if len(files_tdm) > 0 and len(files_geom) > 0:
        qmc_all = qmc_results.from_dqmc_text_file(fname, fname_tdm, fname_geom)


    #################################################3
    # input file parsing test
    #################################################3
    if len(files_in) > 0:
        fname_in = files_in[0]
        input_dict = parse_input_file(fname_in)


    # #################################################
    # test different kinds of output files_out
    # #################################################
    # from convert_qmc_file import *
    # fname1 = "Z://QuestQMC//quest-qmc_avg//EXAMPLE//repulsive_hubbard_dqmc_fitting_results//2018_08_21//out//output7.out"
    # fname2 = "Z:/QuestQMC/quest-qmc_avg/EXAMPLE/UnequalTime/U=+4,period_start=0.1,mu=0/UneqTime.out"
    # fname3 = "Z:/QuestQMC/quest-qmc_avg/EXAMPLE/rfspectroscopy/period_start=0p5,U=0-10,mu=0,npass=20000/outU4/output402.out"
    #
    # settings_dict1, meas_dict1, k_indices1, kvects_reps1 = parse_dqmc_output_file(fname1)
    # settings_inst1 = qmc_settings(settings_dict1)
    # static_meas_inst1 = qmc_static_meas(meas_dict1)
    #
    # settings_dict2, meas_dict2, k_indices2, kvects_reps2 = parse_dqmc_output_file(fname2)
    # settings_inst2 = qmc_settings(settings_dict2)
    # static_meas_inst2 = qmc_static_meas(meas_dict2)
    #
    # settings_dict3, meas_dict3, k_indices3, kvects_reps3 = parse_dqmc_output_file(fname3)
    # settings_inst3 = qmc_settings(settings_dict3)
    # static_meas_inst3 = qmc_static_meas(meas_dict3)

    # from convert_qmc_file import *
    # fname = "Z:/QuestQMC/quest-qmc_avg/EXAMPLE/Conductivity/U=7p5,period_start=0p5,n=0p83,npass=1000x50000/out/output4.tdm.out"
    # tdm_meas_dict, unparsed_sections = parse_dqmc_tdm_output_file(fname)