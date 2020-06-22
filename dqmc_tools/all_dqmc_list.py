import os
import scandir
import glob
import convert_qmc_file as conv_qmc

root_path = os.path.join("Z:", "QuestQMC", "quest-qmc_avg", "EXAMPLE")

settings_list = []
path_list = []
for dir, _, _ in scandir.walk(root_path):
    fnames = glob.glob(os.path.join(dir, '*[!.tdm].out'))
    for fname in fnames:
        settings_dict, meas_dict, k_indices, kvects_reps = conv_qmc.parse_dqmc_output_file(fname)
        settings = conv_qmc.qmc_settings(settings_dict)
        settings_list.append(settings)
        path_list.append(os.path.join(dir, fname))