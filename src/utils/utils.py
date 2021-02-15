import cProfile
import io
import logging
import os
import pathlib
import pstats
import sys
import numpy as np
import yaml
import pandas as pd
from pprint import pprint


def read_pandas(file_path, full=True, nrows=1000):
    if full is True:
        return pd.read_csv(file_path, compression="gzip", sep="\t", low_memory=False)
    if full is False:
        return pd.read_csv(file_path, compression="gzip", sep="\t", low_memory=False, nrows=nrows)


def load_config_file(config_file="src/config.yaml"):
    return yaml.load(open(config_file), Loader=yaml.FullLoader)


def add_files_yaml_config(implicated_class, variable_name, file_path, config_file):
    config_file_data = yaml.load(open(config_file), Loader=yaml.FullLoader)
    if implicated_class not in config_file_data:
        config_file_data[implicated_class] = dict()
    config_file_data[implicated_class][variable_name] = file_path
    yaml.dump(config_file_data, open(config_file, "w"))
    print("Add {} : {} in {} part to config file {}".format(variable_name, file_path, implicated_class, config_file))


def remove_files_yaml_config(implicated_class, variable_name, file_path, config_file):
    implicated_class = implicated_class.__class__.__name__
    config_file_data = yaml.load(open(config_file), Loader=yaml.FullLoader)
    config_file_data[implicated_class].pop(variable_name, None)
    yaml.dump(config_file_data, open(config_file, "w"))
    print("Removed {} : {} in {} part to config file {}".format(variable_name, file_path, implicated_class, config_file))


def mkdir(init_directory):
    """

    Args:
            directory: str
                    Directory to create

    Returns:

    """
    dirs = list(sorted(list(set(["/".join(init_directory.split("/")[: j + 1]) for j in range(len(init_directory.split("/")))]))))
    for directory in dirs:
        if not os.path.exists(directory):
            try:
                pathlib.Path(directory).mkdir(exist_ok=True)
            except FileNotFoundError:
                logging.error("Unable to find or create directory {}".format(directory))
                sys.exit("============\nSee you soon :)\n============")

def mkdirs(output_dir, subdirs_output_dir):
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
		[os.makedirs(output_dir + e) for e in subdirs_output_dir]

# ADD HGNC COLUMN
def add_hgnc(row):
    row = [e.replace("HGNC:HGNC:", "") for e in row.split(",") if "HGNC" in e]
    if row:
        return row[0]
    else:
        return None


def profile(fnc):
    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner


def stats(list_values, elem):
    dict_stats = dict()
    dict_stats["nb"] = len(list_values)
    tmp_list_size = list()
    for e in list_values:
        tmp_list_size.append(int(e[1]) - int(e[0]))
    if tmp_list_size:
        dict_basic_stats = basic_stats(tmp_list_size)
        dict_stats = {**dict_stats, **dict_basic_stats}
        dict_stats["Element"] = elem
    return dict_stats


def basic_stats(tmp_list):
    dict_stats = dict()
    dict_stats["mean"] = np.mean(tmp_list)
    dict_stats["max"] = np.max(tmp_list)
    dict_stats["min"] = np.min(tmp_list)
    dict_stats["median"] = np.median(tmp_list)
    dict_stats["std"] = np.std(tmp_list)
    dict_stats["total_sum"] = np.sum(tmp_list)
    dict_stats["nb"] = len(tmp_list)
    return dict_stats


def update_dict_files(d, dir):
    d = {k: dir + e for k, e in d.items() if dir not in e}
    return d


def create_new_sub_dir_and_update_dict(old_dir, new_output_dir, dict_files):
    if old_dir.endswith("/") is False:
        old_dir += "/"
    if new_output_dir.endswith("/") is False:
        new_output_dir += "/"
    old_dir = old_dir.split("/")
    return_dir = old_dir.copy()
    return_dir[-1] = new_output_dir
    return_dir = "/".join(return_dir)
    mkdir(return_dir)
    return update_dict_files(dict_files, return_dir)


def setup_custom_logger(name):
    formatter = logging.Formatter(fmt="%(asctime)s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    # mkdir('Logging')
    handler = logging.FileHandler(name)
    handler.setFormatter(formatter)
    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.addHandler(screen_handler)
    return logger


#
def prepare_header(vcf, vep_field, vep_separator):

    index_dict = dict()
    if vep_field:
        for h in vcf.header_iter():
            try:
                if h.info()["ID"] == vep_field:
                    csq_header = h.info()["Description"].split(vep_separator)
                    for elem in csq_header:
                        index_dict[elem] = csq_header.index(elem)
            except:
                pass
    return index_dict


# if __name__ == '__main__':
# 	from tqdm import  tqdm
# 	tqdm.pandas()
# 	df = pd.read_csv("../../data/2_processed/GRCh37_RefSeq_lite.csv.gz", compression='gzip', sep='\t')
# 	df['HGNC'] = df['Dbxref'].progress_apply(lambda r: add_hgnc(r))
# 	df.to_csv("../../data/2_processed/GRCh37_RefSeq_lite_hgnc.csv.gz", compression='gzip', sep='\t', index=False)
#
