import yaml
import shutil
import argparse
from tqdm import tqdm
import os


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Copies mutants' short reads from Argos",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-y", "--yaml", metavar="<yaml file with parent and its mutants>",
                        help="yaml file with parent and its mutants")
    parser.add_argument("-o", "--output", metavar="<output directory>",
                        help="output directory: resources/mutants/<parental_strain>")

    return parser.parse_args()


if __name__ == '__main__':
    path_argos = "/home/andrei/Data/Argos/imb_sal_raw/Raw_sequencing_data/Sheida_Heidarian/F20FTSEUHT0780-02_BACzutR/Clean"
    args = get_args()
    # read yaml config
    with open(args.yaml) as f:
        config = yaml.full_load(f)
        mutants = config["mutant"].keys()

    # copy read files for each strain
    for mutant in tqdm(mutants):
        os.makedirs(os.path.join(args.output, mutant), exist_ok=True)
        r1_file_name = mutant + "_1.fq.gz"
        r2_file_name = mutant + "_2.fq.gz"
        shutil.copyfile(src=os.path.join(path_argos, mutant, r1_file_name), dst=os.path.join(args.output, mutant, r1_file_name))
        shutil.copyfile(src=os.path.join(path_argos, mutant, r2_file_name), dst=os.path.join(args.output, mutant, r2_file_name))
    print("Finished")
