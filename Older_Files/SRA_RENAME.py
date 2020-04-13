import pandas, argparse,  os

"""
PURPOSE: Rename SRA files based on species 
"""


class InputArgs:
    def __init__(self, args):
        RenameFun.rename_files(args.name, args.dir)
    pass


class RenameFun:
    @staticmethod
    def rename_files(file, directory):
        directory_files = os.listdir(directory)
        file_paths = list(map (
            lambda x: directory + "/" + x, directory_files
        ))

        with open(file, 'r') as ftr:
            li = ftr.readlines()
            for l in li:
                old_name = l.split('\t')[0]
                for file in file_paths:
                    if old_name in file:
                        print(old_name, file)
                        new_name = l.split("\t")[1].split()[0]
                        new_file_name = "{}/{}.fasta".format(directory, new_name)
                        os.rename(file,new_file_name )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="rip")
    parser.add_argument("-txt", dest="name", help="file with names to replace and current file name")
    parser.add_argument("-d", dest="dir", help="directory with files")

    args = parser.parse_args()

    InputArgs(args)
