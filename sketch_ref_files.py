import os, argparse, subprocess


class InputDir:
    def __init__(self, args):
        self.dirpath = args.dirpath  # Path to directory with all reference files to be used for the database
        self.final = args.final  # final name of the output sketches

        mash_dir = '/'.join(self.dirpath.split('/')[:-1])  # directory where mash sketches for building are left

        if not os.path.isdir("{}/mash_sketches/".format(mash_dir)):
            os.mkdir("{}/mash_sketches/".format(mash_dir))

        self.run()

        subprocess.run(["chmod", "777", "{}/mash_sketches".format(mash_dir)])  # ran into permission issues
        File.concatenate("{}/mash_sketches".format(mash_dir), self.final)

    def run(self):
        files = os.listdir(self.dirpath)

        for f in files:  # would be better to parallelize
            fi = "{}/{}".format(self.dirpath, f)
            File.mash_sketch(fi)


class File:

    @staticmethod
    def mash_sketch(path):
        home_dir = '/'.join(path.split("/")[:-2])
        base = os.path.basename(path).split(".")[0]
        mash_out = "{}_sketch.msh".format(base)

        cmd = ["mash", "sketch", "-s", "10000",
               "-o", "{}/mash_sketches/{}".format(home_dir, mash_out),
               path]  # subprocess call does not re

        subprocess.run(cmd)

    @staticmethod
    def concatenate(path, final):
        all_sketches = list(map(lambda x: "{}/{}".format(path,x), os.listdir(path)))
        with open("subprocc_input.txt", "w") as ftw:
            for sketches in all_sketches:
                ftw.write("{}\n".format(sketches))

        cmd = ["mash", "paste",
               final, "-l",
               "subprocc_input.txt"]

        subprocess.run(cmd)
        os.remove("subprocc_input.txt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a sketch database for mash!")
    parser.add_argument("-d",
                        help="Path to directory with all reference sequences", dest="dirpath",
                        required=True, metavar="/input/folder/")
    parser.add_argument("-final",
                        help="Name of the output prefix for the concatenated sketches",
                        required=True,
                        metavar="str", dest="final")  # set default

    arg = parser.parse_args()
    InputDir(arg)
