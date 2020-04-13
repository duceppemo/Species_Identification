import os, argparse


class File:

    def __init__(self, arg):
        self.ref = args.ref
        self.out = arg.out
        self.path = arg.path

        self.makesym()

    def makesym(self):

        if not os.path.isdir(self.out):
            os.mkdir(self.out)

        with open(self.ref, "r") as files:
            species_names = files.readlines()
            species_names = list(map(lambda x: x.strip(), species_names))
            file_path = list(map(lambda x: "{}/{}".format(self.path, x), species_names))

            for f in file_path:
                try:
                    os.symlink(f, "{}/{}".format(self.out, f.split("/")[-1]))
                except FileExistsError:
                    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="makes symlinks to files in path to output")
    parser.add_argument("-names",
                        help="Input directory with all assembly or read directories.", dest="ref",
                        required=True, metavar="/input/folder/file")
    parser.add_argument("-out",
                        help="Custom reference sketch file. Default is for Mycobacterium", required=True,
                        metavar="/path/tp/", dest="out") # set default
    parser.add_argument("-path",
                        help="where is files",
                        required=True, dest="path")

    args = parser.parse_args()
    File(args)
