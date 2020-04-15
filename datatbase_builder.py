import argparse, os, subprocess
from concurrent import futures
from multiprocessing import cpu_count

class FastaDir:
    def __init__(self, args):
        self.path = args.path
        # Create the output directory if it does not already exist
        if not os.path.isdir("{}/mash_sketches/".format(self.path)):
            os.mkdir("{}/mash_sketches/".format(self.path))
            subprocess.run(["chmod", "777", "{}/mash_sketches".format(self.path)])  # previously ran into permission issues

        self.run()

    def run(self):
        spec_dict = DatabaseBuilding.unique_species(self.path)
        DatabaseBuilding.paste_it(spec_dict, self.path)


class DatabaseBuilding:
    
    @staticmethod
    def delete_files(path, lis):
        delete_path = list(map(lambda x: path + "/" + x, lis))
    
        for f in delete_path:
            print("Say bye to {}\n".format(f))
            os.remove(f)

    @staticmethod
    def unique_species(path):
        # finds all unique species and selects the first of each to be the chosen representative
        spec_dic = dict()
        files = [f for f in os.listdir(path) if os.path.isfile("{}/{}".format(path, f))] # are only allowed to have fasta files in this directory
        # MTC = ["bovis", "orygis", "avium", "microti", "canetti", "caprae", "africanum", "mungi", "pinnipedii", "suricattae", "tuberculosis"]
        # haven't run into this issue of odd naming for the other MTC 

        for f in files:
            try:
                name = f.split('_')[1].lower()
                if "bovis" in f and "bovis" not in name: # find a way to check this with all the MTC, not just bovis
                    name = f.split("_")[3].lower()
            except IndexError:
                name = f.split('-')[1].lower()
                # accomodates naming of bovis, might need to be extended for other MTC
                if "bovis" in f and "bovis" not in name:
                    name = f.split("-")[3].lower()
            if "." in name:
                name = name.split(".")[0]

            if name in spec_dic or "sp" == name:
                pass
            else:
                spec_dic[name] = "{}/{}".format(path, f)

            if len(spec_dic) == len(files):
                print("No duplicates! You have {} unique species\n".format(len(spec_dic)))

        DatabaseBuilding.acc_finder(spec_dic, path)
        return DatabaseBuilding.sketch(spec_dic)


    @staticmethod
    def acc_finder(dic, path):

        output_txt = "{}/reference_accession_species.txt".format(path)

        with open(output_txt, "w") as ftw:
            for name, path in dic.items():
                with open(path, "r") as ftr:
                    acc_num = ftr.readline().split(" ")[0].split(">")[-1]
                    ftw.write("{} \t {} \n".format(name, acc_num))

        print("All the accession numbers in your reference are recorded!\n")

    @staticmethod
    def mash_sketch(path, name, dic):

        home_dir = '/'.join(path.split("/")[:-1])

        base = os.path.basename(path).split(".")[0]
        mash_out = "{}_sketch.msh".format(base)

        cmd = ["mash", "sketch", "-s", "10000",
               "-o", "{}/mash_sketches/{}".format(home_dir, mash_out),
               path]  # subprocess call does not re

        subprocess.run(cmd)
        dic[name] = "{}/mash_sketches/{}".format(home_dir, mash_out)

    @staticmethod
    def sketch(dic):
        # allows for the parallelization of calling mash sketch 
        sketch_loc = dict()
        with futures.ThreadPoolExecutor(max_workers=cpu_count()) as executor:
            paths = (tuple(dic.items()))
            for results in executor.map(
                    lambda x: DatabaseBuilding.mash_sketch(x[1], x[0], sketch_loc), paths):
                        pass

        return sketch_loc

    @staticmethod
    def paste_it(dic, path):
        print("Pasting all the sketches to get one reference file\n!")
        all_sketches = dic.values()
        final_output = "{}/Mycobacterium_reference".format(path)
        with open("subprocc_input.txt", "w") as ftw:
            for sketch in all_sketches:
                if "reference" in sketch:
                    pass
                else:
                    ftw.write("{}\n".format(sketch))

        cmd = ["mash", "paste",
               final_output, "-l",
               "subprocc_input.txt"]

        subprocess.run(cmd)
        os.remove("subprocc_input.txt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identifies all unique MTC subspecies and species and deletes all other files./nIntended to be used on a directory with symlinks. All unqiue files are used to build a mash database!\nFinal output is a .msh file in the input directory."")
    parser.add_argument("-d", metavar="PATH", dest="path", help="Provide the directory with the symlinks to the reference assemblies.)

    args = parser.parse_args()

    FastaDir(args)


