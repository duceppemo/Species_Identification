import subprocess, os, argparse, pandas, gzip
import seaborn as sns
import matplotlib.pyplot as plt
from concurrent import futures
from multiprocessing import cpu_count


class IdentificationMethods:

    @staticmethod
    def collect_paired(path):
        """

        :param path:
        :return:
        """
        files_dict = dict()
        for f in os.listdir(path):
            name = f.split("_")[0]
            if "R1" in f or "R2" in f:
                if name in files_dict:
                    files_dict[name].append("{}/{}".format(path, f))
                else:
                    files_dict[name] = ["{}/{}".format(path, f)]

        return IdentificationMethods.concatenate_pairs(files_dict)

    @staticmethod
    def concatenate_pairs(paired_dict):
        read_dic = dict()
        for sample_name, lis in paired_dict.items():
            file_path = "/".join(lis[0].split('/')[:-1])
            output_fastq = "{}/{}_all.fastq.gz".format(file_path, sample_name)

            cmd = ["cat", lis[0], lis[1]]
            cat_call = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            stdout, stderr = cat_call.communicate()
            with open(output_fastq, "wb") as ftw:
                ftw.write(stdout)

            read_dic[sample_name] = output_fastq

        return read_dic

    @staticmethod
    def collect_unpaired(path):
        """

        :param path: type STR: a full file path for the directory containing fastq files
        :return: dict
        """
        dirs = [sd for sd in os.listdir(path) if os.path.isdir("{}/{}".format(path, sd))]
        files_dict = dict()
        basename = os.path.basename(path)
        # for directories in dirs:
        #     # requires that this is done on fastq reads from GuppyBaseCaller.py
        #     files_dict[directories] = "{}/{}_{}.fastq".format(path, basename, directories)

        if len(dirs) > 1:
            for directories in dirs:
                # requires that this is done on fastq reads from GuppyBaseCaller.py
                files_dict[directories] = "{}/{}_{}.fastq".format(path, basename, directories)
        else:
            for f in os.listdir(path):
                if os.path.isfile("{}/{}".format(path, f)):
                    file_name = f.split(".")[0]
                    files_dict[file_name] = "{}/{}".format(path, f)

        # IdentificationMethods.mash_screen(ref, files_dict)
        return files_dict

    @staticmethod
    def sample_mash_screen(reference, sample_path, tsv_dict, pair):

        sample_name = os.path.basename(sample_path).split('.')[0]

        output_directory = "{}/mash_screen/".format('/'.join(sample_path.split('/')[:-1]))
        output_table = "{}/mash_screen_{}_identification.tsv".format(output_directory, sample_name)

        cmd = ["mash", "screen", reference, sample_path]
        screen_mash = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = screen_mash.communicate()

        with open(output_table, "w") as ftw:
            ftw.write(stdout.decode('ascii'))

        if os.path.isfile(output_table):
            tsv_dict[sample_name] = output_table
        else:
            print("ERROR\nSomething is not alright with the output! check it out girl")

    @staticmethod
    def mash_screen(reference, dic, pair):
        screen_tsv_dict = dict()

        with futures.ThreadPoolExecutor(max_workers=cpu_count()) as executor:
            paths = (tuple(dic.values()))
            for results in executor.map(
                    lambda x: IdentificationMethods.sample_mash_screen(reference, x, screen_tsv_dict, pair),
                    paths):
                pass

        return screen_tsv_dict

    @staticmethod
    def mash_screen_parse(mash_dict):
        tsv_dict = dict()

        for sample_name, screen_path in mash_dict.items():
            output_tsv = "{}/topscreen_{}_identification.tsv".format('/'.join(screen_path.split('/')[:-1]), sample_name)
            df = pandas.read_csv(screen_path,
                                 sep="\t", usecols=[0, 1, 3, 4],
                                 names=['perc_iden', 'hash_match', 'p_value', 'query_match'])

            df = df.sort_values(by=['perc_iden'], kind='mergesort', ascending=False)

            df = df.head(n=5)
            df = df.reset_index(drop=True)
            df.index = df.index + 1
            df.index.name = 'rank'
            df['query_match'] = list(map(
                lambda x: " ".join(x.split('/')[-1].split("_")[:2]).split(".")[0],
                df['query_match'].tolist()))
            # format the sample names

            df.to_csv(output_tsv, sep="\t")

            if os.path.isfile(output_tsv):
                os.remove(screen_path)
                tsv_dict[sample_name] = output_tsv
            else:
                print("ERROR:\n There was an issue collecting the top 5 hits.")

        return tsv_dict

    @staticmethod
    def mash_visualize_rawreads(sample_dict, tsv_dict, name):
        # will only work with ONE sample

        # the two will share a key, and use this to determine bp per fastq, and the time slot it belongs to
        # build a df with this, and visualize with seaborne

        output_graph_name = "{}_ID_per_coverage.png".format(name)
        output_table = dict()  # convert this to a dataframe at the end

        for sample_name, tsv_path in tsv_dict.items():
            bp_count = 0
            counter = 0
            try:
                with open(sample_dict[sample_name], "r") as ftr:
                    for line in ftr:
                        if counter == 0:
                            counter += 1
                        elif counter == 1:
                            bp_count += len(line.split('\n')[0])

                            counter += 1
                        elif counter == 2:
                            counter += 1
                        elif counter == 3:
                            counter = 0

            except UnicodeDecodeError:
                with gzip.open(sample_dict[sample_name], "r") as ftr:
                    for line in ftr:
                        if counter == 0:
                            counter += 1
                        elif counter == 1:
                            bp_count += len(line.decode('ascii').split('\n')[0])
                            counter += 1
                        elif counter == 2:
                            counter += 1
                        elif counter == 3:
                            counter = 0

            coverage = round(bp_count / 100000)
            tsv_df = pandas.read_csv(tsv_path, sep='\t', usecols=['hash_match', 'query_match'])  # get query information
            hash_match = tsv_df['hash_match'].tolist()
            query_match = tsv_df['query_match'].tolist()
            queries = dict()
            index = 0
            while index < 5:  # base on how many samples we want
                queries[os.path.basename(query_match[index].split('.')[0])] = hash_match[index]
                index += 1

            output_table[coverage] = queries

        df = pandas.DataFrame(output_table)
        col_names = ['query_names'] + df.columns.tolist()
        df = df.reset_index()
        df.columns = col_names
        df = pandas.melt(df, id_vars=['query_names'], value_vars=list(output_table.keys()))
        df.columns = ['query names', 'hundred kb sequenced', 'percentage of shared bins']
        df['hundred kb sequenced'] = df['hundred kb sequenced'].astype(float)
        # remove_nan = [0 for f in df['shared bins'].tolist() if str(f) == True]
        df['percentage of shared bins'] = df['percentage of shared bins'].fillna(0)
        df['percentage of shared bins'] = list(
            map(
                lambda x: (int(x.split('/')[0]) / (int(x.split('/')[-1]))) * 100 if isinstance(x, str) else x,
                df['percentage of shared bins'].tolist()
            )
        )
        df['percentage of shared bins'] = df['percentage of shared bins'].astype(float)

        line_plt = sns.lineplot(x='hundred kb sequenced', y='percentage of shared bins', hue='query names', data=df,
                                legend='full')
        box = line_plt.get_position()
        line_plt.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        hnd, lab = line_plt.get_legend_handles_labels()
        lgd = plt.legend(handles=hnd, labels=lab, bbox_to_anchor=(1.05, 1), loc=0, borderaxespad=0, fontsize='x-small',
                         title='Species')

        # line_plt.legend(handles=hnd, labels=lab, loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

        line_plt.get_figure().savefig(output_graph_name, bbox_extra_artists=(lgd,), bbox_inches='tight')

    @staticmethod
    def clear_pair(dic):
        for path in dic.values():
            print("Removing {} \n".format(path))
            os.remove(path)

    @staticmethod
    def mash_triangle(screen_dic, ref):
        all_files = list(screen_dic.values())
        cmd = ["mash", "triangle", ref] + all_files

        file_path = "/".join(all_files[0].split('/')[:-1])
        output = "{}/mash_screen/mash_triangle/distance_matrix.tsv".format(file_path)

        triangle_call = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, sdterr = triangle_call.communicate()

        with open(output, "wb") as ftw:
            ftw.write(stdout)

        return output

    @staticmethod
    def parse_distancematrix(path):
        """
        Parse through the tsv by treating it as a pandas dataframe
        :param path:
        :return:
        """
        path_dict = '/'.join(path.split('/')[:-1])
        tsv_df = pandas.read_csv(path, sep='\t', index=0, header=0)

    # @staticmethod
    # def collect_all(path, typ, ref):
    #     directory_count = [f for f in os.listdir(path) if os.path.isdir("{}/{}".format(path, f))]
    #     if len(directory_count) > 1:
    #         for d in directory_count:
    #             if typ == 'raw-reads':
    #                 IdentificationMethods.collect_samples("{}/{}".format(path, d))
    #             elif typ == 'assemblies':
    #                 IdentificationMethods.collect_assemblies("{}/{}".format(path, d))
    #     else:
    #         if typ == 'raw-reads':
    #             IdentificationMethods.collect_samples("{}/{}".format(path), ref)
    #         elif typ == 'assemblies':
    #             IdentificationMethods.collect_assemblies("{}/{}".format(path), ref)
    #
    # @staticmethod
    # def collect_assemblies(path, ref):  # ADJUST LATER FOR REF / NEW ORGANIZATION
    #     directories = os.listdir(path)
    #     all_assemblies = dict()
    #     for sd in directories:  # should be each sample folder
    #         msh_path = IdentificationMethods.find_assembly("{}/{}".format(path, sd))
    #         all_assemblies[sd] = msh_path
    #
    #     return all_assemblies
    #
    # @staticmethod
    # def find_assembly(directory):
    #     # finds which assembly or polished fasta needs to be given to sketch assemblies
    #     sample_name = os.path.basename(directory)
    #     directories = [sd for sd in os.listdir(directory) if os.path.isdir("{}/{}".format(directory, sd))]
    #
    #     if "pilon" in directories:
    #         final_pilon_int = max(set(map(lambda x: x.split(".")[0].split("_")[-1],
    #                                       os.listdir("{}/pilon".format(directory)))))
    #
    #         mash_input = "{}/pilon/{}_polished_{}.fasta".format(directory, sample_name, final_pilon_int)
    #         output = IdentificationMethods.sketch_assemblies(mash_input)
    #
    #     elif "medaka_consensus" in directories:
    #         mash_input = "{}/medaka_consensus/consensus.fasta".format(directory)
    #         output = IdentificationMethods.sketch_assemblies(mash_input)
    #
    #     elif "spades" in directories:
    #         max_kmer = max(list(map(lambda x: int(x.split("K")[-1]),
    #                                 filter(lambda x: "K" in x,
    #                                        os.listdir("{}/spades".format("{}/spades/".format(directory)))))))
    #
    #         mash_input = "{}/spades/K{}/final_contigs.fasta".format(directory, max_kmer)
    #         output = IdentificationMethods.sketch_assemblies(mash_input)
    #
    #     elif "ShastaRun" in directories:
    #         mash_input = "{}/ShastaRun/Assembly.fasta"
    #         output = IdentificationMethods.sketch_assemblies(mash_input)
    #
    #     assert output, "ERROR: Provided path did not have any assemblies that could be identified"
    #
    #     return output

    # @staticmethod
    # def sketch_assemblies(file):
    #     parent_dir = "/".join(file.split("/")[:-2])
    #     sample_name = "/".join(file.split("/")[-2])
    #     mash_no_ext = "{}/mash/{}_sketch".format(parent_dir, sample_name)
    #     output_mash = mash_no_ext + ".msh"
    #     output_log = mash_no_ext + ".log"
    #
    #     cmd = []
    #     sketching = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #     stdout, stderr = sketching.communicate()
    #
    #     with open(output_log, "wb") as ftw:
    #         ftw.write(stdout)
    #
    #     return output_mash

    #
    # @staticmethod  # Won't be necessary for screen
    # def sketch_samples(file):
    #     """
    #     For fastq.gz. Takes a single file,
    #     :param file:
    #     :return:
    #     """
    #
    #     parent_dir = "/".join(file.split("/")[:-1])
    #     sample_name = os.path.basename(file).split(".")[0]
    #     mash_name_noext = "{}/mash/{}_sketch".format(parent_dir, sample_name)
    #     output_msh = "{}.msh".format(mash_name_noext)
    #     output_log = mash_name_noext + ".log"
    #
    #     cmd = ["mash", "sketch",
    #            "-I", sample_name,
    #            "-r",
    #            "-m", "2",
    #            "-s", "10000",
    #            "-o", mash_name_noext,
    #            file]
    #
    #     sketching = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #     stdout, stderr = sketching.communicate()
    #
    #     with open(output_log, "wb") as ftw:
    #         ftw.write(stdout)
    #     # write the stdout and the sterr into two different log files
    #     # log files can go ... ? in the mash folder
    #
    #     # this might be silly because a sketch should be make!
    #     if os.path.isfile(output_msh):
    #         run_info = output_msh
    #     else:
    #         run_info = False
    #
    #     return run_info
    #
    # @staticmethod  # not necessary with mash
    # def sketch_all(dic):
    #     sketches = dict()
    #     for sample_name, file in dic.items():
    #         if "fasta" in sample_name or "fa" in sample_name:
    #             sketches[sample_name] = IdentificationMethods.sketch_assemblies(file)
    #         else:
    #             sketches[sample_name] = IdentificationMethods.sketch_samples(file)
    #
    #     return sketches
    #
    # @staticmethod
    # def sketch_dist(dic, reference):
    #     """
    #     Create a directory for each sample in the intial path, and
    #     :param dic:
    #     :param reference:
    #     :return:
    #     """
    #
    #     dist_tsv_dict = dict()
    #
    #     for samples_names, sample_paths in dic.items():
    #         output_directory = '/'.join(sample_paths.split('/')[:-1])
    #         output_table = "{}/{}_dist.tsv".format(output_directory, samples_names)
    #
    #         cmd = ["mash", "dist",
    #                "-p", "1",
    #                "-s", "10000",
    #                reference,
    #                sample_paths]  # add output table and inputs and reference
    #         distance_mash = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #         stdout, stderr = distance_mash.communicate()
    #
    #         # can I make this into a pandas dataframe? something marco wants to do
    #         with open(output_table, "w") as ftw:
    #             ftw.write(stdout.decode('ascii'))
    #
    #         if os.path.isfile(output_table):
    #             dist_tsv_dict[samples_names] = output_table
    #         else:
    #             print("ERROR")
    #
    #     # os.remove all mash tables and empty log files - only keep the diagnostic!
    #
    #     return dist_tsv_dict


class Identification(object):
    def run(self):
        if self.type.lower() in ["assembly", "raw-reads"]:
            if not os.path.isdir("{}/mash_screen/".format(self.path)):
                os.mkdir("{}/mash_screen/".format(self.path))
            if self.type == "assembly":
                start = IdentificationMethods.collect_assemblies(self.path)
            else:
                if not self.pair:
                    start = IdentificationMethods.collect_unpaired(self.path)
                else:
                    start = IdentificationMethods.collect_paired(self.path)

        # start = IdentificationMethods.collect_all(self.path, self.type)
        screen_samples = IdentificationMethods.mash_screen(self.ref, start, self.pair)
        final_tsv = IdentificationMethods.mash_screen_parse(screen_samples)
        # IdentificationMethods.mash_visualize_rawreads(start, final_tsv, self.sample_name)
        if self.tree:
            if not os.path.isdir("{}/mash_screen/mash_triangle".format(self.path)):
                os.mkdir("{}/mash_screen/mash_triangle".format(self.path))

            distance_matrix = IdentificationMethods.mash_triangle(start, self.ref)
            IdentificationMethods.parse_distancematrix(distance_matrix)

        if self.clear:
            if self.pair:
                IdentificationMethods.clear_pair(start)
            else:
                print("Careful! You're trying to remove your files without having used Illumina reads. \n ",
                      "We left things were they are")

        # all_sketch = IdentificationMethods.sketch_all(start)
        # tsv_dict = IdentificationMethods.sketch_dist(all_sketch, self.ref)
        # IdentificationMethods.mash_output_parse(tsv_dict)

    def __init__(self, args):

        # Command line arguments
        self.path = args.path
        self.type = args.type
        self.ref = args.ref
        self.pair = args.pair
        self.tree = args.tree
        self.clear = args.clear

        # RUN
        self.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Allows for the identification of raw reads or assemblies via mash")
    parser.add_argument("-d",
                        help="Input directory with all assembly or read directories.", dest="path",
                        required=True, metavar="/input/folder/")
    parser.add_argument("-ref",
                        help="Custom reference sketch file. Default is for Mycobacterium", required=False,
                        metavar="/path/tp/file.msh", dest="ref")  # set default
    parser.add_argument("-type",
                        help="Type of input, either assembly or raw reads. Raw reads must be in file for each sample.",
                        metavar="[assemblies OR raw-reads]", required=True, dest="type")
    # parser.add_argument("-name",
    #                     help="REMOVE AFTER. NOT A SMART THING",
    #                     metavar="full file path", dest="sample_name")

    parser.add_argument("--paired", action='store_true',
                        default=False, dest="pair", help="Provided directory is for paired end Illumina reads.")
    #
    parser.add_argument("--tree", action='store_true',
                        default=False, dest="tree", help="Outputs a tree, with all sequences in input and references.")

    parser.add_argument("--clear", action='store_true',
                        default=False, dest="clear",
                        help="Will remove the concatenated Illumina file, if Illumina reads were used")
    args = parser.parse_args()
    Identification(args)
