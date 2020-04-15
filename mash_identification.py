import subprocess, argparse, pandas, gzip
import os
import numpy as np
import seaborn as sns
from concurrent import futures
from multiprocessing import cpu_count
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform


class Dendro:
    """
    Used for the creation of a basic dendrogram that will be saved to the same directory as output specifies.
    """

    def __init__(self, labels, matrix, output):
        self.labels = labels
        self.matrix = matrix
        self.output = output
        self.linkage_matrix = ""

        self.make_hc_dendrogram()

    def make_hc_dendrogram(self):
        """
        Make hierarchical clustering tree
        :return: None. Results in the creation of a pdf file with the Dendrogram.
        """

        dists = squareform(self.matrix)  # Convert square matrix to condensed matrix
        self.linkage_matrix = linkage(dists, "ward")
        dendrogram(self.linkage_matrix, labels=self.labels, orientation='left', leaf_font_size=9)
        plt.title("Dendrogram")
        plt.tight_layout()
        output_name = self.output + "/Dendrogram.pdf"
        plt.savefig(output_name)


class IdentificationMethods:

    @staticmethod
    def collect_paired(path):
        """
        :param path: type Str: Directory path to all paired Illumina reads.
        :return: type Dict: A dictionary with the keys being the sample name, and the value being the paths for the
        forward and reverse reads.

        :requires: Requires that all samples have their respective read pair present in the directory path provided.
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
        """
        Concatenates the two files in dict

        :param paired_dict: type Dict: A dictionary with sample names as keys and file paths as values.
        :return: type Dict: A dictionary with sample names as keys, and the path to a file where both the forward and
        reverse reads are present.

        :requires: Requires that all samples have two fastq.gz files containing the forward and reverse reads,
        respectively.
        """
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
        :param path: type STR: a full file path for the directory containing fastq.gz files
        :return: type dict: Dictionary with keys being sample names, and values being the file path
        for the fastq/fastq.gz from the sequencing run

        :requires: Requires that all reads are concatenated from the sequencing run for each sample.
        """
        dirs = [sd for sd in os.listdir(path) if os.path.isdir("{}/{}".format(path, sd))]
        files_dict = dict()

        # collects all the fastq/fastq.gz files from the directories. Used for barcoded runs. Assumed that all reads
        # are concatenated
        if len(dirs) > 1:
            for directories in dirs:
                for file in os.listdir("{}/{}".format(path, directories)):
                    if os.path.isfile("{}/{}/{}".format(path, directories, file)):
                        file_name = file.split(".")[0]
                        files_dict[file_name] = "{}/{}/{}".format(path, directories, file)

        # collects all the fastq/fastq.gz for a single sample run
        # assumes that all reads are concatenated
        else:
            for f in os.listdir(path):
                if os.path.isfile("{}/{}".format(path, f)):
                    file_name = f.split(".")[0]
                    files_dict[file_name] = "{}/{}".format(path, f)

        return files_dict

    @staticmethod
    def sample_mash_screen(reference, sample_path, tsv_dict):
        """
        The individual mash screen call for all samples.

        :param reference: type Str: Reference .msh file that is pre-generated.
        :param sample_path: type Str: File path for the fastq/fastq.gz reads.
        :param tsv_dict: type Dict: Dict with sample name as keys, and file path for mash screen output as values.
        :return: None. Updates the tsv_dict input with new file locations for the mash screen output.
        """

        sample_name = os.path.basename(sample_path).split('.')[0]  # finds sample name

        output_directory = "{}/mash_screen/".format('/'.join(sample_path.split('/')[:-1]))  # path to the output dir
        output_table = "{}/mash_screen_{}_identification.tsv".format(output_directory, sample_name)  # output table

        cmd = ["mash", "screen", reference, sample_path]
        screen_mash = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = screen_mash.communicate()  # mash screen of sample against reference

        with open(output_table, "w") as ftw:  # saves the standard output as a tsv document
            ftw.write(stdout.decode('ascii'))

        if os.path.isfile(output_table):
            tsv_dict[sample_name] = output_table  # saves path for each output tsv with respect to the sample
        else:
            print("ERROR\nSomething is not alright with the output! check it out girl")

    @staticmethod
    def mash_screen(reference, dic):
        """
        Parallel call of mash screen against all samples.

        :param reference: type Str: File path to the pre-generated .msh reference
        :param dic: type Dict: Keys are sample names, and the value is the file path for basecalled reads
        :return: type Dict: Keys are the sample name, and the value is the file path for the mash screen output tsv
        """
        screen_tsv_dict = dict()

        # parallelizes the mash screen command
        with futures.ThreadPoolExecutor(max_workers=cpu_count()) as executor:
            paths = (tuple(dic.values()))
            for results in executor.map(
                    lambda x: IdentificationMethods.sample_mash_screen(reference, x, screen_tsv_dict),
                    paths):
                pass

        return screen_tsv_dict

    @staticmethod
    def mash_screen_parse(mash_dict):
        """

        :param mash_dict: type Dict: Keys are the sample name, and the value is the file path for the mash
        screen output tsv
        :return: type Dict: Keys are the sample name, and the value is the file path for the top five mash screen hits,
        as determined by shared bin values.

        """
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
                df['query_match'].tolist()))  # format the sample names

            df.to_csv(output_tsv, sep="\t")

            if os.path.isfile(output_tsv):
                os.remove(screen_path)  # removes the mash screen output with all hits
                tsv_dict[sample_name] = output_tsv  # reassigns the output path to be the top five hits

            else:
                print("ERROR:\n There was an issue collecting the top 5 hits.")

        return tsv_dict

    @staticmethod
    def mash_visualize_rawreads(sample_dict, tsv_dict, name):
        """
        :param sample_dict: type Dict:
        :param tsv_dict: type Dict:
        :param name: type Str:
        :return: None. Results in a pdf file being saved with a graph.

        :requires: Requires that only one sample for all reads.
        """
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
                            bp_count += len(line.decode('ascii').split('\n')[0])

                            counter += 1
                        elif counter == 2:
                            counter += 1
                        elif counter == 3:
                            counter = 0

            except UnicodeDecodeError:

                with gzip.open(sample_dict[sample_name], "r") as ftr:
                    print(sample_dict[sample_name])
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

            coverage = round(bp_count / 1000000)
            tsv_df = pandas.read_csv(tsv_path, sep='\t', usecols=['hash_match', 'query_match'])  # get query information
            hash_match = tsv_df['hash_match'].tolist()
            query_match = tsv_df['query_match'].tolist()
            queries = dict()
            index = 0
            while index < 5:  # based on how many samples we want
                queries[os.path.basename(query_match[index].split('.')[0])] = hash_match[index]
                index += 1

            output_table[coverage] = queries

        df = pandas.DataFrame(output_table)
        col_names = ['query_names'] + df.columns.tolist()
        df = df.reset_index()
        df.columns = col_names
        df = pandas.melt(df, id_vars=['query_names'], value_vars=list(output_table.keys()))
        df.columns = ['query names', 'million bp', 'percentage of shared bins']
        df['million bp'] = df['million bp'].astype(float)
        df['percentage of shared bins'] = df['percentage of shared bins'].fillna(0)
        df['percentage of shared bins'] = list(
            map(
                lambda x: (int(x.split('/')[0]) / (int(x.split('/')[-1]))) * 100 if isinstance(x, str) else x,
                df['percentage of shared bins'].tolist()
            )
        )
        df['percentage of shared bins'] = df['percentage of shared bins'].astype(float)

        line_plt = sns.lineplot(x='million bp', y='percentage of shared bins', hue='query names', data=df,
                                legend='full')
        box = line_plt.get_position()
        line_plt.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        hnd, lab = line_plt.get_legend_handles_labels()
        lgd = plt.legend(handles=hnd, labels=lab, bbox_to_anchor=(1.05, 1), loc=0, borderaxespad=0, fontsize='x-small',
                         title='Species')
        line_plt.get_figure().savefig(output_graph_name, bbox_extra_artists=(lgd,), bbox_inches='tight')

    @staticmethod
    def clear_pair(dic):
        """
        Removes concatenated files for paired Illumina reads

        :param dic: type Dict: Keys are sample names, and values are the file path to the concatenated Illumina reads
        :return: None. Deletes all files in dic.values()
        """
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

        with open(output, "w") as ftw:
            ftw.write(stdout.decode('ascii'))

        return output

    @staticmethod
    def parse_distancematrix(path):
        """
        Parse through the tsv by treating it as a numpy array. Converts the output lower triangle matrix to a square
        matrix.

        :param path: type Str: File path to the tsv output of mash triangle.
        :return: None. Initiates a Dendro object.
        """
        output = "/".join(path.split('/')[:-1])

        with open(path, "r") as ftr:
            lines = ftr.readlines()
            names = list(map(lambda x: [x.split('\t')[0]], lines))
            names = names[1:]
            values = list(map(lambda x: x.split('\t')[1:], lines))
            values = values[2:]
            values = list(map(lambda x: x[:-1] + x[-1].split(), values))
            values = values
        fix_values = [[0]]
        index = 0
        position = 0
        while index <= len(values):
            if index == 0:
                line = list(map(lambda x: x[position], values[index:]))

                fix_values[0] = fix_values[0] + line
                position += 1
                index += 1
            elif index == len(values):
                fix_values.append(values[-1] + [0])
                index += 1
            else:
                line = list(map(lambda x: x[position], values[index:]))
                fix_values.append(values[index - 1] + [0] + line)
                index += 1
                position += 1

        i = 0
        while i < len(fix_values):
            fix_values[i] = list(map(lambda x: float(x), fix_values[i]))
            i += 1

        labels = []
        for name in names:
            labels.append(name[0].split('/')[-1])

        array = np.array(fix_values)
        Dendro(labels, array, output)


class Identification(object):

    def run(self):
        if self.type.lower() in ["assembly", "raw-reads"]:
            if not os.path.isdir("{}/mash_screen/".format(self.path)):
                os.mkdir("{}/mash_screen/".format(self.path))
            else:
                if not self.pair:
                    start = IdentificationMethods.collect_unpaired(self.path)
                else:
                    start = IdentificationMethods.collect_paired(self.path)
        else:
            print("You must have a valid input for the type of sequence information being supplied.")

        screen_samples = IdentificationMethods.mash_screen(self.ref, start)
        final_tsv = IdentificationMethods.mash_screen_parse(screen_samples)

        if self.diagnostic:
            IdentificationMethods.mash_visualize_rawreads(start, final_tsv, self.sample_name)

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

    def __init__(self, args):

        # Command line arguments
        self.path = args.path
        self.ref = args.ref
        self.pair = args.pair
        self.tree = args.tree
        self.clear = args.clear
        self.sample_name = args.sample_name
        self.diagnostic = args.diag

        # Run
        self.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Allows for the species identification of sample raw reads via mash.\nA)
    parser.add_argument("-d",
                        help="Input directory with all reads or directories containing reads. Full path must be provided. Assumes long-reads if --paired is not specified.", dest="path",
                        required=True, metavar="/input/folder/")
    parser.add_argument("-ref",
                        help="Custom reference sketch file. Default is for Mycobacterium", required=True,
                        metavar="/path/tp/file.msh", dest="ref")  # set default
    parser.add_argument("-name",
                        help="Sample name. Used only in conjuction with the diagnostic_run, which assumes a single sample."
                        metavar="PATH", dest="sample_name", required=False)
    parser.add_argument("--paired", action='store_true',
                        default=False, dest="pair",
                        help="Provided directory is for paired end Illumina reads.", required=False)
    parser.add_argument("--diagnostic_run", action='store_true',
                        default=False, dest="diag",
                        help="Used to create a graphical output based on all fastq in directory. Assumes one sample for all different reference reads.", required=False)
    parser.add_argument("--tree", action='store_true',
                        default=False, dest="tree", help="Outputs a dendrogram, with all sequences in input and references.")
    parser.add_argument("--clear", action='store_true',
                        default=False, dest="clear",
                        help="Will remove the concatenated Illumina file, if Illumina reads were used", required=False)
    args = parser.parse_args()
    Identification(args)
