import argparse, os, subprocess
from multiprocessing import cpu_count
from concurrent import futures


class ShortReads:
    def __init__(self, args):
        self.reference = args.reference
        self.path = args.path
        self.output = args.output

        if not os.path.isdir(self.output):
            os.mkdir(self.output)
        if not os.path.isdir("{}/bam_processing".format(self.output)):
            os.mkdir("{}/bam_processing".format(self.output))

    def run(self):

        subprocess.run(["bwa", "index", self.reference])
        short_reads = MethodsCalling.short_read_sorting(self.path)
        bam_files_dict = MethodsCalling.parallel_bam_call(short_reads, self.reference, self.output)
        bam_path = MethodsCalling.bam_file_output(bam_files_dict, self.output)
        # merge_bam_path = MethodsCalling.samtools_merge(bam_path, "{}/bam_processing".format(self.output))
        vcf_path = MethodsCalling.vcf_run(bam_path, self.reference, self.output)

        print("VCF written to {}.".format(vcf_path))

        # Remove all temporary files produced


class MethodsCalling:

    @staticmethod
    def short_read_sorting(path):
        pair_files = [f for f in os.listdir(path) if os.path.isfile("{}/{}".format(path, f))]
        sample_name = set(map(lambda x: x.split("_")[0], pair_files))

        short_dict = dict()

        for sample in sample_name:
            short_dict[sample] = []
            for files in pair_files:
                if sample in files:
                    short_dict[sample].append("{}/{}".format(path, files))
        return short_dict

    @staticmethod
    def bam_formation(sample_name, read1, read2, reference, output, bam_dict):
        # Takes in a dictionary, makes a bam file and indexes it. Returns a dictionary to the full bam file path

        rg = "@RG\\tID:{}\\tSM:{}".format(sample_name, sample_name)  # make a read group header for each set of reads

        bam_output = "{}/bam_processing/{}".format(output, sample_name)
        bwa_mem_cmd = ["bwa", "mem",
                       "-M",
                       "-t", "{}".format(len(os.sched_getaffinity(0))),
                       "-R", rg,
                       reference, read1, read2]

        st_view_cmd = ["samtools", "view",
                       "-h", "-b",
                       "-F", "4"]

        st_sort_cmd = ["samtools", "sort",
                       "-o", "{}_dup.bam".format(bam_output)]

        rm_dup_cmd = ["samtools", "rmdup",
                      "{}_dup.bam".format(bam_output),
                      "{}_final.bam".format(bam_output)]

        st_index = ["samtools", "index", "{}_final.bam".format(bam_output)]

        bwa = subprocess.Popen(bwa_mem_cmd,
                               stdout=subprocess.PIPE)  # pipe out bwa alignment with short reads to assembly
        st_view = subprocess.Popen(st_view_cmd, stdin=bwa.stdout,
                                   stdout=subprocess.PIPE)  # pipe out unmapped reads removed bam file
        bwa.stdout.close()
        st_sort = subprocess.Popen(st_sort_cmd, stdin=st_view.stdout, stdout=subprocess.PIPE)
        st_view.stdout.close()
        st_sort.communicate()
        subprocess.run(rm_dup_cmd)
        subprocess.run(st_index)

        bam_dict[sample_name] = "{}_final.bam".format(bam_output)

    @staticmethod
    def bam_file_output(dictionary, output):
        output_summ = "{}/bam_processing/bam_summ.txt".format(output)

        with open(output_summ, "w") as ftw:
            for path in dictionary.values():
                ftw.write("{}\n".format(path))

        return output_summ

    # @staticmethod
    # def samtools_merge(path, output):
    #     # DO NOT USE THIS!
    #
    #     # Merges all bam files, returns path to output bam
    #     output_bam = "{}/all_short_merged.bam".format(output)
    #     cmd = ["samtools", "merge",
    #            "-n", "-u",
    #            "-@", str(cpu_count()),
    #            "-b", path,
    #            output_bam]
    #
    #     subprocess.run(cmd)
    #     return output_bam

    @staticmethod
    def vcf_run(path, reference, output):
        # Results in a vcf file !

        vcf_final_output = "{}/shortread_variant_call.vcf".format(output)
        bcf_tools_cmd1 = ["bcftools", "mpileup",
                          "--threads", "{}".format(cpu_count()),
                          "-O", "u",
                          "--fasta-ref", reference,
                          "-b", path]

        bcf_tools_cmd2 = ["bcftools", "call",
                          "-Ov", "-m",
                          "-v", "--ploidy", "1",
                          "-o", vcf_final_output]

        p1 = subprocess.Popen(bcf_tools_cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(bcf_tools_cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p2.communicate()

        return vcf_final_output

    @staticmethod
    def parallel_bam_call(dictionary, reference, output):
        bam_dict = dict()

        with futures.ThreadPoolExecutor(max_workers=cpu_count()) as executor:
            bam_files = dictionary.items()  # get sample name and the list

            for results in executor.map(
                    lambda x: MethodsCalling.bam_formation(x[0],
                                                           x[1][0],
                                                           x[1][1],
                                                           reference, output,
                                                           bam_dict), bam_files):
                pass

        return bam_dict


if __name__ == "__main__":
    parse = argparse.ArgumentParser(
        description="Creates a vcf files from multiple illumina samples against the same reference.")
    parse.add_argument("-sd", help="Directory with paired end short reads", dest="path", required=True)
    parse.add_argument("-o", help="Output directory", dest="output", required=True)
    parse.add_argument("-reference", help="FASTA assembly", dest="reference", required=True)

    args = parse.parse_args()

    sr = ShortReads(args)
    sr.run()
