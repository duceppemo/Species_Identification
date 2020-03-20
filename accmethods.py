"""
The structure of the raw reads will be all in one fastq.gz file, for minION

So. We would apply this to a bunch of them! Straight after basecalling!!
"""

import os


class AccMethods(object):
    @staticmethod
    def find_fastq(path):
        """
        ACCOMODATE: fasta assemblies
        results in a list of all files with "fastq" or "fq", and creates a list.

        :param path: type STR: absolute file path containing all fastq files.
        :return: sorted list of all fastq files in the provided path
        """
        fastq_files = list(filter(lambda x: "fastq" in x or "fq" in x, os.listdir(path)))
        fastq_files = list(map(lambda x: path + x, fastq_files))
        fastq_files = sorted(fastq_files)

        assert fastq_files, 'No FASTQ files detected in {pth}'.format(pth=path)

        return fastq_files

    @staticmethod
    def filer_im(lis):
        """

        :param lis:
        :return:
        """
        sample_path = dict()
        sample_name = set(filter(lambda x: x.split('/')[-1].split("_")[0], lis))

        for sample in sample_name:
            for path in lis:
                if sample in path:
                    try:
                        sample_path[sample].add(path)
                    except KeyError:
                        sample_path[sample] = [path]

        return sample_path










