import os, shutil, subprocess, gzip, logging

reference = ""  # hard-coded file path for reference files


# built with the selective primer finder

# Can we convert this reference file into a dictionary with key words and strings?


class IdentifyBactMethods(object):

    @staticmethod
    def ref_to_dict(file_path):
        """
        Converts the fastq file to a dictionary with the sample key and each reference sequence

        returns a dictionary
        """
        pass

    @staticmethod
    def bbduk_sample(file_path, reference):
        """
        Takes a list of full file paths to raw unpaired reads (minION or short) and
        calls bbduk on a reference file and the sequences.
        """

        # The structure of references will be a fasta file OR a dictionary, with key words indicating genus and species
        # Each sample runs once, but the bbduk command

        stringpath = '/'.join(file_path.split('/')[:-1])
        detected_matches = "{}/detected_matches".format(stringpath)
        sample_file_name = os.path.basename(file_path).split(".")[0]
        os.mkdir(
            detected_matches)  # maybe make this before we call the function, so all sampled get stuck together <- I think this is the move

        for k, h in reference.items():
            """
            Run bbduk on each reference file for the input sample, which will be a file path as a string 

            Results in: the output of a log file, and summary .csv files with information per each hit 
            """

            logging.log("We are searching for matches between {} in {} organisms\n".format(sample_file_name, k))
            reference[k] = ",".join(h)  # should this be done before??
            output_file_name = "{}/{}_{}_matches.fastq.gz".format(detected_matches, sample_file_name,
                                                                  k)  # or a fastq.gz !!
            logfile = "{}/{}_{}.log".format(stringpath, k, sample_file_name)
            outfile = open(logfile, "w")

            cmd = ["bbduk.sh", "in={}".format(file_path),
                   "outm={}".format(output_file_name),
                   "literal={}".format(reference[k])]  # make it in /usr/local/bin or something
            p = subprocess.Popen(cmd, stdout=outfile)
            p.communicate()
            outfile.close()

            if os.stat(output_file_name).st_size == 0:
                os.remove(output_file_name)

        detected_match_contents = os.listdir(detected_matches)
        size_of_matches = [f for f in detected_match_contents if
                           os.stat("{}/{}".format(detected_matches, f)).st_size != 0]
        if len(size_of_matches) != 0:
            for f in detected_match_contents:
                # call parser ?
                pass
        else:
            os.rmdir(detected_matches)
            logging.log("No matches were detected from sample {} and our database.".format(os.path.basename(file_path)))

    @staticmethod
    def bbduk_output_parser(file_path):
        # what to look for!!
        file_name = file_path.split('/')[-1].split('.')[0]
        file_dir = '/'.join(file_path.split('/')[:-1])
        csv_summary = "{}/{}_summary.csv".format(file_dir, file_name)
        csv_dict = dict()

        with gzip.open(file_path, "r") as ftr:
            # collect information about fastq file
            # read
            # sample ID
            # length of reads <- produce an average of this? Does this even matter
            pass

        # covert dictionary into our dataframe, saving as a csv


class IdentifyBact(object):

    def __main__(self):
        # Run relevant commands
        pass

    def __init__(self):
        # directory or file paths
        pass

    def bbduk_all(self):
        # calls appropriate static methods for read pulling
        pass
