"""
soapdenovo2_pregraph_sparse.py
A wrapper script for SOAPdenovo2 pregraph sparse module
Copyright   Peter Li - GigaScience and BGI-HK
"""

import optparse
import os
import shutil
import subprocess
import sys
import tempfile
import re
import fnmatch


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def html_report_from_directory(html_out, dir):
    html_out.write('<html>\n<head>\n</head>\n<body>\n<font face="arial">\n<p>Pregraph sparse outputs</p>\n<p/>\n')
    for dirname, dirnames, filenames in os.walk(dir):
        #Link supplementary documents in HTML file
        for file in filenames:
            if fnmatch.fnmatch(file, '*pair_*'):
                continue
            else:
                html_out.write('<p><a href="%s">%s</a></p>\n' % (file, file))
    html_out.write('</font>\n</body>\n</html>\n')


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def main():
    ncpu = 4

    #Parse command line
    parser = optparse.OptionParser()
    parser.add_option('', '--file_source', dest='file_source')
    parser.add_option("", "--config", dest="config")

    parser.add_option("", "--max_read_length", dest="max_read_length")
    #Make list of params
    parser.add_option("", "--avg_ins", action="append", type="string", dest="avg_insert_list")
    parser.add_option("", "--reverse_seq", action="append", type="string", dest="reverse_seq_list")
    parser.add_option("", "--asm_flags", action="append", type="string", dest="asm_flags_list")
    parser.add_option("", "--rd_len_cutoff", action="append", type="string", dest="rd_len_cutoff_list")
    parser.add_option("", "--rank", action="append", type="string", dest="rank_list")
    parser.add_option("", "--pair_num_cutoff", action="append", type="string", dest="pair_num_cutoff_list")
    parser.add_option("", "--map_len", action="append", type="string", dest="map_len_list")

    #Data inputs
    parser.add_option("", "--type_of_data", action="append", type="string", dest="type_of_data_list")
    parser.add_option("", "--format_of_data", action="append", type="string", dest="format_of_data_list")
    parser.add_option("", "--single_fastq_input1", action="append", type="string", dest="single_fastq_input1_list")
    parser.add_option("", "--single_fastq_gzipped_input1", action="append", type="string", dest="single_fastq_gzipped_input1_list")
    parser.add_option("", "--single_fasta_input1", action="append", type="string", dest="single_fasta_input1_list")
    parser.add_option("", "--single_fasta_gzipped_input1", action="append", type="string", dest="single_fasta_gzipped_input1_list")
    parser.add_option("", "--single_bam_input1", action="append", type="string", dest="single_bam_input1_list")

    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")

    parser.add_option("", "--paired_fastq_gzipped_input1", action="append", type="string", dest="paired_fastq_gzipped_input1_list")
    parser.add_option("", "--paired_fastq_gzipped_input2", action="append", type="string", dest="paired_fastq_gzipped_input2_list")
    parser.add_option("", "--paired_fasta_gzipped_input1", action="append", type="string", dest="paired_fasta_gzipped_input1_list")
    parser.add_option("", "--paired_fasta_gzipped_input2", action="append", type="string", dest="paired_fasta_gzipped_input2_list")

    parser.add_option("", "--paired_bam_input1", action="append", type="string", dest="paired_bam_input1_list")
    parser.add_option("", "--paired_bam_input2", action="append", type="string", dest="paired_bam_input2_list")

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")

    #Mandatory params
    parser.add_option("-K", "--kmer_size", dest="kmer_size")
    parser.add_option("-z", "--genome_size", dest="genome_size")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff")
    #Commented out to keep under local control
    #parser.add_option("-p", "--ncpu", dest="ncpu")

    #Optional params
    parser.add_option("-g", "--max_kmer_edge_length", dest="max_kmer_edge_length")
    parser.add_option("-e", "--kmer_edge_freq_cutoff", dest="kmer_edge_freq_cutoff")
    parser.add_option("-R", "--output_extra_info", dest="output_extra_info")
    parser.add_option("-r", "--runmode", dest="runmode")

    #HTML output
    parser.add_option("", "--html_file", dest="html_file")
    parser.add_option("", "--html_file_files_path", dest="html_file_files_path")

    #Outputs
    parser.add_option("", "--kmer_freq", dest='kmer_freq')
    parser.add_option("", "--edge", dest='edge')
    parser.add_option("", "--mark_on_edge", dest='mark_on_edge')
    parser.add_option("", "--path", dest='path')
    parser.add_option("", "--pre_arc", dest='pre_arc')
    parser.add_option("", "--vertex", dest='vertex')
    parser.add_option("", "--pregraph_basic", dest='pregraph_basic')
    parser.add_option("", "--soap_config", dest='soap_config')
    opts, args = parser.parse_args()

    #Create directory to process and store Corrector outputs
    html_file = opts.html_file
    job_work_dir = opts.html_file_files_path

    #Need a temporary directory to perform processing
    dirpath = tempfile.mkdtemp(prefix="tmp-pregraph-sparse-")

    if opts.file_source == "history":
        config_file = opts.config
    else:
        #Create temp file to store soapdenovo2 running configuration
        config_file = tempfile.NamedTemporaryFile(dir=dirpath, prefix="soap_", suffix=".config").name

        try:
            fout = open(config_file, 'w')
            fout.write("max_rd_len=%s\n" % opts.max_read_length)
            #Calculate how many sets of data there are - use avg_ins as a measure of this
            #Separate indices required to keep count of reads
            single_read_index = 0
            paired_read_index = 0
            for index in range(len(opts.avg_insert_list)):
                fout.write("[LIB]\n")
                fout.write("avg_ins=%s\n" % opts.avg_insert_list[index])
                fout.write("reverse_seq=%s\n" % opts.reverse_seq_list[index])
                fout.write("asm_flags=%s\n" % opts.asm_flags_list[index])
                fout.write("rd_len_cutoff=%s\n" % opts.rd_len_cutoff_list[index])
                fout.write("rank=%s\n" % opts.rank_list[index])
                fout.write("pair_num_cutoff=%s\n" % opts.pair_num_cutoff_list[index])
                fout.write("map_len=%s\n" % opts.map_len_list[index])
                #Add data file configuration - needs careful looping due to single and paired reads
                print opts.type_of_data_list[index]
                print opts.format_of_data_list[index]
                if opts.type_of_data_list[index] == "single":  #then only one read
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q=%s\n" % opts.single_fastq_input1_list[single_read_index])
                    elif opts.format_of_data_list[index] == "fastq_gzipped":
                        #Copy file into temp directory and give it a gz suffix
                        print "File: ", opts.single_fastq_gzipped_input1_list[single_read_index]
                        shutil.copy2(dirpath, opts.single_fastq_gzipped_input1_list[single_read_index] + '.gz')
                        fout.write("f=" + dirpath + "%s.gz\n" % opts.single_fastq_gzipped_input1_list[single_read_index])
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f=%s\n" % opts.single_fasta_input1_list[single_read_index])
                    elif opts.format_of_data_list[index] == "fasta_gzipped":
                        print "Ok here!"
                        #Copy file into temp directory and give it a gz suffix
                        print "File: ", opts.single_fasta_gzipped_input1_list[single_read_index]
                        shutil.copy2(opts.single_fasta_gzipped_input1_list[single_read_index], opts.single_fasta_gzipped_input1_list[single_read_index] + '.gz')
                        fout.write("f=" + "%s.fa.gz\n" % opts.single_fasta_gzipped_input1_list[single_read_index])
                    else:
                        fout.write("b=%s\n" % opts.single_bam_input1_list[single_read_index])
                    single_read_index = + 1
                elif opts.type_of_data_list[index] == "paired":
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q1=%s\n" % opts.paired_fastq_input1_list[paired_read_index])
                        fout.write("q2=%s\n" % opts.paired_fastq_input2_list[paired_read_index])
                    elif opts.format_of_data_list[index] == "fastq_gzipped":
                        print "Ok here!"
                        #Copy file into temp directory and give it a gz suffix
                        print "File: ", opts.paired_fastq_gzipped_input1_list[paired_read_index]
                        shutil.copy2(opts.paired_fastq_gzipped_input1_list[paired_read_index], opts.paired_fastq_gzipped_input1_list[paired_read_index] + '.fq.gz')
                        shutil.copy2(opts.paired_fastq_gzipped_input2_list[paired_read_index], opts.paired_fastq_gzipped_input2_list[paired_read_index] + '.fq.gz')
                        fout.write("q1=" + "%s.fq.gz\n" % opts.paired_fastq_gzipped_input1_list[paired_read_index])
                        fout.write("q2=" + "%s.fq.gz\n" % opts.paired_fastq_gzipped_input2_list[paired_read_index])
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("f2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    elif opts.format_of_data_list[index] == "fasta_gzipped":
                        print "Ok here!"
                        #Copy file into temp directory and give it a gz suffix
                        print "File: ", opts.paired_fasta_gzipped_input1_list[paired_read_index]
                        shutil.copy2(opts.paired_fasta_gzipped_input1_list[paired_read_index], opts.paired_fasta_gzipped_input1_list[paired_read_index] + '.fa.gz')
                        shutil.copy2(opts.paired_fasta_gzipped_input2_list[paired_read_index], opts.paired_fasta_gzipped_input2_list[paired_read_index] + '.fa.gz')
                        fout.write("f1=" + "%s.fa.gz\n" % opts.paired_fasta_gzipped_input1_list[paired_read_index])
                        fout.write("f2=" + "%s.fa.gz\n" % opts.paired_fasta_gzipped_input2_list[paired_read_index])
                    else:
                        fout.write("b1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("b2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    paired_read_index = + 1
            fout.close()
        except Exception, e:
            stop_err("File cannot be opened for writing soap.config: " + str(e))

    #Create correct paths to html-linked files
    rex = re.compile('database/(.*)/dataset_')
    #Create replacement text using html_file
    #Split string into tokens
    tokens = html_file.split("/")
    #Get second to last token
    userId = tokens[len(tokens) - 2]
    files_dir = rex.sub("database/files/" + userId + "/dataset_", job_work_dir)
    files_dir = files_dir + "/"
    # print "New html dir: ", files_dir
    #Create directory
    if not os.path.exists(files_dir):
        try:
            os.makedirs(files_dir)
        except:
            pass

    #Set up command line call
    if int(opts.kmer_size) <= 63 and opts.default_full_settings_type == "default":
        cmd = "Pregraph_Sparse_63mer.v1.0.3 -s %s -K %s -z %s -o %s -d %s" % (config_file, opts.kmer_size, opts.genome_size, files_dir + "/out", opts.kmer_freq_cutoff)
    elif int(opts.kmer_size) <= 63 and opts.default_full_settings_type == "full":
        cmd = "Pregraph_Sparse_63mer.v1.0.3 -s %s -K %s -z %s -o %s -d %s -g %s -e %s -R %s -r %s -p %s" % (config_file, opts.kmer_size, opts.genome_size, files_dir + "/out", opts.kmer_freq_cutoff, opts.max_kmer_edge_length, opts.kmer_edge_freq_cutoff, opts.output_extra_info, opts.runmode, ncpu)
    elif int(opts.kmer_size) > 63 and opts.default_full_settings_type == "default":
        cmd = "Pregraph_Sparse_127mer.v1.0.3 -s %s -K %s -z %s -o %s -d %s" % (config_file, opts.kmer_size, opts.genome_size, files_dir + "/out", opts.kmer_freq_cutoff)
    elif int(opts.kmer_size) > 63 and opts.default_full_settings_type == "full":
        cmd = "Pregraph_Sparse_127mer.v1.0.3 -s %s -K %s -z %s -o %s -g %s -d %s -e %s -R %s -r %s -p %s" % (config_file, opts.kmer_size, opts.genome_size, files_dir + "/out", opts.max_kmer_edge_length, opts.kmer_freq_cutoff, opts.kmer_edge_freq_cutoff, opts.output_extra_info, opts.runmode, ncpu)

    # print cmd

    #Perform SOAPdenovo2_pregraph sparse analysis
    buffsize = 1048576
    try:

        tmp_out_file = tempfile.NamedTemporaryFile(dir=files_dir).name
        tmp_stdout = open(tmp_out_file, 'w')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=files_dir).name #Contains pregraph's stdout
        tmp_stderr = open(tmp_err_file, 'w')

        #Call SOAPdenovo2
        #New additional datasets must be placed in the directory provided by $new_file_path__
        proc = subprocess.Popen(args=cmd, shell=True, cwd=files_dir, stdout=tmp_stdout, stderr=tmp_stderr.fileno())
        returncode = proc.wait()

        #Read tool stdout into galaxy stdout
        f = open(tmp_err_file)
        lines = f.readlines()
        for line in lines:
            sys.stdout.write(line)
        f.close()

        #Get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err_file, 'r')
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass

        #Close streams
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        raise Exception, 'Problem performing pregraph sparse process ' + str(e)

    #Read files into their outputs
    kmer_freq_out = open(opts.kmer_freq, 'wb')
    f = open(files_dir + "/out.kmerFreq")
    for line in f:
        kmer_freq_out.write(line)
    kmer_freq_out.close()
    f.close()

    edge_gz_out = open(opts.edge, 'wb')
    with open(files_dir + "/out.edge.gz", mode='rb') as f: # b is important -> binary
        fileContent = f.read()
        edge_gz_out.write(fileContent)
    edge_gz_out.close()
    f.close()

    pre_arc_out = open(opts.pre_arc, 'wb')
    f = open(files_dir + "/out.preArc")
    for line in f:
        pre_arc_out.write(line)
    pre_arc_out.close()
    f.close()

    vertex_out = open(opts.vertex, 'wb')
    f = open(files_dir + "/out.vertex")
    for line in f:
        vertex_out.write(line)
    vertex_out.close()
    f.close()

    pregraph_basic_out = open(opts.pregraph_basic, 'wb')
    f = open(files_dir + "/out.preGraphBasic")
    for line in f:
        pregraph_basic_out.write(line)
    pregraph_basic_out.close()
    f.close()

    config_out = open(opts.soap_config, 'w')
    f = open(config_file)
    for line in f:
        config_out.write(line)
    config_out.close()
    f.close()

    #Delete files not being linked on web page
    os.remove(files_dir + "/out.kmerFreq")
    os.remove(files_dir + "/out.edge.gz")
    os.remove(files_dir + "/out.preArc")
    os.remove(files_dir + "/out.vertex")
    os.remove(files_dir + "/out.preGraphBasic")
    files = os.listdir(files_dir)
    for f in files:
        if f.startswith("tmp"):
            os.remove(os.path.join(files_dir, f))

    #Generate html
    html_report_from_directory(open(html_file, 'w'), files_dir)

    #Clean up temp files
    cleanup_before_exit(dirpath)
    #Check results in output file
    if os.path.getsize(opts.preGraphBasic) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__":
    main()
