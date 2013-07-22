"""
soapdenovo2_pregraph.py
A wrapper script for SOAPdenovo2 pregraph module
Copyright   Peter Li - GigaScience and BGI-HK
"""

import optparse
import os
import shutil
import subprocess
import sys
import tempfile


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def main():
    #Hardcoded ncpu to keep control of thread number
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
    #Commented out to keep control of thread number
    #parser.add_option("-p", "--ncpu", dest="ncpu")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff")

    #Custom params
    parser.add_option("-a", "--init_memory_assumption", dest="init_memory_assumption")
    parser.add_option("-R", "--output_extra_info", dest="output_extra_info")

    #Outputs
    parser.add_option("", "--pregraph_basic", dest='pregraph_basic')
    parser.add_option("", "--vertex", dest='vertex')
    parser.add_option("", "--pre_arc", dest='pre_arc')
    parser.add_option("", "--edge", dest='edge')
    parser.add_option("", "--kmer_freq", dest='kmer_freq')
    parser.add_option("", "--soap_config", dest='soap_config')
    opts, args = parser.parse_args()

    #Need a temporary directory to perform processing
    dirpath = tempfile.mkdtemp(prefix="tmp-pregraph-")

    if opts.file_source == "history":
        config_file = opts.config
    else:
        #Create temp file to store soapdenovo2 running configuration
        config_file = tempfile.NamedTemporaryFile(dir=dirpath, prefix="soap_", suffix=".config").name

        try:
            fout = open(config_file,'w')
            fout.write("max_rd_len=%s\n" % opts.max_read_length)
            #Calculate how many sets of data there are - use avg_ins as a measure of this
            #Variables to keep count of different single file types
            single_read_index = 0
            fastq_single_read_index = 0
            fastq_gz_single_read_index = 0
            fasta_single_read_index = 0
            fasta_gz_single_read_index = 0
            bam_single_read_index = 0

            #Variables to keep count of different paired file types
            paired_read_index = 0
            fastq_paired_read_index = 0
            fastq_gz_paired_read_index = 0
            fasta_paired_read_index = 0
            fasta_gz_paired_read_index = 0
            bam_paired_read_index = 0
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
                if opts.type_of_data_list[index] == "single":  #then only one read
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q=%s\n" % opts.single_fastq_input1_list[fastq_single_read_index])
                        fastq_single_read_index += 1
                    elif opts.format_of_data_list[index] == "fastq_gzipped":
                        #Copy file into temp directory and give it a gz suffix
                        print "File: ", opts.single_fastq_gzipped_input1_list[fastq_gz_single_read_index]
                        shutil.copy2(dirpath, opts.single_fastq_gzipped_input1_list[fastq_gz_single_read_index] + '.gz')
                        fout.write("f=" + dirpath + "%s.gz\n" % opts.single_fastq_gzipped_input1_list[fastq_gz_single_read_index])
                        fastq_gz_single_read_index += 1
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f=%s\n" % opts.single_fasta_input1_list[fasta_single_read_index])
                        fasta_single_read_index += 1
                    elif opts.format_of_data_list[index] == "fasta_gzipped":
                        #Copy file into temp directory and give it a gz suffix
                        shutil.copy2(opts.single_fasta_gzipped_input1_list[fasta_gz_single_read_index], opts.single_fasta_gzipped_input1_list[fasta_gz_single_read_index] + '.gz')
                        fout.write("f=" + "%s.fa.gz\n" % opts.single_fasta_gzipped_input1_list[fasta_gz_single_read_index])
                        fasta_gz_single_read_index += 1
                    else:
                        fout.write("b=%s\n" % opts.single_bam_input1_list[bam_single_read_index])
                        bam_single_read_index = + 1
                    single_read_index += 1
                elif opts.type_of_data_list[index] == "paired":
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q1=%s\n" % opts.paired_fastq_input1_list[fastq_paired_read_index])
                        fout.write("q2=%s\n" % opts.paired_fastq_input2_list[fastq_paired_read_index])
                        fastq_paired_read_index += 1
                    elif opts.format_of_data_list[index] == "fastq_gzipped":
                        #Copy file into temp directory and give it a gz suffix
                        shutil.copy2(opts.paired_fastq_gzipped_input1_list[fastq_gz_paired_read_index], opts.paired_fastq_gzipped_input1_list[fastq_gz_paired_read_index] + '.fq.gz')
                        shutil.copy2(opts.paired_fastq_gzipped_input2_list[fastq_gz_paired_read_index], opts.paired_fastq_gzipped_input2_list[fastq_gz_paired_read_index] + '.fq.gz')
                        fout.write("q1=" + "%s.fq.gz\n" % opts.paired_fastq_gzipped_input1_list[fastq_gz_paired_read_index])
                        fout.write("q2=" + "%s.fq.gz\n" % opts.paired_fastq_gzipped_input2_list[fastq_gz_paired_read_index])
                        fastq_gz_paired_read_index += 1
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f1=%s\n" % opts.paired_fasta_input1_list[fasta_paired_read_index])
                        fout.write("f2=%s\n" % opts.paired_fasta_input2_list[fasta_paired_read_index])
                        fasta_paired_read_index += 1
                    elif opts.format_of_data_list[index] == "fasta_gzipped":
                        #Copy file into temp directory and give it a gz suffix
                        shutil.copy2(opts.paired_fasta_gzipped_input1_list[fasta_gz_paired_read_index], opts.paired_fasta_gzipped_input1_list[fasta_gz_paired_read_index] + '.fa.gz')
                        shutil.copy2(opts.paired_fasta_gzipped_input2_list[fasta_gz_paired_read_index], opts.paired_fasta_gzipped_input2_list[fasta_gz_paired_read_index] + '.fa.gz')
                        fout.write("f1=" + "%s.fa.gz\n" % opts.paired_fasta_gzipped_input1_list[fasta_gz_paired_read_index])
                        fout.write("f2=" + "%s.fa.gz\n" % opts.paired_fasta_gzipped_input2_list[fasta_gz_paired_read_index])
                        fasta_gz_paired_read_index += 1
                    else:
                        fout.write("b1=%s\n" % opts.paired_fasta_input1_list[bam_paired_read_index])
                        fout.write("b2=%s\n" % opts.paired_fasta_input2_list[bam_paired_read_index])
                        bam_paired_read_index += 1
                    paired_read_index += 1
            fout.close()
        except Exception, e:
            stop_err("File cannot be opened for writing soap.config: " + str(e))

    #Set up command line call
    parser.add_option("-K", "--kmer_size", dest="kmer_size")
    #Commented out to keep control of thread number
    #parser.add_option("-p", "--ncpu", dest="ncpu")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff")

    #Custom params
    parser.add_option("-a", "--init_memory_assumption", dest="init_memory_assumption")
    parser.add_option("-R", "--output_extra_info", dest="output_extra_info")

    if int(opts.kmer_size) <= 63 and opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-63mer_v2.0 pregraph -s %s -o %s -K %s -p %s -d %s" % (config_file, dirpath + "/out", opts.kmer_size, ncpu, opts.kmer_freq_cutoff)
    elif opts.kmer_size <= 63 and opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo-63mer_v2.0 pregraph -s %s -o %s -K %s -p %s -d %s -a %s -R %s" % (config_file, dirpath + "/out", opts.kmer_size, ncpu, opts.kmer_freq_cutoff, opts.init_mem_assumption, opts.output_extra_info)
    elif opts.kmer_size > 63 and opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-127mer_v2.0 pregraph -s %s -o %s -K %s -p %s -d %s" % (config_file, dirpath + "/out", opts.kmer_size, ncpu, opts.kmer_freq_cutoff)
    elif opts.kmer_size > 63 and opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo-127mer_v2.0 pregraph -s %s -o %s -K %s -p %s -d %s -a %s -R %s" % (config_file, dirpath + "/out", opts.kmer_size, ncpu, opts.kmer_freq_cutoff, opts.init_mem_assumption, opts.output_extra_info)

    #print cmd

    #Perform SOAPdenovo2_pregraph analysis
    buffsize = 1048576
    try:

        tmp_out_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stdout = open(tmp_out_file, 'w')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=dirpath).name #Contains pregraph's stdout
        tmp_stderr = open(tmp_err_file, 'w')

        #Call SOAPdenovo2
        #New additional datasets must be placed in the directory provided by $__new_file_path__
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr.fileno())
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
        raise Exception, 'Problem performing pregraph process ' + str(e)

    #Read files into their outputs
    kmer_freq_out = open(opts.kmer_freq, 'wb')
    f = open(dirpath + "/out.kmerFreq")
    for line in f:
        kmer_freq_out.write(line)
    kmer_freq_out.close()
    f.close()

    edge_gz_out = open(opts.edge, 'wb')
    with open(dirpath + "/out.edge.gz", mode='rb') as f: # b is important -> binary
        fileContent = f.read()
        edge_gz_out.write(fileContent)
    edge_gz_out.close()
    f.close()

    pre_arc_out = open(opts.pre_arc, 'wb')
    f = open(dirpath + "/out.preArc")
    for line in f:
        pre_arc_out.write(line)
    pre_arc_out.close()
    f.close()

    vertex_out = open(opts.vertex, 'wb')
    f = open(dirpath + "/out.vertex")
    for line in f:
        vertex_out.write(line)
    vertex_out.close()
    f.close()

    pregraph_basic_out = open(opts.pregraph_basic, 'wb')
    f = open(dirpath + "/out.preGraphBasic")
    for line in f:
        pregraph_basic_out.write(line)
    pregraph_basic_out.close()
    f.close()

    config_out = open(opts.soap_config, 'wb')
    f = open(config_file)
    for line in f:
        config_out.write(line)
    config_out.close()
    f.close()

    #Clean up temp files
    cleanup_before_exit(dirpath)
    #Check results in output file
    if os.path.getsize(opts.pregraph_basic) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__":
    main()
