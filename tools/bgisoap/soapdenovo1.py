"""
soapdenovo1.py
A wrapper script for SOAPdenovo-1.0
Copyright   Peter Li - GigaScience and BGI-HK
            Huayan Gao - CUHK
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def main():
    #Parse command line
    parser = optparse.OptionParser()
    parser.add_option("", "--max_read_length", dest="max_read_length")
    parser.add_option("", "--file_source", dest="file_source")
    parser.add_option("", "--configuration", dest="configuration")

    parser.add_option("", "--avg_ins", action="append", type="string", dest="avg_insert_list")
    parser.add_option("", "--reverse_seq", action="append", type="string", dest="reverse_seq_list")
    parser.add_option("", "--asm_flags", action="append", type="string", dest="asm_flags_list")
    parser.add_option("", "--rank", action="append", type="string", dest="rank_list")
    #Data inputs
    parser.add_option("", "--type_of_data", action="append", type="string", dest="type_of_data_list")
    parser.add_option("", "--format_of_data", action="append", type="string", dest="format_of_data_list")
    parser.add_option("", "--single_fastq_input1", action="append", type="string", dest="single_fastq_input1_list")
    parser.add_option("", "--single_fasta_input1", action="append", type="string", dest="single_fasta_input1_list")
    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("", "--kmer_size", dest="kmer_size")
    parser.add_option("", "--ncpu", dest="ncpu")
    parser.add_option("", "--delete_kmers_freq_one", dest="delete_kmers_freq_one")
    parser.add_option("", "--delete_edges_coverage_one", dest="delete_edges_coverage_one")
    parser.add_option("", "--unsolve_repeats", dest="unsolve_repeats")

    #Outputs
    parser.add_option("", "--contig", dest='contig', help="Contig sequence file")
    parser.add_option("", "--scafseq", dest='scafseq', help="Scaffold sequence file")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp(prefix="tmp-soapdenovo1-")
    #Create temp file for configuration
    
    if opts.file_source == "history":
        config_file = opts.configuration
    else:
        #Create temp file for configuration
        config_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        try:
            fout = open(config_file,'w')
            fout.write("max_rd_len=%s\n" % opts.max_read_length)
            #Calculate how sets of data there are - use avg_ins as a measure of this
            #Loop thru this number of times
            #Also use separate index to keep count of reads
            single_read_index = 0
            paired_read_index = 0
            for index in range(len(opts.avg_insert_list)):
            #                print "single_read_index ", single_read_index
            #                print "paired_read_index ", paired_read_index
                fout.write("[LIB]\n")
                fout.write("avg_ins=%s\n" % (opts.avg_insert_list)[index])
                fout.write("reverse_seq=%s\n" % opts.reverse_seq_list[index])
                fout.write("asm_flags=%s\n" % opts.asm_flags_list[index])
                fout.write("rank=%s\n" % opts.rank_list[index])
                #Add data file configuration - needs careful looping due to single and paired reads
                #                print opts.type_of_data_list[index]
                #                print opts.format_of_data_list[index]
                if opts.type_of_data_list[index] == "single":  #then only one read
                    if (opts.format_of_data_list)[index] == "fastq":
                        fout.write("q=%s\n" % (opts.single_fastq_input1_list)[single_read_index])
                    elif opts.format_of_data == "fasta":
                        fout.write("f=%s\n" % opts.single_fasta_input1_list[single_read_index])
                    else:
                        fout.write("b=%s\n" % opts.single_bam_input1_list[single_read_index])
                    single_read_index = single_read_index + 1
                elif opts.type_of_data_list[index] == "paired":
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q1=%s\n" % (opts.paired_fastq_input1_list)[paired_read_index])
                        fout.write("q2=%s\n" % (opts.paired_fastq_input2_list)[paired_read_index])
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("f2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    paired_read_index = paired_read_index + 1
            fout.close()
        except Exception, e:
            stop_err("config file cannot be opened for writing: " + str(e))

    if opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo_v1.0 all -s %s -o %s" % (config_file, tmp_dir + "/result")
    elif opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo_v1.0 all -s %s -o %s -K %s -p %s -d %s -D %s -R %s" % (config_file, tmp_dir + "/result", opts.kmer_size, opts.ncpu, opts.delete_kmers_freq_one, opts.delete_edges_coverage_one, opts.unsolve_repeats)
    #print cmd

    #Perform SOAPdenovo-1.0 analysis
    buffsize = 1048576
    try:
        #Create file in temporary directory
        tmp = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #Open a stream to file
        tmp_stderr = open(tmp, 'wb')
        #Call SOAPdenovo-trans
        proc = subprocess.Popen(args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
        #Close stream
        tmp_stderr.close()
        #Get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp, 'rb')
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        raise Exception, 'Problem performing SOAPdenovo 1 ' + str(e)

    #Read SOAPdenovo 1 results into outputs
    contig_out = open(opts.contig, 'wb')
    file = open(tmp_dir + '/result.contig')
    for line in file:
        #print line
        contig_out.write(line)
    contig_out.close()

    scafseq_out = open(opts.scafseq, 'wb')
    file = open(tmp_dir + '/result.scafSeq')
    for line in file:
        #print line
        scafseq_out.write(line)
    scafseq_out.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.contig) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__":
    main()
