"""
soapdenovo2.py
A wrapper script for SOAPdenovo2
Peter Li - GigaScience and BGI-HK
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
    parser.add_option("", "--single_fasta_input1", action="append", type="string", dest="single_fasta_input1_list")
    parser.add_option("", "--single_bam_input1", action="append", type="string", dest="single_bam_input1_list")

    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")
    parser.add_option("", "--paired_bam_input1", action="append", type="string", dest="paired_bam_input1_list")
    parser.add_option("", "--paired_bam_input2", action="append", type="string", dest="paired_bam_input2_list")

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")

    #Custom params
    parser.add_option("-K", "--kmer_size", dest="kmer_size")
    parser.add_option("-p", "--ncpu", dest="ncpu")
    parser.add_option("-a", "--init_memory_assumption", dest="init_memory_assumption")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff")
    parser.add_option("-R", "--resolve_repeats", dest="resolve_repeats")
    parser.add_option("-D", "--edge_cov_cutoff", dest="edge_cov_cutoff")
    parser.add_option("-M", "--merge_level", dest="merge_level")
    parser.add_option("-m", "--max_k", dest="max_k")
    parser.add_option("-e", "--weight", dest="weight")
    #parser.add_option("-r", "--keep_avail_read", dest="keep_avail_read")
    parser.add_option("-E", "--merge_clean_bubble", dest="merge_clean_bubble")
    #parser.add_option("-f", "--output_gap_related_reads", dest="output_gap_related_reads")
    parser.add_option("-k", "--kmer_r2c", dest="kmer_r2c")
    parser.add_option("-F", "--fill_gaps", dest="fill_gaps")
    parser.add_option("-u", "--unmask_contigs", dest="unmask_contigs")
    parser.add_option("-w", "--keep_contigs_connected", dest="keep_contigs_connected")
    parser.add_option("-G", "--gap_len_diff", dest="gap_len_diff")
    parser.add_option("-L", "--min_contig_len", dest="min_contig_len")
    parser.add_option("-c", "--min_contig_cvg", dest="min_contig_cvg")
    parser.add_option("-C", "--max_contig_cvg", dest="max_contig_cvg")
    parser.add_option("-b", "--insert_size_upper_bound", dest="insert_size_upper_bound")
    parser.add_option("-B", "--bubble_coverage", dest="bubble_coverage")
    parser.add_option("-N", "--genome_size", dest="genome_size")
    parser.add_option("-V", "--ass_visual", dest="ass_visual")

    #Outputs
    parser.add_option("", "--contig", dest="contig")
    parser.add_option("", "--scafseq", dest="scafseq")
    parser.add_option("", "--config", dest="config")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp(prefix="tmp-soapdenovo2-all-")

    if opts.file_source == "history":
        config_file = opts.configuration
    else:
        #Create temp file for configuration
        config_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        try:
            fout = open(config_file,'wb')
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
                fout.write("rd_len_cutoff=%s\n" % opts.rd_len_cutoff_list[index])
                fout.write("rank=%s\n" % opts.rank_list[index])
                fout.write("pair_num_cutoff=%s\n" % opts.pair_num_cutoff_list[index])
                fout.write("map_len=%s\n" % opts.map_len_list[index])
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
                    else:
                        fout.write("b1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("b2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    paired_read_index = paired_read_index + 1
            fout.close()
        except Exception, e:
            stop_err("config file cannot be opened for writing" + str(e))

    if opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-63mer_v2.0 all -s %s -o %s" % (config_file, tmp_dir + "/result")
    elif opts.default_full_settings_type == "full":
        #Check important params
        if int(opts.max_k) <= int(opts.kmer_size):
            sys.stderr.write("Problem: The value for the maximum kmer parameter has to be an odd value and higher than the kmer size.\n\n")
        if int(opts.kmer_r2c) < int(opts.kmer_size):
            sys.stderr.write("Problem: The kmer size used for mapping reads to contigs should be at least equal to the kmer size or higher.\n\n")
        #kmer size is only set in full param analysis
        if int(opts.kmer_size) <= 63:
            cmd = "SOAPdenovo-63mer_v2.0 all -s %s -o %s -K %s -p %s -a %s -d %s -R %s -D %s -M %s -m %s -e %s -E %s -k %s -F %s -u %s -w %s -G %s -L %s -c %s -C %s -b %s -B %s -N %s -V %s" % (config_file, tmp_dir + "/result", opts.kmer_size, opts.ncpu, opts.init_memory_assumption, opts.kmer_freq_cutoff, opts.resolve_repeats, opts.edge_cov_cutoff, opts.merge_level, opts.max_k, opts.weight, opts.merge_clean_bubble, opts.kmer_r2c, opts.fill_gaps, opts.unmask_contigs, opts.keep_contigs_connected, opts.gap_len_diff, opts.min_contig_len, opts.min_contig_cvg, opts.max_contig_cvg, opts.insert_size_upper_bound, opts.bubble_coverage, opts.genome_size, opts.ass_visual)
        elif int(opts.kmer_size) > 63:
            cmd = "SOAPdenovo-127mer_v2.0 all -s %s -o %s -K %s -p %s -a %s -d %s -R %s -D %s -M %s -m %s -e %s -E %s -k %s -F %s -u %s -w %s -G %s -L %s -c %s -C %s -b %s -B %s -N %s -V %s" % (config_file, tmp_dir + "/result", opts.kmer_size, opts.ncpu, opts.init_memory_assumption, opts.kmer_freq_cutoff, opts.resolve_repeats, opts.edge_cov_cutoff, opts.merge_level, opts.max_k, opts.weight, opts.merge_clean_bubble, opts.kmer_r2c, opts.fill_gaps, opts.unmask_contigs, opts.keep_contigs_connected, opts.gap_len_diff, opts.min_contig_len, opts.min_contig_cvg, opts.max_contig_cvg, opts.insert_size_upper_bound, opts.bubble_coverage, opts.genome_size, opts.ass_visual)

        #print cmd

    #Perform SOAPdenovo2 analysis
    buffsize = 1048576
    try:
        #To hold standard output from process
        stdout_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #To hold standard error output from process
        stderr_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #Open streams to files
        stderr_hd = open(stderr_file, 'w')
        stdout_hd = open(stdout_file, 'w')

        #Call SOAPdenovo2
        proc = subprocess.Popen(args=cmd, shell=True, cwd=tmp_dir, stdout=stdout_hd, stderr=stderr_hd.fileno())
        returncode = proc.wait()

        stderr_hd.close()
        stdout_hd.close()

        #Get stderr, allowing for case where it's very large
        stderr_hd = open(stderr_file, 'r')
        stderr = ''
        try:
            while True:
                stderr += stderr_hd.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        stderr_hd.close()

        #Read tool stdout into galaxy stdout
        f = open(stderr_file)
        lines = f.readlines()
        for line in lines:
            sys.stdout.write(line)
        f.close()

        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        raise Exception, 'Problem performing SOAPdenovo2 process ' + str(e)

    #Read SOAPdenovo2 results into outputs
    contig_out = open(opts.contig, 'wb')
    file = open(tmp_dir + '/result.contig')
    for line in file:
        contig_out.write(line)
    contig_out.close()

    scafseq_out = open(opts.scafseq, 'wb')
    file = open(tmp_dir + '/result.scafSeq')
    for line in file:
        scafseq_out.write(line)
    scafseq_out.close()

    config_out = open(opts.config, 'wb')
    file = open(config_file)
    for line in file:
        config_out.write(line)
    config_out.close()
    file.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)

    #Check results in output file
    if os.path.getsize(opts.contig) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("Problem with SOAPdenovo2 process")

if __name__ == "__main__":
    main()
