"""
gapcloser-1.12.py
A wrapper script for GapCloser version 1.12
Peter Li - GigaScience/BGI-HK
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
    parser.add_option("-a", "--scaff_in", dest="scaff_in")

    parser.add_option("-l", "--max_read_length", dest="max_read_length")
    parser.add_option("", "--file_source", dest="file_source")
    parser.add_option("", "--configuration", dest="configuration")

    #Make list of params
    parser.add_option("", "--min_ins", action="append", type="string", dest="min_insert_list")
    parser.add_option("", "--avg_ins", action="append", type="string", dest="avg_insert_list")
    parser.add_option("", "--max_ins", action="append", type="string", dest="max_insert_list")
    parser.add_option("", "--asm_flags", action="append", type="string", dest="asm_flags_list")
    parser.add_option("", "--rank", action="append", type="string", dest="rank_list")

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
    parser.add_option("-p", "--overlap_param", dest="overlap_param")
    parser.add_option("-t", "--thread_num", dest="thread_num")

    #Outputs
    parser.add_option("", "--scaff", dest='scaff')
    parser.add_option("", "--fill_info", dest='fill_info')
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    dirpath = tempfile.mkdtemp(prefix="tmp-gapcloser-")

    if opts.file_source == "history":
        config_file = opts.configuration
    else:
        #Create temp file for configuration
        config_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        try:
            fout = open(config_file,'wb')
            fout.write("max_rd_len=%s\n" % opts.max_read_length)
            #Calculate how sets of data there are - use avg_ins as a measure of this
            #Loop thru this number of times
            #Also use separate index to keep count of reads
            single_read_index = 0
            paired_read_index = 0
            for index in range(len(opts.avg_insert_list)):
                #print "single_read_index ", single_read_index
                #print "paired_read_index ", paired_read_index
                fout.write("[LIB]\n")
                fout.write("min_ins=%s\n" % (opts.min_insert_list)[index])
                fout.write("avg_ins=%s\n" % (opts.avg_insert_list)[index])
                fout.write("max_ins=%s\n" % (opts.max_insert_list)[index])
                fout.write("asm_flags=%s\n" % opts.asm_flags_list[index])
                fout.write("rank=%s\n" % opts.rank_list[index])
                #Add data file configuration - needs careful looping due to single and paired reads
                #print opts.type_of_data_list[index]
                #print opts.format_of_data_list[index]
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
            stop_err("Config file cannot be opened for writing" + str(e))

    #Set up command line call - assumes path to executable has been defined in user's environment
    if opts.default_full_settings_type == "default":
        cmd = "GapCloser -a %s -b %s -o %s" % (opts.scaff_in, config_file, dirpath + "/gapclo.out")
    elif opts.default_full_settings_type == "full":
        cmd = "GapCloser -a %s -b %s -o %s -p %s -t %s" % (opts.scaff_in, config_file, dirpath + "/gapclo.out", opts.overlap_param, opts.thread_num)

    #print cmd

    #Perform GapCloser analysis
    buffsize = 1048576
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stdout = open(tmp_out_file, 'w')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=dirpath).name #Stdout is outputted to here
        tmp_stderr = open(tmp_err_file, 'w')

        #Call SOAPdenovo2
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stderr=tmp_stderr.fileno())
        returncode = proc.wait()

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
        tmp_stderr.close()

        #A return code of 1 is given but result is ok - wierd
        #if returncode > 1:
        if returncode != 0:
            raise Exception, stderr

    except Exception, e:
        raise Exception, 'Problem executing GapCloser ' + str(e)

    #Read results into outputs
    scaff_out = open(opts.scaff, 'wb')
    file = open(dirpath + '/gapclo.out')
    for line in file:
        scaff_out.write(line)
    scaff_out.close()

    fill_info_out = open(opts.fill_info, 'wb')
    file = open(dirpath + '/gapclo.out.fill')
    for line in file:
        fill_info_out.write(line)
    fill_info_out.close()

    #Clean up temp files
    cleanup_before_exit(dirpath)
    #Check results in output file
    if os.path.getsize(opts.scaff) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("Problem with GapCloser process")

if __name__ == "__main__":
    main()
