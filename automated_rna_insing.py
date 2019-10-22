#This is to automate the RNA pipeline 
#Assuming Ensembl formatted reference files 
#Assuming only 1 R1 and 1 R2 file per sample
#Assuming that the data is paired end only 

#Need to import a file that holds all the information that I'll be asking and subprocess to use bash commands from within python
import sys
import subprocess
import time
start_time = time.time()
#data file should include the following: paths to their fastp, pardre, kallisto programs and their fastq files 
#need to decide if I need to ask for the program locations, would make more sense to have a singularity so one doesn't have to worry about version incompatabilities 

#Will require input file, the file will have format such as
#kal_reference: \t path/to/kallisto/reference/file
#after all the reference files and paths have been added, the rest of the file will be as follows:
#path/to/R1/file \t path/to/R2/file \t condition 	and this format repeats until the end of the file for as many samples are provided 
#Need to create a script that will take the input name and make a table of the name and its condition. will output the file to out_dir/name_condition_table.txt

pyther = "/PyRcodes"
data_file = sys.argv[1]
inp_list = []
cond = []
nam_cond = []
kal_tpm_name = []
#Using tic to check if conversion file gets created later
tic = 0 
data = open(data_file, 'r')
for line in data:
	col = line.find(":")
	par = line[:col]
	parr = par.lower()
#The following set of if else statments will separate the lines into their inputs based on the information before the colon, then to separate, clean, and get the file path and or naming for future commands
#Check for comments 
	if(parr.find("#") > -1):
		continue
#Check for empty lines
	elif(parr == ""):
		continue
#Check for Kallisto reference file or an already created index file
	elif(parr.find("kallisto reference") > -1):
		kar_file_check = True
		kalrp = line.strip("\n")
		kalrc = kalrp.split(":")
		kalre = kalrc[1]
		kalrr = kalre.replace(" ","")
		kalr = kalrr.replace("\t", "")
		kalnamm = kalr.split("/")
		kalnama = kalnamm[-1]
#If there is already an index file it will have .idx extension and will be used, otherwise the provided file will be used to make an index file
		if(kalr.find(".idx") > -1):
			ind_file_check = True
		else:
			ind_file_check = False
			kalnam = kalnama.replace(".gz", "")
		print("kallisto reference file found")
#If there already is a conversion file that has the gene to transcript conversion having the first column the transcripts and the second column the genes then it will be used. Otherwise the provided reference file will be used to make a conversion file
	elif(parr.find("conversion file") > -1):
		tic = 1
		con = line.strip("\n")
		conc = con.split(":")
		cone = conc[1]
		convv = cone.replace(" ", "")
		conv = convv.replace("\t", "")
		conv_file_check = True
		conver = conv.split("/")
		conv_nam = conver[-1]
		print("conversion file found")
#Check for the output path that the folders and subsequent files will be outputted to
	elif(parr.find("output") > -1):
		out_file_check = True
		outpp = line.strip("\n")
		outpc = outpp.split(":")
		outpe = outpc[1]
		outpp = outpe.replace(" ","")
		outp = outpp.replace("\t", "")
		print("output folder path found")
#Check for input files, adding the R1 and R2 file paths to a list to be used and a condition that will be saved separately. All these items should be tab separated
	elif(parr.find("input") > -1):
		inp_file_check = True
		linen = line.strip("\n")
		linec = linen.split(":")
		linee = linec[1]
		linet = linee.strip()
		linett = linet.split()
		inp_list.append([linett[0], linett[1]])
		cond.append(linett[2])
	else:
		sys.exit("Error, " + line + " is in unknown format")
#This if else statement checks if a conversion file was created. If it was then tic will equal 1 and pass to the next section. If the file wasn't provided then tic would equal 0 and the bollean value will be set to False 
if (tic == 1):
	pass
else:
	conv_file_check = False

#Check to see if either the reference file was provided or both the index file and conversion file was supplied 
if((ind_file_check and conv_file_check) or ((not ind_file_check) and (not conv_file_check) and kar_file_check)):
	pass
else:
	sys.exit("Error: can only provide the reference file or both the kallisto .idx and the gene transcript conversion file. providing only 1 of these two will result in this error.")
#Take the output path that was provided and making the folders to hold all the output files 
print("Generating output folder")

out_data = outp
if(out_data[-1] == "/"):
	out_dir = out_data + "kallisto_pipe_out"
else:
	out_dir = out_data + "/kallisto_pipe_out"
subprocess.call(["mkdir", out_dir])

dat_out_dir = out_dir + "/Input_data"
subprocess.call(["mkdir", dat_out_dir])

fas_out_dir = out_dir + "/fastp_out"
subprocess.call(["mkdir", fas_out_dir])

par_out_dir = out_dir + "/ParDRe_out"
subprocess.call(["mkdir", par_out_dir])

kal_out_dir = out_dir + "/Kallisto_out"
subprocess.call(["mkdir", kal_out_dir])

des_out_dir = out_dir + "/DESeq2_out"
subprocess.call(["mkdir", des_out_dir])

fq_out_dir = out_dir + "/fastq_size_out"
subprocess.call(["mkdir", fq_out_dir])

#Check to see if index file was provided, if it was then softlink the index file to the output folder. If it wasn't then use the reference file to create one
if(ind_file_check):
	kal_ind = kal_out_dir + "/" + kalnama
else:
	kal_ind = kal_out_dir + "/" + kalnam + "_kal_index.idx"

	print("Generating Kallisto index file")
	#run kallisto to generate the index file
	#command is path/to/kallisto "index" "-i" <,output_file_name> <kallisto reference file>
	subprocess.call(["kallisto", "index", "-i", kal_ind, kalr])

#softlink the reference file to the pipeline folder
subprocess.call(["ln", "-s", kalr, kal_out_dir])

#Check if there is a provide conversion file, if not then one will be created from the provided reference file 
if(conv_file_check):
	tra_gen_tab = out_dir + "/" + conv_nam
	subprocess.call(["ln", "-s", conv, tra_gen_tab])
else:
	print("Generating Gene Transcript conversion table")
#Run tra_gene_table_creator.py script to make the transcript gene conversion file
	tra_gen_tab = out_dir + "/" + kalnam + "_transcript_gene_conversion_table.txt"
	tra_gen_exec = (pyther + "/tra_gene_table_creator.py")
	#tra_gen_exec = (prog_path + "/tra_gene_table_creator.py")
	subprocess.call(["python", tra_gen_exec, "-i", kalr, "-o", tra_gen_tab])
        conv = tra_gen_tab
#Run gene_isoform_checker.py to take the conversion table and make a file that holds all the transcripts that are associated with a single gene, and will also create a histogram of number of transcripts per gene
assoc_tra_exec = pyther + "/gene_isoform_checker.py"
subprocess.call(["python", assoc_tra_exec, conv, out_dir, pyther])


#Taking the list of input file paths and saving the R1 and R2 as variables
for x in range(0, len(inp_list)):
	if(inp_list[x][0].find("R1") > -1 and inp_list[x][1].find("R2")):
		input_R1 = inp_list[x][0]
		input_R2 = inp_list[x][1]
	elif(inp_list[x][0].find("R2") > -1 and inp_list[x][1].find("R1")): 
		input_R1 = inp_list[x][1]
		input_R2 = inp_list[x][0]
	else: 
		sys.exit("Error, R1 and or R2 not found in input names for paired end samples")
#Taking the absolute path out of the inputs and leaving just the file names so the later files can have the same names
	in1sep = input_R1.split("/")
	in1 = in1sep[(len(in1sep) -1)]
	in2sep = input_R2.split("/")
	in2 = in2sep[(len(in2sep) -1)]
	
#Takes the input files and softlinks them to the output folder
	in1sl = (dat_out_dir + "/" + in1)
	in2sl = (dat_out_dir + "/" + in2)
	subprocess.call(["ln", "-s", input_R1, in1sl])
	subprocess.call(["ln", "-s", input_R2, in2sl])
#Takes the original name, removes the .fastq.gz ending and replaces the R1 with nothing to make a folder name for future commands and graphs 
	fol_name = in1[:-9]
	hold = fol_name.find("R1")
#If the file has R1.fastq.gz then this will be used to skip the elif which would call an out of index error
	if(len(fol_name)-1 == hold +1):
            if(fol_name.find("_R1") > -1):
                fol_name1 = fol_name.replace("_R1", "") 
            else:
		fol_name1 = fol_name.replace("R1", "")
	elif(fol_name[(hold-1)] == "_" and fol_name[(hold + 2)] == "_"):
		fol_name1 = fol_name[:hold] + fol_name[hold + 3:]
	else:
            if(fol_name.find("_R1") > -1):
                fol_name1 = fol_name.replace("_R1", "")
            else:
		fol_name1 = fol_name.replace("R1", "")

	fas_overall = fas_out_dir + "/" + fol_name1
	subprocess.call(["mkdir", fas_overall])
	fas_stat_out = fas_overall + "/" + fol_name1
#Adds the name and condition that was saved earlier so that a table of the sample name and condition could be created
	nam_cond.append([fol_name1, cond[x]])

#Creating the output names for the fastp read files and the fastp statistical files 
	output_fastp_R1 = fas_overall + "/" + in1[:-9] + "_fastp.gz"
	out_fas1 = fas_out_dir + "/" + in1[:-9] + "_fastp.gz"
	output_fastp_R2 = fas_overall + "/" + in2[:-9] + "_fastp.gz"
	out_fas2 = fas_out_dir + "/" + in2[:-9] + "_fastp.gz"
	output_fastp_html = fas_stat_out + "_fastp.html"
	output_fastp_json = fas_stat_out + "_fastp.json"
	
#Run the size_dist_counter.py script to calculate the read sizes of the fastq files and make histograms of the read sizes
	sizes_exec = pyther + "/size_dist_counter.py"
	subprocess.call(["python", sizes_exec, fq_out_dir, pyther, input_R1])

	print("Running fastp")

	#run fastp  
	#fastp_exec = (prog_path + "/fastp")
#fastp command is as follows: path/to/fastp "-i" <input R1 file> "-I" <input R2 file> "-o" <output R1 file name> "-O" <output_R2 file name> "-h" <output html file name> "-j" <output json file name>
	subprocess.call(["fastp", "-i", input_R1, "-I", input_R2, "-o", output_fastp_R1, "-O", output_fastp_R2, "-h", output_fastp_html, "-j", output_fastp_json])

#Run the size_dist_counter.py again to make histogram based on the fastp quality checked read files 
	subprocess.call(["python", sizes_exec, fas_overall, pyther, output_fastp_R1])

#Take the fastp output files and use them as ParDRe input files and create the output file names 
	outpar1 = par_out_dir + "/" + in1[:-9] + "_par.gz"
	outpar2 = par_out_dir + "/" + in2[:-9] + "_par.gz"

	print("Running ParDRe")
	#run ParDRe
#ParDRe command is as follows: path/to/ParDRe "-i" <input R1 file> "-p" <input R2 file> "-z" (states that the input is compresssed in .gz format) "-o" <output R1 file> "-r" <output R2 file>
	subprocess.call(["ParDRe", "-i", output_fastp_R1, "-p", output_fastp_R2, "-z", "-o", outpar1, "-r", outpar2])


	print("Running Kallisto quantification")
	#run kallisto quantification step 
	kal_quant_out = kal_out_dir + "/" + fol_name1
#The format for the command are as follows: path/to/kallisto "quant" "-i" <index file> "-o" <output folder to print to> "-b" <number of bootstraps> <input 1> <input 2>
	subprocess.call(["kallisto", "quant", "-i", kal_ind, "-o", kal_quant_out, "-b", "100", outpar1, outpar2])


#Take the Kallisto .tsv output and append the file path to a list to be used to generate a table of tpm's
	kal_tpm_name.append(kal_quant_out + "/abundance.tsv")
#Run kall_tpm_table_creator.py script to take the list of Kallisto output and combine all the tpm's into a single tpm table	
kal_name = " ".join(kal_tpm_name)
kal_tpm_exec = (pyther + "/kall_tpm_table_creator.py")
subprocess.call(["python", kal_tpm_exec, kal_out_dir, pyther, kal_name])


#At this point all the files will be gathered and used together after this point

#Taking the previously created name condition list, this will iterate and print a table that holds the name of the sample and the condition of the sample
namc_file = (out_dir + "/name_condition_table.txt")
nam_file = open(namc_file, "w")
nam_file.write("Name\tCondition\n")
for z in range(0,len(nam_cond)):
	nam_file.write(nam_cond[z][0] + "\t" + nam_cond[z][1] + "\n")
nam_file.close()


#Linux based commands have finished and now importing the information to the connected Rscript to run DESEq2 and generate other tables and graphs
print("Kallisto quantification finished, entering R commands")
R_exec = (pyther + "/rna_kal_automated.R")
subprocess.call(["Rscript", R_exec, out_dir, conv])



tot_time = time.time() - start_time
print("Total run time: " + str(tot_time) + " seconds.")

