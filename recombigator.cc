#include "recombigator.h"

const char* info = "\nrecombigator version 2.0\nPeter Reifenstein June 2024\n            _____        ____\n           / (^) \\______/ (^) \\ \n          /      ________      \\ \n/\\__/\\____|__o___\\______/__o___|____/\\/\\/\\ \n\nA tool to classify Nanopore or HiFi DNA reads as belonging to Genome A, Genome B, or Recombinant\n\ncompile:\ng++ recombigator.cc -std=c++11 -O3 -o recombigator\n\nusage: \n./recombigator alignment_A.sam alignment_B.sam variations.tsv [output_name]\n\n-- alignment_A.sam\nincludes all reads aligned to genome A, with match/mismatch information in the CIGAR string\nfor example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):\n./minimap2 --eqx -a --secondary=no genomeA.fasta reads.fasta\n\n-- variations.tsv\na tab seperated value file with the format specified by SYRI (https://schneebergerlab.github.io/syri/fileformat.html):\n  genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nregion_type beginning with \"SYN\" and variation_type = \"SNP\" indicates a SNP in a syntenic region\ngenomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\nfields marked with ~ are not used and can be any string\n";

bool*** SNP_coordinates;
int num_chromosomes;

ofstream*** output_SAM_files;
ofstream** output_TXT_files;
ofstream summary_file;
ofstream log_file;
ofstream recombinantinfo_file;

ifstream sam_files[2];

map<string, int> chromosome_numbers;
map<int, string> chromosome_names[2];

static string SNP_symbols[] = {"A","B","!","?"};
static string SNP_meanings[] = {"GenotypeA","GenotypeB","BothMismatch","NeitherMismatch"};

static string CATEGORY_NAMES[] = {"GenotypeA", "GenotypeB", "LowQualityRecombinant", "Recombinant","ChompBadAlign", "ChompNoSNPs", "ChompBadGT"};
static int NUM_CATEGORIES = 7;
static int CHOMP_CATEGORIES_INDEX = 3;
int num_reads_in_category[7];

int main(int argc, char* argv[]) {
    if(argc < 4) {
        printProgramInfo();
    }
    string dir_string = "recombigator_output";
    if(argc > 4) {
        string append = argv[4];
        dir_string += "_on_" + append;
    }

    makeDirectories(dir_string);
    readChromosomes(argv[1], argv[2]);
    readSNPs(argv[3]);

    auto begin = chrono::high_resolution_clock::now();
    analyzeReads();
    cout << "[recombigator] Done analyzing reads" << endl;

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    printCategorizationResults((double)elapsed / (double) 1000000);

    memoryCleanup();
}

void categorizeRead(int categorization, string& read_name, string& sam_line_A, string& sam_line_B, string& textoutput, bool silenced) {

    // 1. print textoutput to standard output, if not silenced
    textoutput.append(" ***** Read is " + CATEGORY_NAMES[categorization] + "\n");
    if(!silenced) cout << textoutput << endl;

    // 2. increment counter based on categorization
    num_reads_in_category[categorization]++;

    // 3. write the SAM files with the read
    *output_SAM_files[0][categorization] << read_name << sam_line_A << endl;
    *output_SAM_files[1][categorization] << read_name << sam_line_B << endl;

    // 4. write the text file with the read and the all reads text file
    *output_TXT_files[categorization] << textoutput << endl;
    *output_TXT_files[NUM_CATEGORIES] << textoutput << endl;

    return;
}

void analyzeReads() {
    
    // the header fields are already consumed when finding the chromosome sizes

    cout << "[recombigator] Analyzing reads" << endl;

    string current_read_name = "";
    string next_read_name = "";

    sam_files[0] >> next_read_name;
    sam_files[1] >> next_read_name;

    int read_number = 0;
    
    while(sam_files[0]) {

        // 1. get mismatch locations for both alignments
        vector<int> mmlocations[2];
        int align_len[2];
        unsigned long align_loc[2];
        unsigned long align_end_loc[2];
        int align_chrom[2];
        string align_chrom_name[2];
        int n_alternate_alignments[2];
        string SAM_line[2];

        vector<int> expectedSNPlocations[2];
        vector<bool> genotypeMatrix[2];
        vector<int> concensusMatrix;
        int num_gt_switches = -1;

        string textoutput;

        bool silenced = true;
        if(read_number % 100 == 0)
            silenced = false;
        read_number++;

        current_read_name = next_read_name;
        textoutput.append("#### [" + to_string(read_number) + "] Processing Read " + current_read_name + '\n');

        string line;
        for(int i = 0; i < 2; i++) {
            getline(sam_files[i], SAM_line[i]);
            stringstream read(SAM_line[i]);

            parseCIGAR(read, (mmlocations[i]), align_len[i], align_loc[i], align_chrom_name[i]);
            align_end_loc[i] = align_loc[i] + align_len[i];
            align_chrom[i] = chromosome_numbers[align_chrom_name[i]];

            n_alternate_alignments[i] = 0;
            sam_files[i] >> next_read_name;
            while(next_read_name == current_read_name && sam_files[i]) { // stop looking for alternates if a new read name comes up or end of file reached
                // alternate alignment
                n_alternate_alignments[i]++;
                getline(sam_files[i], line);
                sam_files[i] >> next_read_name;
            }
        }

        // 2. create SNP matrix for both alignments
        for(int i = 0; i < 2; ++i) {
            auto mismatch_location_it = mmlocations[i].begin();
            // stores the last checked index for a mismatch
            //cout << "Finding SNP matrix for i = " << i << " and align location = " << align_loc[i] << " and end location = " << align_end_loc[i] << " chromosome of " << align_chrom[i] << endl;
            for(unsigned long pos = align_loc[i]; pos < align_end_loc[i]; ++pos) {
                if(SNP_coordinates[i][align_chrom[i]][pos] == 1) {
                    int relative_align_pos = pos - align_loc[i];
                    // add this to the list of SNP locations
                    expectedSNPlocations[i].push_back(relative_align_pos);
                    bool SNP_exists = false;
                    // linear scan the mismatch locations
                    while(true) {
                        if(mismatch_location_it == mmlocations[i].end() || *mismatch_location_it > relative_align_pos) {
                            break;
                        }
                        //cout << "checking mismatch location if it matches SNP location, mm location is: " << *mismatch_location_it << endl;
                        if(*mismatch_location_it == relative_align_pos) {
                            SNP_exists = true;
                            break;
                        }
                        mismatch_location_it++;
                    }
                    //cout << "Finished" << endl;
                    if(SNP_exists)
                        genotypeMatrix[i].push_back(i ^ true);
                    else
                        genotypeMatrix[i].push_back(i);
                }
            }
        }
        // TO TEST / DEBUG

        // 3. create consensus SNP matrix
        vector<bool>::iterator matrix_it[2];
        int numSNPs = 0;
        int numIgnoredSNPs = 0;
        for(int i = 0; i < 2; i++) matrix_it[i] = genotypeMatrix[i].begin();
        int previous_gt = -1, breakpoint = -1, index = 0, start_gt = -1;
        while(matrix_it[0] != genotypeMatrix[0].end() && matrix_it[1] != genotypeMatrix[1].end()) {

            if(*matrix_it[0] == *matrix_it[1]) {
                // both agree on categorization
                int SNP_categorization = *matrix_it[0];
                concensusMatrix.push_back(SNP_categorization);
                if(start_gt == -1)
                    start_gt = SNP_categorization;
                if(previous_gt != SNP_categorization) {
                    previous_gt = SNP_categorization;
                    num_gt_switches++;
                    breakpoint = numSNPs;
                }
                ++numSNPs;
            } else if(*matrix_it[0] == 1 && *matrix_it[1] == 0) {
                // both have a mismatch, ignore this SNP
                ++numIgnoredSNPs;
                concensusMatrix.push_back(2);
            } else {
                // neither have a mismatch, ignore this SNP
                ++numIgnoredSNPs;
                concensusMatrix.push_back(3);
            }
            ++matrix_it[0]; ++ matrix_it[1]; ++ index;
        }

        // Write to file
        for(int i = 0; i < 2; i++) {
            textoutput.append("  >> Aligned to chromosome " + to_string(align_chrom[i] + 1)  + " (called " + align_chrom_name[i] + ") from position " + to_string(align_loc[i]) + " to " + to_string(align_end_loc[i]) + " (Total Length " + to_string(align_len[i]) + ")");
            if(n_alternate_alignments[i] > 0)
                textoutput.append(" Note: " + to_string(n_alternate_alignments[i]) + " alternate alignment(s)");
            textoutput.append("\n");

            textoutput.append("Mismatch Locations     (Total " + to_string(mmlocations[i].size()) + "):");
            for(auto it = mmlocations[i].begin(); it != mmlocations[i].end(); ++it) {
                textoutput.append('\t' + to_string(*it));
            } textoutput.append("\n");

            textoutput.append("Expected SNP Locations (Total " + to_string(expectedSNPlocations[i].size()) + "):");
            for(auto it = expectedSNPlocations[i].begin(); it != expectedSNPlocations[i].end(); ++it) {
               textoutput.append('\t' + to_string(*it));
            } textoutput.append("\n");

            textoutput.append("Genotype Matrix:");
            for(auto it = genotypeMatrix[i].begin(); it != genotypeMatrix[i].end(); ++it) {
                textoutput.append("\t[" + (SNP_symbols[*it]) + "]");
            } textoutput.append("\n");
        }
        textoutput.append(" **** Concensus Matrix ");
        if(genotypeMatrix[0].size() != genotypeMatrix[1].size()) {
            log_file << "[anomolous case] read " << current_read_name << " has genotype matrixes with different number of SNPs" << endl;
            textoutput.append("(Warning: different GT Matrix sizes)");
        }
        textoutput.append(":");
        for(auto it = concensusMatrix.begin(); it != concensusMatrix.end(); ++it) {
            textoutput.append("\t" + (SNP_symbols[*it]) + "");
        } textoutput.append("\n");

        // 4. categorization
        // #CHOMP CATEGORY 4
        if(align_chrom[0] != align_chrom[1]) {
            categorizeRead(4, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
            continue;
        }
        // #CHOMP CATEGORY 5
        if(numSNPs < 1) {
            categorizeRead(5, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
            continue;
        }
        // #CHOMP CATEGORY 6
        if(num_gt_switches > 1) {
            categorizeRead(6, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
            continue;
        }

        // #GENOTYPE A or GENOTYPE B
        if(num_gt_switches == 0) {
            categorizeRead(start_gt, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
            continue;
        }
        
        // #LOW QUALITY RECOMBINANT
        // TO COMPLETE
        if(breakpoint == 1 || breakpoint == numSNPs - 1) {
            categorizeRead(2, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
            continue;
        }

        // #RECOMBINANT
        double relativeBreakpoint = (double) breakpoint / (double) numSNPs;
        // write to recombinant info
        // name numSNPs(A:B)    relativeBreakpoint  numIgnoredSNPs  numAlternateAlignments(A:B)   firstSNP    alignDir(A:B)   chromosome(A:B)   len(A:B)    alignloc(A:B)
        recombinantinfo_file << current_read_name << '\t' << genotypeMatrix[0].size() << ':' << genotypeMatrix[1].size() << relativeBreakpoint
            << numIgnoredSNPs << '\t' << n_alternate_alignments[0].size() << '\t' << n_alternate_alignments[1].size()

        categorizeRead(3, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced);
        
    }

    return;
}

void parseCIGAR(stringstream& read, vector<int>& mm_locs, int& len, unsigned long& loc, string& chrom_name) {
    // parse the CIGAR string of the next read in the sam file, and return the name of the next read
    
    // Lsat_Sal__Chrom2_(7016910..7020531)	2064	SERH_U_Chr_2	35850229	1	4=1X82=1X79=1X60=1X111=1X20=1X3=1X97=1X194=1X12=1X28=1X12=1X21=1X15=1X3=1I49=1X10=1X1=1X144=1D76=1X9=1X1=1X5=2565H
    // readname                             flag    chrom_aligned   position    qual    CIGARstring
    
    int num;
    char token;

    read >> num;
    read >> chrom_name;
    read >> loc;
    read >> num;

    int pos = 0; // current position on reference

    // cout << "Parsing CIGAR" << endl;
    while(read >> num) {
        read >> token;
        // cout << num << token << endl;
        if(token == '=' || token == 'D') { // = match // D deletion
            pos += num;
        } else if(token == 'X') {
            for(int i = 0; i < num; i++) {
                mm_locs.push_back(pos);
                // cout << "Found mismatch at position " << pos << endl;
                pos++;
            }
        } else {
            // ignore
        }
        
        // X mismatch
        // H hard clip (ignore, does not consume reference)
        // S soft clip (ignore, does not consume reference)
        // D deletion in reference (consume the reference)
        // I insertion (ignore, does not consume reference)
    }
    // cout << "Finished" << endl;
    len = pos;
    return;
}

void printProgramInfo() {
    cout << info << endl;
    exit(1);
}

void printCategorizationResults(double time) {
    stringstream summary;
    long total_reads = 0;
    long total_discarded_reads = 0;

    summary << left;
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        if(i > CHOMP_CATEGORIES_INDEX) {
            if(i == CHOMP_CATEGORIES_INDEX + 1)
                summary << endl;
            total_discarded_reads += num_reads_in_category[i];
            summary << "Total " << setw(16) << (CATEGORY_NAMES[i] + ":") << num_reads_in_category[i] << endl;
        } else {
            summary << "Total " << setw(16) << (CATEGORY_NAMES[i] + ":") << num_reads_in_category[i] << endl;
        }
        total_reads += num_reads_in_category[i];
    }

    double attrition_rate = (double)total_discarded_reads / (double)total_reads;

    summary << endl << "Processed [" << total_reads << "] reads, " << 100 * attrition_rate << "% were thrown out (chomped)" << endl;

    
    summary << "Elapsed time: " << time << " seconds" << endl;

    string summary_str = summary.str();
    cout << summary_str << endl;
    summary_file << summary_str << endl;

    return;
}

void makeDirectories(string base_dir) {
    mkdir((base_dir).c_str(), 0777);

    string alignments_dir = base_dir + "/sorted_alignments";
	mkdir((alignments_dir).c_str(), 0777);

    // open ofstreams for writing SAM files
    output_SAM_files = (ofstream***) malloc(2 * sizeof(ofstream**));

    for(int i = 0; i < 2; i++) {
        string SAM_names_prefix = "toGenome";
        SAM_names_prefix += SNP_symbols[i];
        output_SAM_files[i] = (ofstream**) malloc(NUM_CATEGORIES * sizeof(ofstream*));
        string gt_dir = alignments_dir + "/" + SAM_names_prefix;
        mkdir((gt_dir).c_str(), 0777);
        for(int j = 0; j < NUM_CATEGORIES; j++) {
            output_SAM_files[i][j] = new ofstream(gt_dir + "/" + SAM_names_prefix + "_" + CATEGORY_NAMES[j] + ".SAM");
        }
    }

    // open ofstreams for writing to text files
    string txt_dir = base_dir + "/text_files";
	mkdir((txt_dir).c_str(), 0777);
    output_TXT_files = (ofstream**) malloc((NUM_CATEGORIES + 1) * sizeof(ofstream*));
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        output_TXT_files[i] = new ofstream(txt_dir + "/" + CATEGORY_NAMES[i] + ".txt");
    }
    output_TXT_files[NUM_CATEGORIES] = new ofstream(txt_dir + "/" + "all" + ".txt");

    // open ofstream for writing summary
    summary_file.open(base_dir + "/summary.txt");

    // open ofstream for writing log
    log_file.open(base_dir + "/log.txt");

    // open ofstream for writing recombinant info
    recombinantinfo_file.open(base_dir + "/recombinantinfo.tsv");
    recombinantinfo_file << "name numSNPs(A:B)    relativeBreakpoint  numIgnoredSNPs(A:B)  numAlternateAlignments(A:B)   firstSNP    alignDir(A:B)   chromosome(A:B)   len(A:B)    alignloc(A:B)" << endl;

    return;
}

void createSNPfile(string polymorphism_file) {
    ifstream p_file(polymorphism_file);
    ofstream p_SNP_file(polymorphism_file + ".recombigator.SNPs");
    ofstream p_info_file(polymorphism_file + ".recombigator.info");
    ofstream p_VCF_files[2];
    for(int i = 0; i < 2; i++) {
        p_VCF_files[i].open(polymorphism_file + ".recombigator.SNPs" + ".genome" + SNP_symbols[i] + ".VCF");
        p_VCF_files[i] << "##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
    }
    p_SNP_file << "chromosome_num\tgenome_A_location\tgenome_B_location" << endl;

    if(!p_file) { cout << "[ERROR] cannot open " << polymorphism_file << endl; exit(1); }
    string chrBloc_str, unique_id, parent_id, type, prev_parent_id, str;
    string seq[2];
    string chrName[2];
    //int chrAloc, chrBloc, chrAnum, chrBnum, chrAname, chrBname;
    int chrLoc[2];
    int chrNum[2];
    unsigned long totalcount = 0;
    int num_SNPs_on_chrom [num_chromosomes] = {0};
    int num_syntenic_regions_on_chrom [num_chromosomes] = {0};
    while(p_file >> chrName[0]) {
        p_file >> chrLoc[0]; p_file >> chrLoc[0]; p_file >> seq[0]; p_file >> seq[1];
        p_file >> chrName[1]; p_file >> chrBloc_str;  p_file >> chrBloc_str;
        p_file >> unique_id; p_file >> parent_id; p_file >> type; p_file >> str;

        // SALH_C_Chr_1	1	26851	-	-	-	-	-	NOTAL1	-	NOTAL	-
        if(parent_id.substr(0, 3) == "SYN" && type == "SNP") {
            // found SNP in syntenic region
            chrLoc[1] = stoi(chrBloc_str);
            totalcount++;
            chrNum[0] = chromosome_numbers[chrName[0]];
            chrNum[1] = chromosome_numbers[chrName[1]];
            num_SNPs_on_chrom[chrNum[0]] ++;
            if(chrNum[0] != chrNum[1]) {
                cout << "[warning] Differnent chromosomes in 'syntenic region'" << endl;
            }

            // write to SNP file
            p_SNP_file << chrNum[0] + 1 << '\t' << chrLoc[0] << '\t' << chrLoc[1] << endl;

            // write to both VCF files
            for(int i = 0; i < 2; i++) {
                // CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                p_VCF_files[i] << (chromosome_names[i][chrNum[i]] + "\t" + to_string(chrLoc[i]) + "\t" + unique_id + "\t" + seq[i] + "\t" + seq[i ^ true] + "\t.\tPASS\t." ) << endl;
            }

            if(parent_id != prev_parent_id) {
                num_syntenic_regions_on_chrom[chrNum[0]]++;
            }
            prev_parent_id = parent_id;
            //if(totalcount % 400000 == 0)
            //    cout << "[SNP parse] Found SNP in syntenic region on chromosome " << chrAnum << " A coordinate: " << chrAloc << " B coordinate: " << chrBloc << endl;
        }
    }
    for(int i = 0; i < num_chromosomes; i++) {
        p_info_file << "Chromosome " << i + 1 << " has total " << num_SNPs_on_chrom[i] << " SNPs among " << num_syntenic_regions_on_chrom[i] << " syntenic regions" << endl;
    }
    
}

void readSNPs(string polymorphism_file) {

    // Lsat_PI251246_v14_Chr1  3435181 3435181 C       T       Lsat_PI251246_v14_Chr4  360187218       360187218       SNP13644        INVTR7135       SNP     -
    // ref_chrom               ref_pos ref_pos ref_seq qry_seq qry_chrom               qry_pos         qry_pos         unique_id       parent_id       type    copy-status
    //genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nwhere the region_type beginning with "SYN" indicates a syntenic region and variation_type = "SNP" indicates a SNP\nand genomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\n(the other fields marked with ~ are not used and can be any string)
    //
    //
    // (the other fields marked with ~ are not used and can be any string)
    // look for parent_id = SNYx (i.e. the variation occurs in the xth syntenic alignment region) and type = SNP   

    cout << "[recombigator] Reading in Polymorphism file" << endl;

    ifstream p_SNP_file_test(polymorphism_file + ".recombigator.SNPs");
    if(!p_SNP_file_test) {
        cout << "[SNP parse] No SNP file created yet from polymorphisms file, creating SNP file and VCF files..." << endl;
        createSNPfile(polymorphism_file);
    }

    cout << "[recombigator] Reading in SNP file" << endl;
    ifstream p_SNP_file(polymorphism_file + ".recombigator.SNPs");
    string header;
    getline(p_SNP_file, header);
    unsigned long chrAloc, chrBloc;
    int chromosome;
    int count = 0;
    while(p_SNP_file >> chromosome >> chrAloc >> chrBloc) {
        SNP_coordinates[0][chromosome - 1][chrAloc] = 1;
        SNP_coordinates[1][chromosome - 1][chrBloc] = 1;
        count ++;
    }
    //cout << "[SNP parse] Done reading in SNP file" << endl;
    
}

void readChromosomes(string alignment_genome_Afile, string alignment_genome_Bfile) {
    sam_files[0].open(alignment_genome_Afile);
    sam_files[1].open(alignment_genome_Bfile);

    if(!sam_files[0] || !sam_files[1]) { cout << "[ERROR] cannot open one or more SAM files" << endl; exit(1); }

    cout << "[recombigator] Reading chromosomes from SAM files" << endl;
    vector<int> genome_sizes[2];
    string chrom_name;

    char c;
    string str;
    int len;
    int n = 0;
    for(int i = 0; i < 2; i++) {
        n = 0;
        while(true) {
            sam_files[i] >> c;
            if(c!= '@') {
                sam_files[i].unget();
                break;
            }
            
            getline(sam_files[i], str);
            str = "@" + str;
            stringstream header_line(str);
            header_line << str;
            for(int k = 0; k < NUM_CATEGORIES; k++) {
                *output_SAM_files[i][k] << str << endl;
            }

            header_line >> str;
            if(str == "@SQ") { // genomic alignment reference sequence (chromosome)
                header_line >> str;
                str = str.substr(3);
                header_line >> c;header_line >> c;header_line >> c;
                header_line >> len;

                log_file << "[chromosome] " << str << " len: " << len << endl;
                chromosome_numbers[str] = n;
                chromosome_names[i][n] = str;
                genome_sizes[i].push_back(len);
                n++;
            } else { // consume the rest of the field
                getline(header_line, str); 
            }
            
        }
    }
    num_chromosomes = genome_sizes[1].size();

    cout << "[recombigator] Allocating Memory for SNP locations" << endl;

    SNP_coordinates    = (bool***) malloc(2 * sizeof(bool**));
    SNP_coordinates[0] = (bool**)  malloc(num_chromosomes * sizeof(bool*));
    SNP_coordinates[1] = (bool**)  malloc(num_chromosomes * sizeof(bool*));

    for(int i = 0; i < num_chromosomes; i++) {
        // memory allocated with calloc will be zero initialized
        SNP_coordinates[0][i] = (bool*) calloc((genome_sizes[0][i] + 1), sizeof(bool));
        SNP_coordinates[1][i] = (bool*) calloc((genome_sizes[1][i] + 1), sizeof(bool));
        log_file << "[calloc] Allocated memory for chromosome " << i + 1 << " SNPs (total " << genome_sizes[0][i] + 1 << " bytes)" << endl;
    }
    return;
}

void memoryCleanup() {
    log_file << "[recombigator] Deallocating heap memory" << endl;
    for(int i = 0; i < num_chromosomes; i++) {
        free(SNP_coordinates[0][i]);
        free(SNP_coordinates[1][i]);
    }
    free(SNP_coordinates[0]);
    free(SNP_coordinates[1]);
    free(SNP_coordinates);

    for(int i = 0; i < NUM_CATEGORIES; i++) {
        free(output_SAM_files[0][i]);
        free(output_SAM_files[1][i]);
    }
    free(output_SAM_files[0]);
    free(output_SAM_files[1]);
    free(output_SAM_files);

    for(int i = 0; i < NUM_CATEGORIES + 1; i++) {
        free(output_TXT_files[i]);
    }
    free(output_TXT_files);

    log_file << "[recombigator] Deallocated heap memory" << endl;
    return;
}