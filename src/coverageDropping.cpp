/*
    FRC: computes the FRC curve starting from alignments
    Copyright (C) 2011  F. Vezzi(vezi84@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
 
#include <stdio.h> 
#include <time.h>
#include <string>
#include <vector>
#include "samtools/sam.h"
#include "options/Options.h"

#include "data_structures/Contig.h"


#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)

#define EXIT_IF_NULL(P) \
  if (P == NULL) \
    return 1;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core)
{

  if (core->flag&BAM_FUNMAP) {
    return false;
  }

  return true;
}


/**
 * Open a .sam/.bam file.
 * @returns NULL is open failed.
 */
samfile_t * open_alignment_file(std::string path)
{
  samfile_t * fp = NULL;
  std::string flag = "r";
  if (path.substr(path.size()-3).compare("bam") == 0) {
    //BAM file!
    flag += "b";
  }
  if ((fp = samopen(path.c_str(), flag.c_str() , 0)) == 0) {
    fprintf(stderr, "coverageDropping: Failed to open file %s\n", path.c_str());
  }
  return fp;
}

static int fetch_func(const bam1_t *b, void *data)
{


	const bam1_core_t* core =  &b->core;
	uint32_t* cigar = bam1_cigar(b);
	int32_t alignmentLength = bam_cigar2qlen(core,cigar);
//	cout << core->tid << " " << core->pos <<  " " << alignmentLength << " a\n";

	vector<bam1_t> *buf = (vector<bam1_t>*)data;
	buf->push_back(*b);

	uint32_t size = buf->size();
	const bam1_t* bAF = &buf->at(0);
	const bam1_core_t* coreAF =  &bAF->core;
	uint32_t* cigarAF = bam1_cigar(bAF);
	int32_t alignmentLengthAF = bam_cigar2qlen(coreAF,cigarAF);
	//	cout << coreAF->tid << " " << coreAF->pos <<  " " << alignmentLengthAF << " b\n";


    return 0;
}



void computeLibraryStats(samfile_t *fp, unsigned int minInsert, unsigned int maxInsert, uint64_t genomeLength);


int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string assemblyFile = "";
	unsigned int numLibraries = 0;

	vector<string> libraries;
	vector<unsigned int> minInserts;
	vector<unsigned int> maxInserts;

	unsigned int WINDOW = 1000;
	uint64_t estimatedGenomeSize;

	string outputFile="";

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
			("libraries", po::value< vector < string > >(), "bam files, one for each different library (assumed ordered from the shortest to the longest)")
			("min-insert",  po::value<vector< unsigned int> >(), "minimum allowed insert size, one for each library (assumed ordered from the shortest to the longest). Used in order to filter outliers. Insert size goeas from beginning of first read to end of second read")
			("max-insert",  po::value<vector< unsigned int> >(), "maximum allowed insert size, one for each library (assumed ordered from the shortest to the longest). Used in order to filter outliers. Insert size goeas from beginning of first read to end of second read")
			("window",  po::value<unsigned int>(), "window size for coverage drop computation")
			("output",  po::value<string>(), "Header output file names (default FRC.txt and Features.txt)")
			("genome-size", po::value<unsigned long int>(), "genome size (estimation)")
			;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(2);
	}
	if (vm.count("help")) {
		DEFAULT_CHANNEL << desc << endl;
		exit(0);
	}


	if (!vm.count("libraries")  or !vm.count("min-insert") or !vm.count("max-insert")) {
		DEFAULT_CHANNEL << "At least one library must be specified.\n"
				"libraries, min-insert, and max-insert are all mandatory parameters.\n"
				"Run with --help for more details" << endl;
		exit(0);
	}


	libraries = vm["libraries"].as<vector<string> >();
	numLibraries = libraries.size();
	minInserts = vm["min-insert"].as<vector<unsigned int> >();
	maxInserts = vm["max-insert"].as<vector<unsigned int> >();

	cout << numLibraries << " " << minInserts.size() << " " << maxInserts.size() << "\n";

	if(!(numLibraries == minInserts.size() && numLibraries == maxInserts.size()) ) {
		DEFAULT_CHANNEL << "different number of libraries, min-inserts and max-inserts!!!" << endl;
		exit(0);
	}

	for(unsigned int i = 0; i < libraries.size(); i++) {
		cout << libraries.at(i) << "\t" << minInserts.at(i) <<  "\t" << maxInserts.at(i) << "\n";
	}

	if (vm.count("window")) {
		WINDOW = vm["window"].as<unsigned int>();
	}

	if (vm.count("output")) {
		string header = vm["output"].as<string>();
		outputFile = header;
	}

	if (vm.count("genome-size")) {
		estimatedGenomeSize = vm["genome-size"].as<unsigned long int>();
	} else {
		estimatedGenomeSize = 0;
	}


// NOW start to compute statics
	samfile_t *firstBAMFile;
	firstBAMFile = open_alignment_file(libraries.at(0)); // Open the first file in order to memorize contigs
	EXIT_IF_NULL(firstBAMFile);

	bam_header_t* head = firstBAMFile->header; // sam header
	unsigned int numSequences = head->n_targets;
	unsigned long int genomeLength = 0;

	for(unsigned int i=0; i< numSequences ; i++) {
		genomeLength += head->target_len[i];
	}
	if (estimatedGenomeSize == 0) {
		estimatedGenomeSize = genomeLength; // if user has not provided genome size, approaximated it with assembly size
	}
	samclose(firstBAMFile);

	cout << "number of libraries " << numLibraries << endl;
	cout << "total number of contigs " << numSequences << endl;
	cout << "assembly length " << genomeLength << "\n";
	cout << "estimated length " << estimatedGenomeSize << "\n";


	vector<samfile_t *> librariesBAM;
	vector<const bam_index_t*> librariesBAMindex;

	for(unsigned int i = 0; i< libraries.size() ; i++) {
		string path2bam = libraries.at(i);
		cout << "now opening " << path2bam << " ... ";
		samfile_t *BAMFile;
		BAMFile = open_alignment_file(path2bam); // Open the first file in order to memorize contigs
		EXIT_IF_NULL(BAMFile);
		const bam_index_t *bamIndex;
		bamIndex = bam_index_load(path2bam.c_str());
		EXIT_IF_NULL(bamIndex);
		librariesBAM.push_back(BAMFile);
		librariesBAMindex.push_back(bamIndex);
		cout << "done\n";
	}


	vector<bam1_t> buffer;
	int beg = 0;
	int end = 0x7fffffff;
	int ref;
	unsigned int contigSize;
	//const bam1_t *b = bam_init1();

	//computeLibraryStats(librariesBAM.at(0) , minInserts.at(0), maxInserts.at(0), estimatedGenomeSize);
	for(unsigned int i=0; i< numSequences ; i++) {
			beg = 0;
			end = contigSize = head->target_len[i];
			Contig *currentContig =  new Contig(contigSize, numLibraries);
			for(unsigned int lib = 0; lib < numLibraries; lib++) {
				currentContig->setLibraryLimits(minInserts.at(lib), maxInserts.at(lib), lib);
			}

			bam_parse_region(librariesBAM.at(0)->header, head->target_name[i] , &ref, &beg, &end);
			if (ref < 0) {
				fprintf(stderr, "Invalid region %s\n", head->target_name[i]);
				return 1;
			}


			for(unsigned int lib = 0; lib < numLibraries; lib++) {
				bam_fetch(librariesBAM.at(lib)->x.bam, librariesBAMindex.at(lib) , ref, beg, end, &buffer, fetch_func);
				unsigned int sizeBuffer = buffer.size();
//				cout << "\tnumber of alignments on contig " <<  head->target_name[i] <<  " with library " << lib << " is " << sizeBuffer << "\n";
				for(unsigned int j = 0; j < sizeBuffer; j++ ) {
					const bam1_t* b = &buffer.at(j);
				//	const bam1_core_t* core =  &b->core;
				//	int32_t alignmentLength = core->l_qseq; // bam_cigar2qlen(core,cigar);
				//	cout << core->tid << " " << head->target_name[core->tid] << " " << core->pos <<  " " << alignmentLength << "\n";
					currentContig->updateContig(b, lib);
				}
				buffer.clear();
			}


			currentContig->computeContigStats();
			currentContig->printStats();

			currentContig->computeCoverageDrops();

			//currentContig->print();
			delete currentContig;
		}


}



void computeLibraryStats(samfile_t *fp, unsigned int minInsert, unsigned int maxInsert, uint64_t genomeLength) {
	//Initialize bam entity
	bam1_t *b = bam_init1();

	//All var declarations
	unsigned int contigs = 0; // number of contigs/scaffolds
// total reads
	uint32_t reads = 0;
	uint64_t readsLength = 0;   // total length of reads
	uint32_t unmappedReads = 0;
	uint32_t mappedReads = 0;
	uint32_t zeroQualityReads = 0;
	uint32_t duplicates = 0;

	uint64_t contigSize = 0;
	uint64_t insertsLength = 0; // total inserts length
	float insertStd;

// mated reads (not necessary correctly mated)
	uint32_t matedReads = 0;        // length of reads that align on a contig with the mate
	uint64_t matedReadsLength = 0;  // total length of mated reads
 // correctly aligned mates
	uint32_t correctlyMatedReads = 0; // total number of correctly mated reads
	uint64_t correctlyMatedReadsLength = 0; // length of correctly mated reads
// wrongly oriented reads
	uint32_t wronglyOrientedReads = 0; // number of wrongly oriented reads
	uint64_t wronglyOrientedReadsLength = 0; // length of wrongly oriented reads
// wrongly distance reads
	uint32_t wronglyDistanceReads       = 0; // number of reads at the wrong distance
	uint64_t wronglyDistanceReadsLength = 0;  // total length of reads placed in different contigs
// singletons
	uint32_t singletonReads = 0; // number of singleton reads
	uint64_t singletonReadsLength = 0;     // total length of singleton reads
// mates on different contigs
	uint32_t matedDifferentContig = 0; // number of contig placed in a different contig
	uint64_t matedDifferentContigLength = 0; // total number of reads placed in different contigs

	float C_A = 0; // total read coverage
	float S_A = 0; // total span coverage
	float C_M = 0; // coverage induced by correctly aligned pairs
	float C_W = 0; // coverage induced by wrongly mated pairs
	float C_S = 0; // coverage induced by singletons
	float C_C = 0; // coverage induced by reads with mate on a diferent contif

// compute mean and std on the fly
	float Mk = 0;
	float Qk = 0;
	uint32_t counterK = 1;
//Keep header for further reference
    bam_header_t* head = fp->header;
    int32_t currentTid = -1;
    int32_t iSize;

    while (samread(fp, b) >= 0) {
	      //Get bam core.
	      const bam1_core_t *core = &b->core;
	      if (core == NULL) {  //There is something wrong with the read/file
	    	  printf("Input file is corrupt!");
	    	  return;
	      }
	      ++reads; // otherwise one more read is readed

	      if (!is_mapped(core)) {
	    	  ++unmappedReads;
	      } else {
	    	  if (core->tid != currentTid) {
	    		  //Get length of next section
	    		  contigSize = head->target_len[core->tid];
	    		  contigs++;
	    		  if (contigSize < 1) {//We can't have such sizes! this can't be right
	    			  fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid], contigSize);
	    		  }
	    		  currentTid = core->tid;
	    	  }

	    	  if(!(core->flag&BAM_FUNMAP) && !(core->flag&BAM_FDUP) && !(core->flag&BAM_FSECONDARY) && !(core->flag&BAM_FQCFAIL)) { // if read has been mapped and it is not a DUPLICATE or a SECONDARY alignment
	    		  uint32_t* cigar = bam1_cigar(b);
	    		  ++mappedReads;
	    		  uint32_t alignmentLength = bam_cigar2qlen(core,cigar);
	    		  readsLength += alignmentLength;
	    		  uint32_t startRead = core->pos; // start position on the contig
	    		  uint32_t startPaired;
	    		  //Now check if reads belong to a proper pair: both reads aligned on the same contig at the expected distance and orientation
	    		  if ((core->flag&BAM_FREAD1) //First in pair
	    				  && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
	    				  && (core->tid == core->mtid) /*Mate on the same chromosome*/
	    		  ) {
	    			  //pair is mapped on the same contig and I'm looking the first pair
	    			  startPaired = core->mpos;
	    			  if(startRead < startPaired) {
	    				  iSize = (startPaired + core->l_qseq -1) - startRead; // insert size, I consider both reads of the same length
	    				  if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
	    					  //here reads are correctly oriented
	    					  if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
	    						  if(counterK == 1) {
	    							  Mk = iSize;
	    							  Qk = 0;
	    							  counterK++;
	    						  } else {
	    							  float oldMk = Mk;
	    							  float oldQk = Qk;
	    							  Mk = oldMk + (iSize - oldMk)/counterK;
	    							  Qk = oldQk + (counterK-1)*(iSize - oldMk)*(iSize - oldMk)/(float)counterK;
	    							  counterK++;
	    						  }
	    						  insertsLength += iSize;
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  wronglyDistanceReads++;
	    						  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  //pair is wrongly oriented
	    					  wronglyOrientedReads++;
	    					  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    				  }
	    			  } else {
	    				  iSize = (startRead + alignmentLength - 1) - startPaired;
	    				  if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
	    					  //here reads are correctly oriented
	    					  //here reads are correctly oriented
	    					  if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
	    						  if(counterK == 1) {
	    							  Mk = iSize;
	    							  Qk = 0;
	    							  counterK++;
	    						  } else {
	    							  float oldMk = Mk;
	    							  float oldQk = Qk;
	    							  Mk = oldMk + (iSize - oldMk)/counterK;
	    							  Qk = oldQk + (counterK-1)*(iSize - oldMk)*(iSize - oldMk)/(float)counterK;
	    							  counterK++;
	    						  }
	    						  insertsLength += iSize;
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  wronglyDistanceReads++;
	    						  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  //pair is wrongly oriented
	    					  wronglyOrientedReads++;
	    					  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    				  }
	    			  }
	    		  } else  if ((core->flag&BAM_FREAD2) //Second in pair
	    				  && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
	    				  && (core->tid == core->mtid) /*Mate on the same chromosome*/
	    		  )
	    	// if I'm considering the second read in a pair I must check it is is a correctly mated read and if this is the case update the right variables
	    		  {
	    			  startPaired = core->mpos;
	    			  if(startRead > startPaired) {
	    				  iSize = (startRead + alignmentLength -1) - startPaired;
	    				  if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
	    					  //here reads are correctly oriented
	    					  if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert, no need to update insert coverage
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  } else {
	    						  wronglyDistanceReads++;
	    						  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  //pair is wrongly oriented
	    					  wronglyOrientedReads++;
	    					  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    				  }
	    			  } else {
	    				  iSize = (startPaired + core->l_qseq -1) - startRead;
	    				  if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
	    					  //here reads are correctly oriented
	    					  if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert, no need to update insert coverage
	    						  correctlyMatedReads++;
	    						  correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
	    					  }else {
	    						  wronglyDistanceReads++;
	    						  wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
	    					  }
	    				  } else {
	    					  //pair is wrongly oriented
	    					  wronglyOrientedReads++;
	    					  wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
	    				  }
	    			  }
	    		  } else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
	    			  //Count inter-chrom pairs
	    			  matedDifferentContig++;
	    			  matedDifferentContigLength += bam_cigar2qlen(core,cigar);
	    		  } else if(core->flag&BAM_FMUNMAP) {
		    		  // if mate read is unmapped
	    			  singletonReads++;
	    			  singletonReadsLength =+ bam_cigar2qlen(core,cigar);
	    		  }


	    		  if (core->flag&BAM_FPROPER_PAIR) {
	    			  //Is part of a proper pair
	    			  matedReads ++; // increment number of mated reads
	    			  matedReadsLength += bam_cigar2qlen(core,cigar); // add the length of the read aligne as proper mate (not necessary correctly mated)
	    		  }

	    		  if (core->flag&BAM_FDUP) {   //This is a duplicate. Don't count it!.
	    			  ++duplicates;
	    		  }
	    	  } else {
	    		  ++zeroQualityReads;

	    	  }
	      }
    }

    cout << "LIBRARY STATISTICS\n";
    cout << "\ttotal reads number " << reads << "\n";
    cout << "\ttotal mapped reads " << mappedReads << "\n";
    cout << "\ttotal unmapped reads " << unmappedReads << "\n";
    cout << "\tproper pairs " << matedReads << "\n";
    cout << "\tzero quality reads " << zeroQualityReads << "\n";
    cout << "\tcorrectly oriented " << correctlyMatedReads << "\n";
    cout << "\twrongly oriented " << wronglyOrientedReads << "\n";
    cout << "\twrongly distance " << wronglyDistanceReads << "\n";
    cout << "\twrongly contig " <<  matedDifferentContig << "\n";
    cout << "\tsingletons " << singletonReads << "\n";

    uint32_t total = correctlyMatedReads + wronglyOrientedReads + wronglyDistanceReads + matedDifferentContig + singletonReads;
    cout << "\ttotal " << total << "\n";
    cout << "\tCoverage statistics\n";

    C_A = readsLength/(float)genomeLength;
     S_A = insertsLength/(float)genomeLength;
    C_M = correctlyMatedReadsLength/(float)genomeLength;
    C_W = (wronglyDistanceReadsLength + wronglyOrientedReadsLength)/(float)genomeLength;
     C_S = (singletonReadsLength)/(float)genomeLength;
    C_C = matedDifferentContigLength/(float)genomeLength;

    Qk = sqrt(Qk/counterK);
    insertStd = Qk;

    cout << "\tC_A = " << C_A << endl;
    cout << "\tS_A = " << S_A << endl;
    cout << "\tC_M = " << C_M << endl;
    cout << "\tC_W = " << C_W << endl;
    cout << "\tC_S = " << C_S << endl;
    cout << "\tC_C = " << C_C << endl;
    cout << "\tMean Insert length = " << Mk << endl;
    cout << "\tStd Insert length = " << Qk << endl;
    cout << "----------\n";



}


