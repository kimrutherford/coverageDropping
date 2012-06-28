/*
 * Contig.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#include "Contig.h"


Contig::Contig() {

}

Contig::Contig(unsigned int contigLength, unsigned int libraries) {
	this->contigLength = contigLength;
	this->libraries = libraries;
	this->readCoverage = new unsigned int*[libraries+1];
	this->spanCoverage = new unsigned int*[libraries+1];
	this->maxInserts   = new unsigned int[libraries];
	this->minInserts   = new unsigned int[libraries];
	this->meanInserts  = new unsigned int[libraries];

	this->edgeCutoff = new unsigned int[libraries + 1];
	this->dropCutoff = new unsigned int[libraries +1 ];

	this->readCov 	   = new float[libraries + 1];
	this->spanCov 	   = new float[libraries + 1];

	this->totalReadLength 	   = new unsigned long int[libraries + 1];
	this->totalSpanLength 	   = new unsigned long int[libraries + 1];

	for(unsigned int i=0; i < libraries; i ++) {
		this->readCoverage[i] = new unsigned int[contigLength];
		this->spanCoverage[i] = new unsigned int[contigLength];
		this->minInserts[i] = 0;
		this->maxInserts[i] = 1000000;
		this->meanInserts[i] = 10000;
	}
	this->readCoverage[libraries] = new unsigned int[contigLength];
	this->spanCoverage[libraries] = new unsigned int[contigLength];


	// reset memory location WRITE 0 everywhere
	for(unsigned int i = 0; i < libraries + 1 ; i ++) {
		for(unsigned int j = 0; j < contigLength ; j++) {
			this->readCoverage[i][j] = 0;
			this->spanCoverage[i][j] = 0;
		}
		this->edgeCutoff[i]  = 0;
		this->dropCutoff[i]  = 0;
		this->readCov[i]	 = 0;
		this->spanCov[i]	 = 0;
		this->totalReadLength[i] = 0;
		this->totalSpanLength[i] = 0;
	}



}

Contig::~Contig() {
	if(readCoverage != NULL) {
		for(unsigned int i=0 ; i < libraries + 1 ; i++) {
			delete [] this->readCoverage[i];
			delete [] this->spanCoverage[i];
		}
		delete [] this->readCoverage;
		delete [] this->spanCoverage;
	}
	delete [] this->maxInserts ;
	delete [] this->minInserts;
	delete [] this->meanInserts;

	delete [] this->edgeCutoff;
	delete [] this->dropCutoff;

	delete [] this->readCov;
	delete [] this->spanCov;

	delete [] this->totalReadLength;
	delete [] this->totalSpanLength;
}


void Contig::setLibraryLimits(unsigned int library, unsigned int minInsert, unsigned int maxInsert, unsigned int meanInsert, unsigned int edgeCutoff,  unsigned int dropCutoff) {
	if(library >= 0 and library < this->libraries) {
		this->minInserts[library]	= minInsert;
		this->maxInserts[library]	= maxInsert;
		this->meanInserts[library]	= meanInsert;
		this->edgeCutoff[library]	= edgeCutoff;
		this->dropCutoff[library]	= dropCutoff;

	} else {
		cout << "in Contig::setLibrariesLimits passed a library identifier larger than the number of libraries\n";
		exit(1);
	}
}


void Contig::updateReadCov(unsigned int start, unsigned int end, unsigned int library) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
		end = this->contigLength;
	}
	// now update
	for(unsigned int i = start; i< end ; i ++ ) {
		this->readCoverage[library][i]++;
		this->totalReadLength[library]++;

		this->readCoverage[this->libraries][i]++;
		this->totalReadLength[this->libraries]++;
	}

}

void Contig::updateSpanCov(unsigned int start, unsigned int end, unsigned int library) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
		end = this->contigLength;
	}
	// now update
	for(unsigned int i = start; i< end ; i ++ ) {
		this->spanCoverage[library][i]++;
		this->totalSpanLength[library]++;

		this->spanCoverage[this->libraries][i]++;
		this->totalSpanLength[this->libraries]++;

	}

}



void Contig::computeContigStats() {


	for(unsigned int i = 0; i < this->libraries + 1; i++) {
		this->readCov[i] = ((float)this->totalReadLength[i]/this->contigLength);
	}

	for(unsigned int i = 0; i < this->libraries + 1; i++) {
		this->spanCov[i] = ((float)this->totalSpanLength[i]/this->contigLength);
	}
}


void  Contig::computeCoverageDrops() {

	float covCutOff = 0;
	unsigned int dropsLength = 0;


	unsigned int drops = 0;
	for(unsigned int lib = 0; lib < this->libraries + 1; lib++) {
		bool dropFound = false;
		unsigned int dropStart = 0;
		unsigned int dropEnd   = 0;

		for(unsigned int pos = 0; pos < this->contigLength; pos++) {
			if(this->spanCoverage[lib][pos] <= covCutOff) {
				if(dropFound) {
					dropEnd++;
					dropsLength++;
				} else {
					dropStart = dropEnd = pos;
					dropFound = true;
					dropsLength++;
				}
			} else {
				if(dropFound) { // I am closing the drop
					//cout << "\t" << dropStart << "\t" << dropEnd << "\n";
					dropFound = false;
					drops++;
				}
			}
		}
		if(dropFound) { // I am closing the drop
			//cout << "\t" << dropStart << "\t" << dropEnd << "\n";
			dropFound = false;
			drops++;
		}
		float meanDropLength = (float)dropsLength/drops;
		cout << "\t" << lib << "\t" << drops << "\t" << meanDropLength << "\n";
	}
	drops = 0;
	dropsLength = 0;

}



void  Contig::computeCoverageDrops(string ctgID, unsigned int lib, ofstream &output) {

	unsigned int drops = 0;
	unsigned int dropsLength = 0;
	float meanDropLength;

	bool dropFound = false;
	unsigned int dropStart = 0;
	unsigned int dropEnd   = 0;

	unsigned int startCov = this->edgeCutoff[lib];
	unsigned int endCov	  = this->contigLength - this->edgeCutoff[lib];

	if(this->contigLength > this->meanInserts[lib] and endCov > startCov) { // contig MUST be at least long as the expected insert size otherwise "succhia"
		for(unsigned int pos = startCov; pos < endCov; pos++) {
			if(this->spanCoverage[lib][pos] <= this->dropCutoff[lib]) {
				if(dropFound) {
					dropEnd++;
					dropsLength++;
				} else {
					dropStart = dropEnd = pos;
					dropsLength++;
					dropFound = true;
				}
			} else {
				if(dropFound) { // I am closing the drop
					//cout << "\t" << dropStart << "\t" << dropEnd << "\n";
					dropFound = false;
					drops++;
				}
			}
		}
		if(dropFound) { // I am closing the drop
			//cout << "\t" << dropStart << "\t" << dropEnd << "\n";
			dropFound = false;
			drops++;
		}
		// store number of drops
		if(drops > 0) {
			meanDropLength = (float)dropsLength/drops;
			output << ctgID  << "\t" << this->contigLength  << "\t" << drops << "\t" << meanDropLength << "\n";
		} else {
			output << ctgID << "\t" << this->contigLength << "\t" << 0 << "\t" << 0 << "\n";
		}
	} else {
		output << ctgID  << "\t" << this->contigLength  << "\t" << -1 << "\t" << 0 << "\n";
	}




}



void Contig::setUpTotal() {
	unsigned int min = 99999999;
	for(unsigned int i = 0; i < this->libraries; i++) {
		if(this->edgeCutoff[i] < min) {
			min = this->edgeCutoff[i];
		}
	}
	this->edgeCutoff[this->libraries] = min;

	min = 99999999;
	for(unsigned int i = 0; i < this->libraries; i++) {
		if(this->dropCutoff[i] < min) {
			min = this->dropCutoff[i];
		}
	}
	this->dropCutoff[this->libraries] = min;




}


void Contig::updateContig(const bam1_t* b, unsigned int library) {
	const bam1_core_t* core =  &b->core;
	//uint32_t* cigar = bam1_cigar(b);
	int32_t alignmentLength = 0;
	int32_t startRead=0;
	int32_t endRead=0;
	int32_t startPaired=0;
	int32_t startInsert=0;
	int32_t endInsert=0;
	uint32_t iSize=0;

	uint32_t minInsert = this->minInserts[library];
	uint32_t maxInsert = this->maxInserts[library];

	if(!(core->flag&BAM_FUNMAP) && !(core->flag&BAM_FDUP) && !(core->flag&BAM_FSECONDARY) && !(core->flag&BAM_FQCFAIL)) { // if read has been mapped and it is not a DUPLICATE or a SECONDARY alignment
		alignmentLength = core->l_qseq; // bam_cigar2qlen(core,cigar);
		startRead = core->pos; // start position on the contig
		endRead = startRead + alignmentLength ; // position where reads ends

		updateReadCov(startRead, endRead, library); //  update read coverage

		iSize = abs(core->isize);

		if ((core->flag&BAM_FREAD1) //First in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			startPaired = core->mpos;
			if(startRead < startPaired) {
				iSize = (startPaired + core->l_qseq -1) - startRead; // insert size, I consider both reads of the same length
				startInsert = startRead;
				endInsert = startRead + iSize;

				// DO NOT CHECK EITHER ORIENTATION NOR DISTANCE
				//updateSpanCov(startInsert,endInsert, library);

				// DO NOT CHECK ORIENTATION, CHECK DISTANCE
				//if (minInsert <= iSize && iSize <= maxInsert) {
				//	updateSpanCov(startInsert,endInsert, library); // update spanning coverage
				//}

				// CHECK ORIENTATION, CHECK DISTANCE
				if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
						updateSpanCov(startInsert, endInsert, library); // update good read coverage
					}
				}
			} else {
				iSize = (startRead + alignmentLength - 1) - startPaired;
				startInsert = startPaired;
				endInsert = startInsert + iSize;

				// DO NOT CHECK EITHER ORIENTATION NOR DISTANCE
				//updateSpanCov(startInsert,endInsert, library); // update spanning coverage

				// DO NOT CHECK ORIENTATION, CHECK DISTANCE
				//if (minInsert <= iSize && iSize <= maxInsert) { //
				//	updateSpanCov(startInsert,endInsert, library); // update spanning coverage
				//}

				// CHECK ORIENTATION, CHECK DISTANCE
				if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
						updateSpanCov(startInsert,endInsert, library); // update good read coverage
					}
				}
			}
		} // other case are not considered here

	}
}



void Contig::print() {
	for(unsigned int i = 0; i < this->contigLength; i++) {
		cout << "(";
		for(unsigned int lib = 0; lib < this->libraries; lib++) {
			cout << this->spanCoverage[lib][i] << ",";
		}
		cout << this->spanCoverage[this->libraries][i] << ") ";
	}
	cout << "\n";
}

void Contig::printStats() {
	cout << "read coverage stats: (";
	for(unsigned int lib = 0; lib < this->libraries; lib++) {
		cout << this->readCov[lib] << ",";
	}
	cout << this->readCov[this->libraries] << ") ";

	cout << " (";
	for(unsigned int lib = 0; lib < this->libraries; lib++) {
		cout << this->spanCov[lib] << ",";
	}
	cout << this->spanCov[this->libraries] << ") ";

	cout << "\n";

}

void Contig::printLibrariesThresholds() {

	cout << "minLibraryInsert ";
	for(unsigned int lib = 0; lib < this->libraries ; lib++) {
		cout << this->minInserts[lib] << " ";
	}
	cout << "\n";


	cout << "meanLibraryInsert ";
	for(unsigned int lib = 0; lib < this->libraries ; lib++) {
		cout << this->meanInserts[lib] << " ";
	}
	cout << "\n";

	cout << "maxLibraryInsert ";
	for(unsigned int lib = 0; lib < this->libraries; lib++) {
		cout << this->maxInserts[lib] << " ";
	}
	cout << "\n";

	cout << "edgeCutOff ";
	for(unsigned int lib = 0; lib < this->libraries +1 ; lib++) {
		cout << this->edgeCutoff[lib] << " ";
	}
	cout << "\n";

	cout << "dropCutoff ";
	for(unsigned int lib = 0; lib < this->libraries +1 ; lib++) {
		cout << this->dropCutoff[lib] << " ";
	}
	cout << "\n";

}

