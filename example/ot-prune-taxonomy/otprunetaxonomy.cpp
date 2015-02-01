//	Copyright (C) 2007-2008 Mark T. Holder
//
//	This file is part of NCL (Nexus Class Library).
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc.,
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
#include <fstream>
#include "ncl/othelpers.h"

using namespace std;
extern long gStrictLevel;
bool gVerbose = false;
bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
const NxsSimpleNode * findMRCAFromIDSet(const map<long, const NxsSimpleNode *> & ref, const set<long> & idSet, long trigger);

/* use some globals, because I'm being lazy... */
NxsSimpleTree * gTaxonTree = 0;
map<long, const NxsSimpleNode *> gOttID2RefNode;
map<const NxsSimpleNode *, string> gRefTipToName;
map<long, const NxsSimpleNode *> gOttID2TaxNode;
map<const NxsSimpleNode *, long> gTaxNode2ottID;
set<const NxsSimpleNode *> gSupportedNodes;
string gCurrentFilename;
string gCurrTmpFilepath;
ostream * gCurrTmpOstream = 0L;
bool gReadingTxtFile = false;
map<long, set<long> > gNonMono;
const bool gTrustNamedNodes = true;
map<const NxsSimpleNode *, long> gExpanded;
map<long, const NxsSimpleNode *> gTabooLeaf;
set<long> gTaxLeafOTTIDs;
map<const NxsSimpleNode *, set<long> > gAPrioriProblemNode;
bool gNoAprioriTests = true;
int gRefTreeNumNamedInternalsNodes = 0;
int gExitCode = 0;

void summarize(ostream & out) {
}


void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafOTTIDs.insert(ottID);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}

void processSourceTree(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {

}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gVerbose)
		cerr << "Read tree " <<  gTreeCount<< '\n';
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	vector<string> parNames;
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		processSourceTree(taxa, nst);
	}
	if (gTaxonTree != nst) {
		delete nst;
	}
	return false;
}


void printHelp(ostream & out) {
	out << "otprunetaxonomy take a newick form of OTT and a file listing the trees to be synthesized. It prints out a newick form of the taxonomy with any unused clades pruned.\n    NEXUStosplits <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -v verbose output\n\n";
	vector<string> fmtNames =  MultiFormatReader::getFormatNames();
	for (vector<string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n) {
		out << "            "<< *n << "\n";
	}
}

int main(int argc, char *argv[]) {
	MultiFormatReader::DataFormatType f(MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT);
	bool readfile = false;
	bool el = false;
	bool depth = false;
	bool brief = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v')
			gVerbose = true;
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel))) {
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
			}
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f') {
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2) {
				string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
			}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT) {
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				return 2;
			}
		} else {
			readfile = true;
			const string filepathstr(filepath);
			const size_t sp = filepathstr.find_last_of('/');
			if (sp == string::npos) {
				gCurrentFilename = filepathstr;
			} else {
				gCurrentFilename = filepathstr.substr(sp + 1);
			}
			gCurrTmpFilepath = string("tmp/") + gCurrentFilename;
			cerr << "gCurrTmpFilepath = " << gCurrTmpFilepath << '\n';
			std::ofstream tostream(gCurrTmpFilepath.c_str());
			gCurrTmpOstream = &tostream;
			const size_t fnl = gCurrentFilename.length();
			if (gCurrentFilename.substr(fnl - 4) == string(".txt")) {
				gReadingTxtFile = true;
			} else {
				gReadingTxtFile = false;
			}
			try {
				if (gReadingTxtFile) {
					gNoAprioriTests = false;
					cerr << "found text file. Treating each line as the name of a tree file.\n";
					readFilesListedInFile(filepath, f, newTreeHook, 0L);
				} else {
					readFilepathAsNEXUS(filepath, f, newTreeHook, 0L);
				}
				tostream.close();
			} catch (...) {
				tostream.close();
				cerr << "\nExiting due to and exception.\n";
				return 1;
			}
			gCurrTmpOstream = 0L;
		}
	}
	if (!readfile) {
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
	}
	summarize(cout);
	return gExitCode;
}

