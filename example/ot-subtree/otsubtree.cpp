#include <fstream>
#include "ncl/othelpers.h"
using namespace std;

class OTCLI {
	public:
		OTCLI(const char *title, const char *usage)
			:exitCode(0),
			verbose(false),
			strictLevel(2),
			fmt(MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT),
			titleStr(title),
			usageStr(usage) {
			}
		int exitCode;
		bool verbose;
		long strictLevel;

		int readFilepath(const std::string &fp,
						  ProcessedTreeValidationFunction func=0L,
						  void * blob=0L) {
			return readFilepathAsNEXUS(fp.c_str(), this->fmt, func, blob);
		}
		bool parseArgs(int argc, char *argv[], std::vector<std::string> args);
		void printHelp(ostream & out);
	private:
		MultiFormatReader::DataFormatType fmt;
		std::string titleStr;
		std::string usageStr;
};

inline void OTCLI::printHelp(ostream & out) {
	out << this->titleStr << "\n";
	out << "\nThe most common usage is simply:\n";
	out << this->usageStr << "\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -v verbose output\n\n";
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -frelaxedphyliptree  newick(this is also the default)\n";
	out << "        The complete list of format names that can follow the -f flag is:\n";
	vector<string> fmtNames =  MultiFormatReader::getFormatNames();
	for (vector<string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n) {
		out << "            "<< *n << "\n";
	}
}

inline bool OTCLI::parseArgs(int argc, char *argv[], std::vector<std::string> args) {
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h') {
			this->printHelp(cout);
			this->exitCode = 1;
			return false;
		} else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v') {
			this->verbose = true;
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &(this->strictLevel)))) {
				cerr << "Expecting an integer after -s\n" << endl;
				this->printHelp(cerr);
				this->exitCode = 2;
				return false;
			}
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f') {
			fmt = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2) {
				string fmtName(filepath + 2, slen - 2);
				fmt =  MultiFormatReader::formatNameToCode(fmtName);
			}
			if (fmt == MultiFormatReader::UNSUPPORTED_FORMAT) {
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				this->exitCode = 2;
				return false;
			}
		} else {
			const string filepathstr(filepath);
			args.push_back(filepathstr);
		}
	}
	this->exitCode = 0;
	return true;
}

OTCLI gOTCLI("otsubtree: takes a newick tree file and a list of (numeric) ott IDs that are MRCA designators. Prints the subtree defined by the MRCA of the identifiers",
				"otsubtree mytree.tre 514 775241 > subtree.tre");



bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
set<long> gMRCADesignatorSet;


bool processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	map<const NxsSimpleNode *, set<long> > refNdp2mrca;
	map<const NxsSimpleNode *, string > leafNode2name;
	const unsigned numMRCADesignators = gMRCADesignatorSet.size();
	assert(numMRCADesignators > 1);
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		if (nd->IsTip()) {
			long ottID = getOTTIndex(tb, *nd);
			assert(ottID >= 0);
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const string tn = tb->GetTaxonLabel(ind);
			leafNode2name[nd] = tn;
			if (gMRCADesignatorSet.find(ottID) != gMRCADesignatorSet.end()) {
				gMRCADesignatorSet.erase(ottID);
				refNdp2mrca[nd].insert(ottID);
			}
		} else {
			useChildrenToFillMRCASet(nd, refNdp2mrca);
			if (refNdp2mrca[nd].size() == numMRCADesignators) {
				writeNewickSubtree(cout, nd, leafNode2name);
				cout << ";\n";
				return true;
			}
		}
	}
	if (gMRCADesignatorSet.empty()) {
		std::cerr << "Very odd: no node found that is an ancestor of all MRCA designators, but all designators found.\n";
	} else {
		std::cerr << "There following MRCA designator(s) not found (they all have to be leaf nodes):\n";
		for (set<long>::const_iterator mIt = gMRCADesignatorSet.begin(); mIt != gMRCADesignatorSet.end(); ++mIt) {
			std::cerr << *mIt << "\n";
		}
	}
	return false;
}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB) {
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gOTCLI.verbose) {
		cerr << "Read tree " <<  gTreeCount<< '\n';
	}
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	vector<string> parNames;
	NxsSimpleTree nst = NxsSimpleTree(ftd, 0.0, 0, true);
	if (!processRefTree(taxa, &nst)) {
		gOTCLI.exitCode = 1;
	}
	return false;
}


int main(int argc, char *argv[]) {
	std::vector<std::string> args;
	if (!gOTCLI.parseArgs(argc, argv, args)) {
		return 1;
	}
	if (args.size() < 3) {
		cerr << "Expecting a file name and at least 2 MRCA designators.\n";
	}
	std::string filepath = args[0];
	for (int j = 1; j < args.size(); ++j) {
		long mott;
		if (!NxsString::to_long(argv[j], &mott) || mott < 1) {
			cerr << "Expecting positive integer for an OTT ID as a MRCA designators.\n";
			return 3;
		}
		gMRCADesignatorSet.insert(mott);
	}
	gOTCLI.readFilepath(filepath, 0L, 0L);
	return gOTCLI.exitCode;
}
