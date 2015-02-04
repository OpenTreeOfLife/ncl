#include <fstream>
#include "ncl/othelpers.h"
using namespace std;

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
		writeOttSet(std::cerr, "  ", gMRCADesignatorSet, "\n");
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
	gOTCLI.readFilepath(filepath, newTreeHook);
	return gOTCLI.exitCode;
}
