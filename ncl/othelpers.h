#if !defined OT_HELPERS_H
#define OT_HELPERS_H
#include <ostream>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <cassert>
#include "ncl/nxstreesblock.h"
#include "ncl/nxstreesblock.h"
#include "ncl/nxsstring.h"
#include "ncl/nxsmultiformat.h"
template<typename T>
const std::string & getStrOrThrow(const T & nd, const std::map<T, std::string> & tipNameMap);
long getOTTIndex(const NxsTaxaBlockAPI * taxa,
						const NxsSimpleNode & nd,
						std::map<const NxsSimpleNode *, long> *expandedNodesPtr=0);
const NxsSimpleNode * findNextSignificantNode(const NxsSimpleNode * node,
											  const std::map<const NxsSimpleNode *,std::set<long> > & ndp2mrca);
void describeUnnamedNode(const NxsSimpleNode *nd,
						 std::ostream & out,
						 unsigned int anc,
						 const std::map<const NxsSimpleNode *, std::string> & tipNameMap,
						 bool useNdNames);
void writeSet(std::ostream & out, const char *indent, const std::set<long> &fir, const char * sep);
void fillTipOTTIDs(const std::map<long, const NxsSimpleNode *> &taxonomy,
					long ottID,
					std::set<long> & tipOTTIDs,
					std::map<const NxsSimpleNode*, long> & taxNode2ottID);
bool isTheMrcaInThisLeafSet(const NxsSimpleNode * nd,
							const std::map<const NxsSimpleNode *, std::set<long> > & refNdp2mrcaThisLeafSet);
const std::string & getLeftmostDesName(const NxsSimpleNode *nd,
									   const std::map<const NxsSimpleNode *, std::string> & tipNameMap,
									   bool useNdNames);
const std::string &  getRightmostDesName(const NxsSimpleNode *nd,
										 const std::map<const NxsSimpleNode *, std::string> & tipNameMap,
										 bool useNdNames);
const NxsSimpleNode * findFirstBranchingAnc(const NxsSimpleNode *nd);
void writeNewickOTTIDs(std::ostream &out,
					   const NxsSimpleTree * tree,
					   const std::map<const NxsSimpleNode *, std::set<long> > & ndp2mrca);

int readFilepathAsNEXUS(const char *filename,
						MultiFormatReader::DataFormatType fmt,
						ProcessedTreeValidationFunction func, 
						void * blob=0L);

int readFilesListedInFile(const char *masterFilepath,
						  MultiFormatReader::DataFormatType fmt,
						  ProcessedTreeValidationFunction func, /*!< your pointer to your callback function */
	  				  	  void * blob=0L);
template<typename T>
inline const std::string & getStrOrThrow(const T & nd, const std::map<T, std::string> & tipNameMap) {
	typename std::map<T, std::string>::const_iterator tnIt = tipNameMap.find(nd);
	if (tnIt == tipNameMap.end()) {
		throw NxsException("AssertionError: Key not found");
	}
	return tnIt->second;
}

/// impl


inline void writeNewickSubtree(std::ostream & out,
								const NxsSimpleNode * sr,
								std::map<const NxsSimpleNode *, std::string > & leafNode2name) {
	assert(sr != 0);
	if (sr->IsTip()) {
		out << NxsString::GetEscaped(leafNode2name[sr]);
	} else {
		out << '(';
		bool first = true;
		const std::vector<NxsSimpleNode *> children = sr->GetChildren();
		for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			if (first) {
				first = false;
			} else {
				out << ',';
			}
			writeNewickSubtree(out, child, leafNode2name);
		}
		assert(!first);
		out << ')';
	}
}


void useChildrenToFillMRCASet(const NxsSimpleNode * nd,
				 std::map<const NxsSimpleNode *,std::set<long> > &n2m) {
	const std::vector<NxsSimpleNode *> children = nd->GetChildren();
	std::set<long> & mrca = n2m[nd];
	for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		assert(n2m.find(c) != n2m.end());
		std::set<long> & csl = n2m[c];
		mrca.insert(csl.begin(), csl.end());
	}
}


inline long ottIDFromName(const std::string & n) {
	//cout << "name \"" << n << "\"\n";
	if (n.empty()) {
		return -1;
	}
	const unsigned lastInd = n.length() - 1;
	unsigned currInd = lastInd;
	const char * c = n.c_str();
	if (strchr("0123456789", c[currInd]) == 0) {
		return -2;
	}
	while (currInd > 1) {
		--currInd;
		if (strchr("0123456789", c[currInd]) == 0) {
			++currInd;
			break;
		}
	}
	long conv = -2;
	NxsString::to_long(c + currInd, &conv);
	return conv;
}

inline const NxsSimpleNode * findNextSignificantNode(const NxsSimpleNode * node,
													 const std::map<const NxsSimpleNode *,std::set<long> > & ndp2mrca) {
	const NxsSimpleNode * currNode = node;
	for (;;) {
		std::map<const NxsSimpleNode *, std::set<long> >::const_iterator mIt = ndp2mrca.find(currNode);
		assert(mIt != ndp2mrca.end());
		const std::set<long> & oset = mIt->second;
		std::vector<NxsSimpleNode *> children = currNode->GetChildren();
		const NxsSimpleNode * sc = 0L;
		for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			const std::map<const NxsSimpleNode *, std::set<long> >::const_iterator nIt = ndp2mrca.find(child);
			if (nIt != ndp2mrca.end()) {
				if (sc == 0L) {
					sc = child;
				} else {
					return currNode; // more than one child is in ndp2mrca, so this node is significant
				}
			}
		}
		if (sc == 0L) {
			const char * msg = "Failing. Node found with ottIDs marked, but no children with ottIDs marked";
			std::cerr << msg << ":\n";
			for (std::set<long>::const_iterator oIt = oset.begin(); oIt != oset.end(); ++oIt) {
				if (oIt != oset.begin()) {
					std::cerr << ", ";
				}
				std::cerr << *oIt;
			}
			std::cerr << std::endl;
			assert(false);
			throw NxsException(msg);
		}
		const std::set<long> & dset = ndp2mrca.find(sc)->second;
		if (dset != oset) {
			const char * msg = "Failing. Internal node found with an ottID assignment. At this point the ottID should map to leaves.";
			std::cerr << msg << " Par ottIDs:\n";
			for (std::set<long>::const_iterator oIt = oset.begin(); oIt != oset.end(); ++oIt) {
				if (oIt != oset.begin()) {
					std::cerr << ", ";
				}
				std::cerr << *oIt;
			}
			std::cerr << std::endl;
			std::cerr << "child ottIDs:\n";
			for (std::set<long>::const_iterator oIt = dset.begin(); oIt != dset.end(); ++oIt) {
				if (oIt != dset.begin()) {
					std::cerr << ", ";
				}
				std::cerr << *oIt;
			}
			std::cerr << std::endl;
			assert(false);
			throw NxsException(msg);
		}
		currNode = sc;
	}
}

// suppresses nodes of out-degree=1
inline void writeSubtreeNewickOTTIDs(std::ostream &out,
									 const NxsSimpleNode * node,
									 const std::map<const NxsSimpleNode *, std::set<long> > & ndp2mrca) {
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator nIt = ndp2mrca.find(node);
	assert(nIt != ndp2mrca.end());
	const std::set<long> & ottIDSet = nIt->second;
	if (nIt->second.size() == 1) {
		const long ottID = *ottIDSet.begin();
		out << "ott" << ottID;
	} else {
		assert(nIt->second.size() > 1);
		const NxsSimpleNode * nsn = findNextSignificantNode(node, ndp2mrca);
		out << '(';
		unsigned numcwritten = 0;
		std::vector<NxsSimpleNode *> children = nsn->GetChildren();
		for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			nIt = ndp2mrca.find(child);
			if (nIt != ndp2mrca.end()) {
				if (numcwritten > 0) {
					out << ',';
				}
				writeSubtreeNewickOTTIDs(out, child, ndp2mrca);
				++numcwritten;
			}
		}
		assert(numcwritten > 1);
		out << ')';
	}
}

inline long getOTTIndex(const NxsTaxaBlockAPI * taxa,
						const NxsSimpleNode & nd,
						std::map<const NxsSimpleNode *, long> *expandedNodesPtr) {
	const NxsSimpleNode * ndp = &nd;
	if (expandedNodesPtr) {
		std::map<const NxsSimpleNode *, long>::const_iterator expIt = expandedNodesPtr->find(ndp);
		if (expIt != expandedNodesPtr->end()) {
			return expIt->second;
		}
	}
	const std::string & name = nd.GetName();
	if (name.empty()) {
		const unsigned ind = nd.GetTaxonIndex();
		if (ind < taxa->GetNumTaxonLabels()) {
			const std::string tn = taxa->GetTaxonLabel(ind);
			return ottIDFromName(tn);
		}
		return -1;
	}
	return ottIDFromName(name);
}


inline void describeUnnamedNode(const NxsSimpleNode *nd,
						 std::ostream & out, unsigned int anc,
						 const std::map<const NxsSimpleNode *, std::string> & tipNameMap,
						 bool useNdNames) {
	if (useNdNames && nd->GetName().length() > 0) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before \"" << nd->GetName() << "\"" << std::endl;
		} else {
			out << "the node \"" << nd->GetName() << "\"" << std::endl;
		}
		return;
	}
	std::vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 0) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before the leaf \"" << getStrOrThrow(nd, tipNameMap)  << "\"" << std::endl;
		} else {
			out << "the leaf \"" << getStrOrThrow(nd, tipNameMap)  << "\"" << std::endl;
		}
	} else if (outDegree == 1U) {
		describeUnnamedNode(children[0], out, anc + 1, tipNameMap, useNdNames);
	} else {
		std::string left = getLeftmostDesName(children[0], tipNameMap, useNdNames);
		std::string right = getRightmostDesName(children[outDegree - 1], tipNameMap, useNdNames);
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right <<'\"' << std::endl;
		} else {
			out <<  "MRCA of \"" << left << "\" and " << "\"" << right <<'\"' << std::endl;
		}
	}
}

const NxsSimpleNode * findMRCAFromIDSet(const std::map<long, const NxsSimpleNode *> & ref,
	                           			const std::set<long> & idSet,
	                           			long trigger) {
	std::map<const NxsSimpleNode *, unsigned int> n2c;
	long shortestPathLen = -1;
	const NxsSimpleNode * shortestPathNode = 0;
	for (std::set<long>::const_iterator toIt = idSet.begin(); toIt != idSet.end(); ++toIt) {
		const NxsSimpleNode * nd = 0;
		std::map<long, const NxsSimpleNode *>::const_iterator rIt = ref.find(*toIt);
		if (rIt == ref.end()) {
			std::cerr << "tip " << *toIt << " a descendant of " << trigger << " not found.\n";
			assert(false);
		}
		nd = rIt->second;
		long currPathLen = 0;
		while (nd != 0) {
			n2c[nd] += 1;
			currPathLen += 1;
			nd = nd->GetEdgeToParentRef().GetParent();
		}
		if (shortestPathLen < 0 || currPathLen < shortestPathLen) {
			shortestPathLen = currPathLen;
			shortestPathNode = rIt->second;
		}
	}
	const unsigned nTips = idSet.size();
	const NxsSimpleNode * cn = shortestPathNode;
	while (cn != 0) {
		if (n2c[cn] == nTips) {
			return cn;
		}
		cn = cn->GetEdgeToParentRef().GetParent();
	}
	assert(false);
	return 0L;
}

inline void writeOttSet(std::ostream & out,
						const char *indent,
						const std::set<long> &fir,
						const char * sep) {
	for (std::set<long>::const_iterator rIt = fir.begin(); rIt != fir.end(); ++rIt) {
		if (rIt != fir.begin()) {
			out << sep;
		}
		out << indent << "ott" << *rIt;
	}
}

inline void writeOttSetDiff(std::ostream & out,
							const char *indent,
							const std::set<long> &fir,
							const char *firN,
							const std::set<long> & sec,
							const char *secN) {
	for (std::set<long>::const_iterator rIt = fir.begin(); rIt != fir.end(); ++rIt) {
		if (sec.find(*rIt) == sec.end()) {
			out << indent << "ott" << *rIt << " is in " << firN << " but not " << secN << "\n";
		}
	}
	for (std::set<long>::const_iterator rIt = sec.begin(); rIt != sec.end(); ++rIt) {
		if (fir.find(*rIt) == fir.end()) {
			out << indent << "ott" << *rIt << " is in " << secN << " but not " << firN << "\n";
		}
	}
}

inline const std::string & getLeftmostDesName(const NxsSimpleNode *nd,
									   const std::map<const NxsSimpleNode *,
									   std::string> & tipNameMap, bool useNdNames) {
	if (useNdNames) {
		const std::string & name = nd->GetName();
		if (!name.empty()) {
			return name;
		}
	}
	if (nd->IsTip()) {
		return getStrOrThrow(nd, tipNameMap);
	}
	return getLeftmostDesName(nd->GetFirstChild(), tipNameMap, useNdNames);
}

inline const std::string &  getRightmostDesName(const NxsSimpleNode *nd,
									const std::map<const NxsSimpleNode *, std::string> & tipNameMap,
									bool useNdNames) {
	if (useNdNames) {
		const std::string & name = nd->GetName();
		if (!name.empty()) {
			return name;
		}
	}
	if (nd->IsTip()) {
		return getStrOrThrow(nd, tipNameMap);
	}
	return getRightmostDesName(nd->GetLastChild(), tipNameMap, useNdNames);
}

// find most recent anc of nd with out-degree > 1
inline const NxsSimpleNode * findFirstBranchingAnc(const NxsSimpleNode *nd) {
	const NxsSimpleNode * anc = nd->GetEdgeToParentRef().GetParent();
	assert(anc);
	while (anc->GetOutDegree() == 1) {
		anc = anc->GetEdgeToParentRef().GetParent();
	}
	return anc;
}


inline void writeNewickOTTIDs(std::ostream &out,
					   const NxsSimpleTree * tree,
					   const std::map<const NxsSimpleNode *, std::set<long> > & ndp2mrca) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	const NxsSimpleNode * root = nodes[0];
	writeSubtreeNewickOTTIDs(out, root, ndp2mrca);
	out << ";\n";
}


inline bool IsRedundantNodeAroundTip(const NxsSimpleNode * nd) {
	if (nd->IsTip()) {
		return true;
	}
	if (nd->GetOutDegree() == 1) {
		return IsRedundantNodeAroundTip(nd->GetFirstChild());
	}
	return false;
}

class OTCLI {
	public:
		OTCLI(const char *title, const char *usage)
			:exitCode(0),
			verbose(false),
			strictLevel(2),
			currReadingDotTxtFile(false),
			fmt(MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT),
			titleStr(title),
			usageStr(usage) {
			}
		int exitCode;
		bool verbose;
		long strictLevel;
		bool currReadingDotTxtFile;
		std::string currentFilename;
		std::string currTmpFilepath;

		int readFilepath(const std::string &fp,
						  ProcessedTreeValidationFunction func=0L,
						  void * blob=0L);
		bool parseArgs(int argc, char *argv[], std::vector<std::string> & args);
		void printHelp(std::ostream & out);
		bool isDotTxtFile(const std::string &fp);
	private:
		MultiFormatReader::DataFormatType fmt;
		std::string titleStr;
		std::string usageStr;
};


#endif