#include "ncl/othelpers.h"
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
using namespace std;
long gStrictLevel = 2;
////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void processFilepath(
	const char * filename, // name of the file to be read
	ostream *out, // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
	ProcessedTreeValidationFunction func, /*!< your pointer to your callback function */
	void * blob) {
	assert(filename);
	try {
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		if (gStrictLevel < 2)
			treesB->SetAllowImplicitNames(true);
		treesB->setValidationCallbacks(func, blob);
		if (gStrictLevel < 2) {
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);
		}
		cerr << "Executing" <<endl;
		try {
			nexusReader.ReadFilepath(filename, fmt);
		} catch(...) {
			nexusReader.DeleteBlocksFromFactories();
			throw;
		}
		nexusReader.DeleteBlocksFromFactories();
	} catch (const NxsException &x) {
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

int readFilepathAsNEXUS(const char *filename,
						MultiFormatReader::DataFormatType fmt,
						ProcessedTreeValidationFunction func, 
						void * blob) {
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		outStream = &cout;
		processFilepath(filename, outStream, fmt, func, blob);
	} catch (...) {
		cerr << "Parsing of " << filename << " failed (with an exception)" << endl;
		return 1;
	}
	return 0;
}


/*! \returns 0 on success*/
int readFilesListedInFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt,
						ProcessedTreeValidationFunction func, /*!< your pointer to your callback function */
	  				  	void * blob)
	{
	ifstream masterStream(masterFilepath, ios::binary);
	if (masterStream.bad())
		{
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
		}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good())
		{
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			{
			int rc = readFilepathAsNEXUS(filename, fmt, func, blob);
			if (rc != 0)
				return rc;
			}
		}
	return 0;
	}
