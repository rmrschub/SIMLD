#include <stdexcept>
#include <iostream>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <string>
#include <serd/serd.h>
#include <set>
#include <vector>
#include <tuple>
#include <sys/stat.h>
#include "marisa.h"


namespace SIMLD
{
	static uint64_t getSize(std::istream &in);
	static uint64_t getSize(const char *file);

	enum RDFNotation 
	{
		XML, // No longer supported.
		NTRIPLES,
		TURTLE,
		N3, // Not supported.
		NQUAD, // Deprecated: use `NQUADS' instead.
		JSON, // Not supported.
		NQUADS,
		TRIG,
	};

	class ParseException: public std::exception 
	{
		protected:
			uint64_t byte;
			uint64_t line;
			uint32_t column;
			std::string reason;
		public:
			ParseException(uint64_t byte, uint64_t line, uint32_t column, std::string reason) :
				byte(byte),
				line(line),
				column(column),
				reason(reason) {  }
			ParseException(uint64_t line, uint32_t column, std::string reason) :
				// byte(byte),
				line(line),
				column(column),
				reason(reason) {  }
			ParseException(uint64_t line, std::string reason) :
				// byte(byte),
				line(line),
				// column(column),
				reason(reason) {  }
			ParseException(std::string reason) :
				// byte(byte),
				// line(line),
				// column(column),
				reason(reason) {  }

			/** Returns a C-style character string describing the general cause
			 *  of the current error.  */
			virtual const char* what() const throw() {
				return reason.c_str();
			}

			virtual ~ParseException() throw() {}
	};
	

    class NTriplesParser 
    {
        private:
	        SerdEnv *env;
			std::vector<std::tuple<std::string, std::string, std::string>> *container;
			marisa::Keyset *keyset;
	        std::string error;
	        uint64_t numByte;

	        std::string getString(const SerdNode *term);
	        std::string getStringObject(const SerdNode *term, const SerdNode *dataType, const SerdNode *lang);
	        SerdSyntax getParserType(SIMLD::RDFNotation notation);

        public:
            NTriplesParser();
	        virtual ~NTriplesParser();

//	        void doParse(const char *fileName, const char *baseUri, RDFNotation notation, bool ignoreErrors, RDFCallback *callback);
			void doParse(const char *fileName, const char *baseUri, SIMLD::RDFNotation notation, bool ignoreErrors, std::vector<std::tuple<std::string, std::string, std::string>>* container, marisa::Keyset * keyset);

	        friend SerdStatus hdtserd_on_statement(void               *handle,
	                                               SerdStatementFlags  flags,
	                                               const SerdNode     *graph,
	                                               const SerdNode     *subject,
	                                               const SerdNode     *predicate,
	                                               const SerdNode     *object,
	                                               const SerdNode     *datatype,
	                                               const SerdNode     *lang);

	        friend SerdStatus hdtserd_on_prefix(void           *handle,
                                                const SerdNode *name,
	                                            const SerdNode *uri);

	        friend SerdStatus hdtserd_on_base(void *handle, const SerdNode *uri);
	        friend SerdStatus hdtserd_on_error(void *handle, const SerdError *e);

    };

}   // namespace SIMLD