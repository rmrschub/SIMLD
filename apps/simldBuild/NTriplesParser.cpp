#include "SIMLD/NTriplesParser.hpp"

namespace SIMLD
{

	uint64_t getSize(std::istream &in) 
	{
		long long begin = in.tellg();
		in.seekg(0, std::ios::end);
		long long end = in.tellg();
		in.seekg(begin, std::ios::beg);

		return end-begin;
	}

	uint64_t getSize(const char *file) 
	{
		#ifdef WIN32
		ifstream in(file);
		in.seekg(0, std::ios::end);
		uint64_t size = in.tellg();
		in.close();
		return size;
		#else
			struct stat fileStat;
			if(stat(file, &fileStat)==0) {
				return fileStat.st_size;
			}
			return 0;
		#endif
	}



   std::string NTriplesParser::getString(const SerdNode *term) 
    {
	    std::string out;
	    out.reserve(term->n_bytes + 2);

	    if(term->type==SERD_URI) 
        {
		    out.append((const char *)term->buf, term->n_bytes);
	    } 
        else if(term->type==SERD_BLANK) 
        {
		    out.append("_:");
		    out.append((const char *)term->buf, term->n_bytes);
	    } else if(term->type==SERD_CURIE) 
        {
		    SerdChunk uri_prefix, uri_suffix;
		    if (serd_env_expand(env, term, &uri_prefix, &uri_suffix)) 
            {
			    // ERROR BAD Curie / Prefix
		    }
		    out.append((const char *)uri_prefix.buf, uri_prefix.len);
		    out.append((const char *)uri_suffix.buf, uri_suffix.len);
	    }

	    return out;
    }

    std::string NTriplesParser::getStringObject(const SerdNode *term, const SerdNode *dataType, const SerdNode *lang) 
    {
	    if(term->type!=SERD_LITERAL) 
        {
		    return getString(term);
	    }

	    std::string out;
	    out.reserve(term->n_bytes + 2 +
	                (dataType ? dataType->n_bytes + 4 : 0) +
	                (lang     ? lang->n_bytes + 1     : 0));

	    out.push_back('\"');
	    out.append((const char *)term->buf, term->n_bytes);
	    out.push_back('\"');
	    if(lang!=NULL){
		    out.push_back('@');
		    out.append((const char *)lang->buf, lang->n_bytes);
	    }
	    if(dataType!=NULL) {
		    out.append("^^<");
		    out.append(getString(dataType));
		    out.push_back('>');
        }

	    return out;
    }

    SerdStatus hdtserd_on_error(void *handle, const SerdError *error) 
    {
	    fprintf(stderr, "error: %s:%u:%u: ",
	            error->filename, error->line, error->col);
	    vfprintf(stderr, error->fmt, *error->args);
	    throw std::runtime_error("Error parsing input.");
	    return error->status;
    }

    // Callback for base URI changes (@base directives)
    SerdStatus hdtserd_on_base(void *handle, const SerdNode *uri) 
    {
	    NTriplesParser *serdParser = reinterpret_cast<NTriplesParser *>(handle);

	    return serd_env_set_base_uri(serdParser->env, uri);
    }

    // Callback for namespace definitions (@prefix directives)
    SerdStatus hdtserd_on_prefix(void *handle, const SerdNode *name, const SerdNode *uri) 
    {
	    NTriplesParser *serdParser = reinterpret_cast<NTriplesParser *>(handle);

	    return serd_env_set_prefix(serdParser->env, name, uri);
    }

    // Callback for statements
    SerdStatus hdtserd_on_statement(void               *handle,
                                    SerdStatementFlags  flags,
                                    const SerdNode     *graph,
                                    const SerdNode     *subject,
                                    const SerdNode     *predicate,
                                    const SerdNode     *object,
                                    const SerdNode     *datatype,
                                    const SerdNode     *lang) 
    {
	    NTriplesParser *serdParser = reinterpret_cast<NTriplesParser *>(handle);

		serdParser->container->emplace_back(
			std::make_tuple(
				serdParser->getString(subject),
				serdParser->getString(predicate),
				serdParser->getStringObject(object, datatype, lang))
		);

		serdParser->keyset->push_back(serdParser->getString(subject), 1.0);
		serdParser->keyset->push_back(serdParser->getString(predicate), 1.0);
		serdParser->keyset->push_back(serdParser->getStringObject(object, datatype, lang), 1.0);

	    return SERD_SUCCESS;
    }

	NTriplesParser::NTriplesParser() : numByte(0) { }

	NTriplesParser::~NTriplesParser() { }

	SerdSyntax NTriplesParser::getParserType(SIMLD::RDFNotation notation) 
	{
		switch(notation)
		{
			case NQUAD: // Deprecated: use `NQUADS' instead.
				return SERD_NQUADS;
			case NQUADS:
				return SERD_NQUADS;
			case NTRIPLES:
				return SERD_NTRIPLES;
			case TRIG:
				return SERD_TRIG;
			case TURTLE:
				return SERD_TURTLE;
			default:
				throw ParseException("Serd parser only supports N-Triples, N-Quads, TriG, and Turtle.");
		}
	}

	void NTriplesParser::doParse(const char *fileName, const char *baseUri, SIMLD::RDFNotation notation, bool ignoreErrors, std::vector<std::tuple<std::string, std::string, std::string>>* container, marisa::Keyset *keyset)
	{
		this->container = container;
		this->keyset = keyset;
		this->numByte = getSize(fileName);

		// Create Base URI and environment
		SerdURI  base_uri = SERD_URI_NULL;
		SerdNode base = serd_node_new_file_uri((const uint8_t *)fileName, NULL, &base_uri, false);
		env = serd_env_new(&base);

		SerdReader* reader = serd_reader_new(
					getParserType(notation), this, NULL,
					(SerdBaseSink)hdtserd_on_base,
					(SerdPrefixSink)hdtserd_on_prefix,
					(SerdStatementSink)hdtserd_on_statement,
					NULL);

		serd_reader_set_error_sink(reader, hdtserd_on_error, NULL);

		const uint8_t* input=serd_uri_to_path((const uint8_t *)fileName);

		FILE *in_fd = fopen((const char*)input, "r");
		// TODO: fadvise sequential
		if (in_fd==NULL) 
		{
			throw ParseException("Could not open input file for parsing");
		}

		serd_reader_read_file_handle(reader, in_fd, (const uint8_t *)fileName);
		
		fclose(in_fd);

		serd_reader_free(reader);

		serd_env_free(env);
		
		serd_node_free(&base);
	}

}   // namespace SIMLD