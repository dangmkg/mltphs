module;

#include <ostream>
#include <istream>
#include <vector>
#include <map>
#include <string>

export module Json;

#define MAXSTRLEN 0x100000
#define MAXNUMLEN 32

export namespace Json {
	enum Type { string, number, object, array, boolean, null };
	enum Symbol { nosym, symString, symNumber, symTrue, symFalse, symNull, comma, colon, lbracket, rbracket, lbrace, rbrace };

	class Object {
	public:
		Type type;
		std::map<std::string, Object*> objfields;
		std::map<std::string, std::string> strfields;
		std::map<std::string, double> numfields;
		std::map<std::string, bool> boolfields;
		std::vector<Type> idx;

		~Object() {
			std::map<std::string, Object*>::iterator it;
			Object* obj;

			for (it = objfields.begin(); it != objfields.end(); it++) {
				obj = (*it).second;
				delete obj;
			}
		};
	};

	class Scope {
	public:
		Type type;
		int format, lev, nfields;
		std::string tabs;
		Scope* desc;
	};	

	export class Reader {
	private:
		std::istream* stream;
		Scope* topscope;
		Symbol sym; char ch; bool error;
		std::string str;
		double num;
		bool boolval;
		void getSym();
		void parseString();
		void parseNumber();
		void parseIdent();
		Object* parseObject();
		Object* parseArray();
	public:
		Reader();
		void setInput(std::istream* input);
		Object* parse();
	};

	export class Writer {
	private:
		std::ostream* stream;
		Scope* topscope;
		void writeKey(std::string const& name);
	public:
		Writer();
		void setOutput(std::ostream* output);

		void openObject(int format, std::string const& name);
		void closeObject();
		void openArray(int format, std::string const& name);
		void closeArray();
		void writeInt(std::string const& name, int value);
		void writeDouble(std::string const& name, double value);
		void writeArrayOfDouble(std::string const& name, double const arr[], size_t len);
	};

	void Json::Writer::writeKey(std::string const& name) {
		if (topscope != nullptr) {
			if (topscope->nfields > 0) *stream << ',';
			if (topscope->format == 1) *stream << "\n" << topscope->tabs;
			else *stream << ' ';
			topscope->nfields++;
		}
		if (name.compare("") != 0) *stream << '"' << name << "\": ";
	}

	void Writer::openObject(int format, std::string const& name) {
		int i;
		Scope* scope = new Scope();

		scope->type = object; scope->nfields = 0;
		scope->format = format; scope->tabs = "\t";
		if (topscope == nullptr) {
			scope->lev = 0; scope->desc = nullptr;
		}
		else {
			scope->lev = topscope->lev + 1; scope->desc = topscope;
			for (i = scope->lev; i > 0; i--) scope->tabs.append("\t");
		}
		writeKey(name); *stream << "{";
		topscope = scope;
	}

	void Writer::closeObject() {
		Scope* scope;
		if (topscope != nullptr) {
			scope = topscope; topscope = topscope->desc;
			if (scope->format == 1) {
				*stream << "\n";
				if (topscope != nullptr) *stream << topscope->tabs;
			}
			*stream << "}";
			delete scope;
		}
	}

	void Writer::openArray(int format, std::string const& name) {
		int i;
		Scope* scope = new Scope();

		scope->type = array; scope->nfields = 0;
		scope->format = format; scope->tabs = "\t";
		if (topscope == nullptr) {
			scope->lev = 0; scope->desc = nullptr;
		}
		else {
			scope->lev = topscope->lev + 1; scope->desc = topscope;
			for (i = scope->lev; i > 0; i--) scope->tabs.append("\t");
		}
		writeKey(name); *stream << "[";
		topscope = scope;
	}

	void Writer::closeArray() {
		Scope* scope;
		if (topscope != nullptr) {
			scope = topscope; topscope = topscope->desc;
			if (scope->format == 1) {
				*stream << "\n";
				if (topscope != nullptr) *stream << topscope->tabs;
			}
			*stream << "]";
			delete scope;
		}
	}

	void Json::Writer::writeInt(std::string const& name, int value) {
		if (topscope != nullptr) {
			writeKey(name); *stream << value;
		}
	}

	void Json::Writer::writeDouble(std::string const& name, double value) {
		if (topscope != nullptr) {
			writeKey(name); *stream << value;
		}
	}

	void Json::Writer::writeArrayOfDouble(std::string const& name, double const arr[], size_t const len) {
		size_t i;
		if (topscope != nullptr) {
			writeKey(name); *stream << "[";
			for (i = 0; i < len - 1; i++) *stream << arr[i] << ", ";
			if (len > 0) *stream << arr[len - 1];
			*stream << "]";
		}
	}

	void Json::Writer::setOutput(std::ostream* output) {
		stream = output;
	}

	Json::Writer::Writer() {
		stream = nullptr;
		topscope = nullptr;
	}
	
	void Reader::getSym() {
		while (!stream->eof() && ch <= ' ') stream->get(ch);
		if (ch < '0') {
			if (ch == 0x22) this->parseString();
			else if (ch == '-') this->parseNumber();
			else if (ch == ',') {
				sym = comma; stream->get(ch);
			}
			else {
				sym = nosym; stream->get(ch);
			}
		}
		else if (ch <= '9') this->parseNumber();
		else if (ch < 'a') {
			if (ch == ':') sym = colon;
			else if (ch == '[') sym = lbracket;
			else if (ch == ']') sym = rbracket;
			else sym = nosym;
			stream->get(ch);
		}
		else if (ch <= 'z') this->parseIdent();
		else {
			if (ch == '{') sym = lbrace;
			else if (ch == '}') sym = rbrace;
			else sym = nosym;
			stream->get(ch);
		}
	}

	void Reader::parseString() {
		int i = 0;
		bool escape = false;

		sym = symString; str.clear(); stream->get(ch);
		while (i < MAXSTRLEN && !stream->eof() && ch != '"') {
			if (!escape) {
				if (ch != '\\') {
					str.push_back(ch); i++;
				}
				else escape = true;
				stream->get(ch);
			}
			else {
				if (ch == '"' || ch == '\\' || ch == '/') str.push_back(ch);
				else if (ch == 'b') str.push_back('\b');
				else if (ch == 'f') str.push_back('\f');
				else if (ch == 'n') str.push_back('\n');
				else if (ch == 'r') str.push_back('\r');
				else if (ch == 't') str.push_back('\t');
				else if (ch == 'u') {
					error = true;
				}
				else error = true;
				escape = false;
				i++; stream->get(ch);
			}
		}
		if (stream->eof() || ch != '"') error = true;
		stream->get(ch);
	}

	void Reader::parseNumber() {
		int i = 0;
		std::string numstr;

		sym = symNumber;
		if (ch == '-') {
			numstr.push_back(ch);
			i++; stream->get(ch);
		}
		while (!stream->eof() && ch >= '0' && ch <= '9') {
			error = (i >= MAXNUMLEN);
			if (!error) numstr.push_back(ch);
			i++; stream->get(ch);
		}
		if (ch == '.') {
			error = (i >= MAXNUMLEN);
			if (!error) numstr.push_back(ch);
			i++; stream->get(ch);
			if (stream->eof() || ch < '0' || ch > '9') error = true;
			while (!stream->eof() && ch >= '0' && ch <= '9') {
				error = (i >= MAXNUMLEN);
				if (!error) numstr.push_back(ch);
				i++; stream->get(ch);
			}
		}
		if (ch == 'E' || ch == 'e') {
			error = (i >= MAXNUMLEN);
			if (!error) numstr.push_back(ch);
			i++; stream->get(ch);
			if (ch == '-' || ch == '+') {
				error = (i >= MAXNUMLEN);
				if (!error) numstr.push_back(ch);
				i++; stream->get(ch);
			}
			if (stream->eof() || ch < '0' || ch > '9') error = true;
			while (!stream->eof() && ch >= '0' && ch <= '9') {
				error = (i >= MAXNUMLEN);
				if (!error) numstr.push_back(ch);
				i++; stream->get(ch);
			}
		}
		if (!error) num = std::stod(numstr);
	}

	void Reader::parseIdent() {
		int i = 0;
		std::string ident;

		while (!stream->eof() && i < 10) {
			ident.push_back(ch); stream->get(ch); i++;
		}
		if (ident.compare("true") == 0) sym = symTrue;
		else if (ident.compare("false") == 0) sym = symFalse;
		else if (ident.compare("null") == 0) sym = symNull;
		else sym = nosym;
	}

	Object* Reader::parseObject() {
		std::string key;
		Object* obj = new Object();

		obj->type = Type::object;
		getSym();
		while (sym == symString) {
			key = this->str; getSym();
			if (sym == colon) getSym(); else error = true;
			if (sym == symNumber) {
				obj->numfields[key] = this->num;
				getSym();
			}
			else if (sym == symString) {
				obj->strfields[key] = this->str;
				getSym();
			}
			else if (sym == lbrace) {
				obj->objfields[key] = parseObject();
			}
			else if (sym == lbracket) {
				obj->objfields[key] = parseArray();
				getSym();
			}
			else if (sym == symTrue) {
				obj->boolfields[key] = true;
				getSym();
			}
			else if (sym == symFalse) {
				obj->boolfields[key] = false;
				getSym();
			}
			else if (sym == symNull) {
				obj->objfields[key] = nullptr;
				getSym();
			}
			else error = true;
			while (sym == comma) getSym();
		}
		if (sym == rbrace) getSym(); else error = true;
		return obj;
	}

	Object* Reader::parseArray() {
		std::string key;
		Object* obj = new Object();
		int i;

		obj->type = Type::array; i = 0;
		do {
			getSym();
			if (sym == symNumber) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::number);
				obj->numfields[key] = this->num;
				getSym();
			}
			else if (sym == symString) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::string);
				obj->strfields[key] = this->str;
				getSym();
			}
			else if (sym == lbrace) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::object);
				obj->objfields[key] = parseObject();
			}
			else if (sym == lbracket) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::array);
				obj->objfields[key] = parseArray();
			}
			else if (sym == symTrue) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::boolean);
				obj->boolfields[key] = true;
				getSym();
			}
			else if (sym == symFalse) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::boolean);
				obj->boolfields[key] = false;
				getSym();
			}
			else if (sym == symNull) {
				key = std::to_string(i); i++;
				obj->idx.push_back(Type::null);
				obj->objfields[key] = nullptr;
				getSym();
			}
			else if (sym == comma) {
				error = true;
			}
		} while (sym == comma);
		if (sym == rbracket) getSym(); else error = true;
		return obj;
	}

	Reader::Reader() {
		stream = nullptr;
		topscope = nullptr;
		ch = 0;
	}

	void Reader::setInput(std::istream* input) {
		stream = input; ch = 0;
	}

	Object* Reader::parse() {
		Object* obj = nullptr;

		error = false; sym = nosym;
		while (sym != lbrace && !stream->eof()) getSym();
		if (sym == lbrace && !stream->eof()) {
			obj = parseObject();
		}
		if (error) {
			delete obj; obj = nullptr;
		}

		return obj;
	}
}