/**
 * Copyright (C) 2010 Johannes Weißl <jargon@molb.org>
 * License: your favourite BSD-style license
 *
 * See OptionParser.h for help.
 */

#include "OptionParser.h"

#include "utils/mardyn_assert.h"

#include <cstdlib>
#include <algorithm>
#include <complex>
#include <sstream>

#ifdef ENABLE_MPI
#include <mpi.h>  // For MPI_Finalize()
#endif

#if defined(ENABLE_NLS) && ENABLE_NLS
# include <libintl.h>
# define _(s) gettext(s)
#else
# define _(s) ((const char *) (s))
#endif



namespace optparse {

////////// auxiliary (string) functions { //////////
struct str_wrap {
	str_wrap(const std::string& l, const std::string& r) : lwrap(l), rwrap(r) {}
	str_wrap(const std::string& w) : lwrap(w), rwrap(w) {}
	std::string operator() (const std::string& s) { return lwrap + s + rwrap; }
	const std::string lwrap, rwrap;
};
template<typename InputIterator, typename UnaryOperator>
static std::string str_join_trans(const std::string& sep, InputIterator begin, InputIterator end, UnaryOperator op) {
	std::string buf;
	for (InputIterator it = begin; it != end; ++it) {
		if (it != begin)
			buf += sep;
		buf += op(*it);
	}
	return buf;
}
template<class InputIterator>
static std::string str_join(const std::string& sep, InputIterator begin, InputIterator end) {
	return str_join_trans<InputIterator, str_wrap>(sep, begin, end, str_wrap(""));
}
static std::string& str_replace(std::string& s, const std::string& patt, const std::string& repl) {
	size_t pos = 0, n = patt.length();
	while (true) {
		pos = s.find(patt, pos);
		if (pos == std::string::npos)
			break;
		s.replace(pos, n, repl);
		pos += repl.size();
	}
	return s;
}
static std::string str_replace(const std::string& s, const std::string& patt, const std::string& repl) {
	std::string tmp = s;
	str_replace(tmp, patt, repl);
	return tmp;
}
static std::string str_format(const std::string& s, size_t pre, size_t len, bool indent_first = true) {
	std::stringstream ss;
	std::string p;
	if (indent_first)
		p = std::string(pre, ' ');

	size_t pos = 0, linestart = 0;
	size_t line = 0;
	while (true) {
		bool wrap = false;

		size_t new_pos = s.find_first_of(" \n\t", pos);
		if (new_pos == std::string::npos)
			break;
		if (s[new_pos] == '\n') {
			pos = new_pos + 1;
			wrap = true;
		}
		if (line == 1)
			p = std::string(pre, ' ');
		if (wrap || new_pos + pre > linestart + len) {
			ss << p << s.substr(linestart, pos - linestart - 1) << std::endl;
			linestart = pos;
			line++;
		}
		pos = new_pos + 1;
	}
	ss << p << s.substr(linestart) << std::endl;
	return ss.str();
}
static std::string str_inc(const std::string& s) {
	std::stringstream ss;
	std::string v = (s != "") ? s : "0";
	long i;
	std::istringstream(v) >> i;
	ss << i+1;
	return ss.str();
}
static unsigned int cols() {
	unsigned int n = 80;
	const char *s = getenv("COLUMNS");
	if (s)
		std::istringstream(s) >> n;
	return n;
}
static std::string basename(const std::string& s) {
	std::string b = s;
	size_t i = b.find_last_not_of('/');
	if (i == std::string::npos) {
		if (b[0] == '/')
			b.erase(1);
		return b;
	}
	b.erase(i+1, b.length()-i-1);
	i = b.find_last_of("/");
	if (i != std::string::npos)
		b.erase(0, i+1);
	return b;
}
////////// } auxiliary (string) functions //////////


////////// class OptionParser { //////////
OptionParser::OptionParser() :
	_usage(_("%prog [options]")),
	_add_help_option(true),
	_add_version_option(true),
	_interspersed_args(true) {}

Option& OptionParser::add_option(const std::string& opt) {
	const std::string tmp[1] = { opt };
	return add_option(std::vector<std::string>(tmp, tmp + 1));
}
Option& OptionParser::add_option(const std::string& opt1, const std::string& opt2) {
	const std::string tmp[2] = { opt1, opt2 };
	return add_option(std::vector<std::string>(tmp, tmp + 2));
}
Option& OptionParser::add_option(const std::string& opt1, const std::string& opt2, const std::string& opt3) {
	const std::string tmp[3] = { opt1, opt2, opt3 };
	return add_option(std::vector<std::string>(tmp, tmp + 3));
}
Option& OptionParser::add_option(const std::vector<std::string>& v) {
	_opts.resize(_opts.size()+1);
	Option& option = _opts.back();
	std::string dest_fallback;
	for (std::vector<std::string>::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it->substr(0,2) == "--") {
			const std::string s = it->substr(2);
			if (option.dest() == "")
				option.dest(str_replace(s, "-", "_"));
			option._long_opts.insert(s);
			_optmap_l[s] = &option;
		} else {
			const std::string s = it->substr(1,1);
			if (dest_fallback == "")
				dest_fallback = s;
			option._short_opts.insert(s);
			_optmap_s[s] = &option;
		}
	}
	if (option.dest() == "")
		option.dest(dest_fallback);
	return option;
}

OptionParser& OptionParser::add_option_group(const OptionGroup& group) {
	for (std::list<Option>::const_iterator oit = group._opts.begin(); oit != group._opts.end(); ++oit) {
		const Option& option = *oit;
		for (std::set<std::string>::const_iterator it = option._short_opts.begin(); it != option._short_opts.end(); ++it)
			_optmap_s[*it] = &option;
		for (std::set<std::string>::const_iterator it = option._long_opts.begin(); it != option._long_opts.end(); ++it)
			_optmap_l[*it] = &option;
	}
	_groups.push_back(&group);
	return *this;
}

const Option& OptionParser::lookup_short_opt(const std::string& opt) const {
	optMap::const_iterator it = _optmap_s.find(opt);
	if (it == _optmap_s.end())
		error(_("no such option") + std::string(": -") + opt);
	return *it->second;
}

void OptionParser::handle_short_opt(const std::string& opt, const std::string& arg) {

	_remaining.pop_front();
	std::string value;

	const Option& option = lookup_short_opt(opt);
	if (option._nargs == 1) {
		value = arg.substr(2);
		if (value == "") {
			if (_remaining.empty())
				error("-" + opt + " " + _("option requires an argument"));
			value = _remaining.front();
			_remaining.pop_front();
		}
	} else {
		if (arg.length() > 2)
			_remaining.push_front(std::string("-") + arg.substr(2));
	}

	process_opt(option, std::string("-") + opt, value);
}

const Option& OptionParser::lookup_long_opt(const std::string& opt) const {

	std::list<std::string> matching;
	for (optMap::const_iterator it = _optmap_l.begin(); it != _optmap_l.end(); ++it) {
		if (it->first == opt)
			matching.push_back(it->first);
	}
	if (matching.size() > 1) {
		std::string x = str_join(", ", matching.begin(), matching.end());
		error(_("ambiguous option") + std::string(": --") + opt + " (" + x + "?)");
	}
	if (matching.size() == 0)
		error(_("no such option") + std::string(": --") + opt);

	return *_optmap_l.find(matching.front())->second;
}

void OptionParser::handle_long_opt(const std::string& optstr) {

	_remaining.pop_front();
	std::string opt, value;

	size_t delim = optstr.find("=");
	if (delim != std::string::npos) {
		opt = optstr.substr(0, delim);
		value = optstr.substr(delim+1);
	} else
		opt = optstr;

	const Option& option = lookup_long_opt(opt);
	if (option._nargs == 1 and delim == std::string::npos) {
		if (not _remaining.empty()) {
			value = _remaining.front();
			_remaining.pop_front();
		}
	}

	if (option._nargs == 1 and value == "")
		error("--" + opt + " " + _("option requires an argument"));

	process_opt(option, std::string("--") + opt, value);
}

Values& OptionParser::parse_args(const int argc, char const* const* const argv) {
	if (prog() == "")
		prog(basename(argv[0]));
	return parse_args(&argv[1], &argv[argc]);
}
Values& OptionParser::parse_args(const std::vector<std::string>& v) {

	_remaining.assign(v.begin(), v.end());

	if (add_version_option() and version() != "") {
		add_option("--version") .action("version") .help(_("show program's version number and exit"));
		_opts.splice(_opts.begin(), _opts, --(_opts.end()));
	}
	if (add_help_option()) {
		add_option("-h", "--help") .action("help") .help(_("show this help message and exit"));
		_opts.splice(_opts.begin(), _opts, --(_opts.end()));
	}

	while (not _remaining.empty()) {
		const std::string arg = _remaining.front();

		if (arg == "--") {
			_remaining.pop_front();
			break;
		}

		if (arg.substr(0,2) == "--") {
			handle_long_opt(arg.substr(2));
		} else if (arg.substr(0,1) == "-" and arg.length() > 1) {
			handle_short_opt(arg.substr(1,1), arg);
		} else {
			_remaining.pop_front();
			_leftover.push_back(arg);
			if (not interspersed_args())
				break;
		}
	}
	while (not _remaining.empty()) {
		const std::string arg = _remaining.front();
		_remaining.pop_front();
		_leftover.push_back(arg);
	}

	for (strMap::const_iterator it = _defaults.begin(); it != _defaults.end(); ++it) {
		if (not _values.is_set(it->first))
			_values[it->first] = it->second;
	}

	for (std::list<Option>::const_iterator it = _opts.begin(); it != _opts.end(); ++it) {
		if (it->get_default() != "" and not _values.is_set(it->dest()))
		    _values[it->dest()] = it->get_default();
	}

	return _values;
}

void OptionParser::process_opt(const Option& o, const std::string& opt, const std::string& value) {
	if (o.action() == "store") {
		std::string err = o.check_type(opt, value);
		if (err != "")
			error(err);
		_values[o.dest()] = value;
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "store_const") {
		_values[o.dest()] = o.get_const();
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "store_true") {
		_values[o.dest()] = "1";
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "store_false") {
		_values[o.dest()] = "0";
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "append") {
		std::string err = o.check_type(opt, value);
		if (err != "")
			error(err);
		_values[o.dest()] = value;
		_values.all(o.dest()).push_back(value);
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "append_const") {
		_values[o.dest()] = o.get_const();
		_values.all(o.dest()).push_back(o.get_const());
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "count") {
		_values[o.dest()] = str_inc(_values[o.dest()]);
		_values.is_set_by_user(o.dest(), true);
	}
	else if (o.action() == "help") {
		print_help();
		#ifdef ENABLE_MPI
			MPI_Finalize();
		#endif
		std::exit(EXIT_SUCCESS);
	}
	else if (o.action() == "version") {
		print_version();
		#ifdef ENABLE_MPI
			MPI_Finalize();
		#endif
		std::exit(EXIT_SUCCESS);
	}
	else if (o.action() == "callback" && o.callback()) {
		(*o.callback())(o, opt, value, *this);
	}
}

std::string OptionParser::format_option_help(unsigned int indent /* = 2 */) const {
	std::stringstream ss;

	if (_opts.empty())
		return ss.str();

	for (std::list<Option>::const_iterator it = _opts.begin(); it != _opts.end(); ++it) {
		if (it->help() != SUPPRESS_HELP)
			ss << it->format_help(indent);
	}

	return ss.str();
}

std::string OptionParser::format_help() const {
	std::stringstream ss;

	if (usage() != SUPPRESS_USAGE)
		ss << get_usage() << std::endl;

	if (description() != "")
		ss << str_format(description(), 0, cols()) << std::endl;

	ss << _("Options") << ":" << std::endl;
	ss << format_option_help();

	for (std::list<OptionGroup const*>::const_iterator it = _groups.begin(); it != _groups.end(); ++it) {
		const OptionGroup& group = **it;
		ss << std::endl << "  " << group.title() << ":" << std::endl;
		if (group.group_description() != "")
			ss << str_format(group.group_description(), 4, cols()) << std::endl;
		ss << group.format_option_help(4);
	}

	if (epilog() != "")
		ss << std::endl << str_format(epilog(), 0, cols());

	return ss.str();
}
void OptionParser::print_help() const {
	std::cout << format_help();
}

void OptionParser::set_usage(const std::string& u) {
	std::string lower = u;
	transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
	if (lower.compare(0, 7, "usage: ") == 0)
		_usage = u.substr(7);
	else
		_usage = u;
}
std::string OptionParser::format_usage(const std::string& u) const {
	std::stringstream ss;
	ss << _("Usage") << ": " << u << std::endl;
	return ss.str();
}
std::string OptionParser::get_usage() const {
	if (usage() == SUPPRESS_USAGE)
		return std::string("");
	return format_usage(str_replace(usage(), "%prog", prog()));
}
void OptionParser::print_usage(std::ostream& out) const {
	std::string u = get_usage();
	if (u != "")
		out << u << std::endl;
}
void OptionParser::print_usage() const {
	print_usage(std::cout);
}

std::string OptionParser::get_version() const {
	return str_replace(_version, "%prog", prog());
}
void OptionParser::print_version(std::ostream& out) const {
	out << get_version() << std::endl;
}
void OptionParser::print_version() const {
	print_version(std::cout);
}

void OptionParser::exit() const {
	std::ostringstream error_message;
	error_message << "OptionParser::exit() called" << std::endl;
	MARDYN_EXIT(error_message.str());
}
void OptionParser::error(const std::string& msg) const {
	print_usage(std::cerr);
	std::ostringstream error_message;
	error_message << prog() << ": " << _("error") << ": " << msg << std::endl;
	MARDYN_EXIT(error_message.str());
}
////////// } class OptionParser //////////

////////// class Values { //////////
const std::string& Values::operator[] (const std::string& d) const {
	strMap::const_iterator it = _map.find(d);
	static const std::string empty = "";
	return (it != _map.end()) ? it->second : empty;
}
void Values::is_set_by_user(const std::string& d, bool yes) {
	if (yes)
		_userSet.insert(d);
	else
		_userSet.erase(d);
}
////////// } class Values //////////

////////// class Option { //////////
std::string Option::check_type(const std::string& opt, const std::string& val) const {
	std::istringstream ss(val);
	std::stringstream err;

	if (type() == "int" || type() == "long") {
		long t;
		if (not (ss >> t))
			err << _("option") << " " << opt << ": " << _("invalid integer value") << ": '" << val << "'";
	}
	else if (type() == "float" || type() == "double") {
		double t;
		if (not (ss >> t))
			err << _("option") << " " << opt << ": " << _("invalid floating-point value") << ": '" << val << "'";
	}
	else if (type() == "choice") {
		if (find(choices().begin(), choices().end(), val) == choices().end()) {
			std::list<std::string> tmp = choices();
			transform(tmp.begin(), tmp.end(), tmp.begin(), str_wrap("'"));
			err << _("option") << " " << opt << ": " << _("invalid choice") << ": '" << val << "'"
				<< " (" << _("choose from") << " " << str_join(", ", tmp.begin(), tmp.end()) << ")";
		}
	}
	else if (type() == "complex") {
		std::complex<double> t;
		if (not (ss >> t))
			err << _("option") << " " << opt << ": " << _("invalid complex value") << ": '" << val << "'";
	}

	return err.str();
}

std::string Option::format_option_help(unsigned int indent /* = 2 */) const {

	std::string mvar_short, mvar_long;
	if (nargs() == 1) {
		std::string mvar = metavar();
		if (mvar == "") {
			mvar = type();
			transform(mvar.begin(), mvar.end(), mvar.begin(), ::toupper);
		 }
		mvar_short = " " + mvar;
		mvar_long = "=" + mvar;
	}

	std::stringstream ss;
	ss << std::string(indent, ' ');

	if (not _short_opts.empty()) {
		ss << str_join_trans(", ", _short_opts.begin(), _short_opts.end(), str_wrap("-", mvar_short));
		if (not _long_opts.empty())
			ss << ", ";
	}
	if (not _long_opts.empty())
		ss << str_join_trans(", ", _long_opts.begin(), _long_opts.end(), str_wrap("--", mvar_long));

	return ss.str();
}

std::string Option::format_help(unsigned int indent /* = 2 */) const {
	std::stringstream ss;
	std::string h = format_option_help(indent);
	unsigned int width = cols();
	unsigned int opt_width = std::min(width*3/10, 36u);
	bool indent_first = false;
	ss << h;
	// if the option list is too long, start a new paragraph
	if (h.length() >= (opt_width-1)) {
		ss << std::endl;
		indent_first = true;
	} else {
		ss << std::string(opt_width - h.length(), ' ');
		if (help() == "")
			ss << std::endl;
	}
	if (help() != "") {
		std::string help_str = (get_default() != "") ? str_replace(help(), "%default", get_default()) : help();
		ss << str_format(help_str, opt_width, width, indent_first);
	}
	return ss.str();
}

Option& Option::action(const std::string& a) {
	_action = a;
	if (a == "store_const" || a == "store_true" || a == "store_false" ||
	    a == "append_const" || a == "count" || a == "help" || a == "version")
		nargs(0);
	return *this;
}
////////// } class Option //////////

}
