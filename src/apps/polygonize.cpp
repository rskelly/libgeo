/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

void usage() {
	std::cerr << "Usage: polygonize <input raster> <output vector> <output layer> [field name (dn)]\n"
			<< " -f <format>     The output format. Any of the available\n"
			<< "                 OGR drivers. Default SQLite.\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -b <band>       The band. Default 1.\n"
			<< " -d              Remove dangles.\n"
			<< " -h              Remove holes.\n"
			<< " -t <t>          The number of threads. Default 3. There are 2 threads for\n"
			<< "                 reading and writing and n threads for processing polygons.\n"
			<< " -m <mask image> A mask image to determine the bounds of the input.\n";
}

class DummyStatus : public Status {
private:
	geo::util::Callbacks *m_callbacks;
	float m_start, m_end;
public:
	DummyStatus() : Status(nullptr, 0, 1),
		m_callbacks(nullptr), m_start(0), m_end(0) {}
	void update(float s, const std::string& msg = "") {
		std::cerr << (s * 100) << "% " << msg << "\n";
	}
	float start() const {
		return 0;
	}
	float end() const {
		return 1;
	}
	geo::util::Callbacks* callbacks() const {
		return nullptr;
	}
};

int main(int argc, char** argv) {

	using namespace geo::raster;

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string driver = "SQLite";
	uint16_t srid = 0;
	uint16_t band = 1;
	uint16_t threads = 1;
	bool holes = false;
	bool dangles = false;
	std::vector<std::string> args;
	std::string mask;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-f") {
			driver = argv[++i];
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-h") {
			holes = true;
		} else if(v == "-d") {
			dangles = true;
		} else if(v == "-t") {
			threads = atoi(argv[++i]);
		} else if(v == "-m") {
			mask = argv[++i];
		} else {
			args.push_back(argv[i]);
		}
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	if(args.size() < 3)
		args.push_back("dn");

	if(threads <= 0) {
		std::cerr << "Illegal thread value. Defaulting to 1.\n";
		threads = 1;
	}

	DummyStatus status;
	bool cancel;

	Raster test(args[0]);
	test.polygonize(args[1], args[2], driver, srid, band, holes, dangles, "", 0, threads, cancel, status);

	return 0;
}


