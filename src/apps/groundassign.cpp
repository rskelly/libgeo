#include <fstream>
#include <iomanip>
#include <liblas/liblas.hpp>

int main(int argc, char** argv) {
	
	std::string infile = argv[1];
	std::string gndfile = argv[2];
	std::string outfile = argv[3];

	liblas::ReaderFactory rf;

	std::cout << "opening " << infile << "\n";
	std::ifstream instr(infile, std::ios::binary|std::ios::in);
	liblas::Reader inrdr = rf.CreateWithStream(instr);

	std::cout << "opening " << gndfile << "\n";
	std::ifstream gndstr(gndfile, std::ios::binary|std::ios::in);
	liblas::Reader gndrdr = rf.CreateWithStream(gndstr);

	std::cout << "opening " << outfile << "\n";
	std::ofstream outstr(outfile, std::ios::binary|std::ios::out);
	liblas::Header outhdr(inrdr.GetHeader());
	liblas::Writer outwtr(outstr, outhdr);

	long ground = 0;
	long a = 0;
	long b = 0;
	bool found = false;
	liblas::Classification cls(2, false, false, false);
	while(gndrdr.ReadNextPoint()) {
		++a;
		const liblas::Point& gpt = gndrdr.GetPoint();
		double gx = (int) (gpt.GetX() * 100);
		double gy = (int) (gpt.GetY() * 100);
		double gz = (int) (gpt.GetZ() * 100);
		while(inrdr.ReadNextPoint()) {
			++b;
			liblas::Point pt(inrdr.GetPoint());
			double x = (int) (pt.GetX() * 100);
			double y = (int) (pt.GetY() * 100);
			double z = (int) (pt.GetZ() * 100);
			//std::cout << std::setprecision(12) << a << " - " << b << " " << pt.GetX() << " " << gpt.GetX() << " " << pt.GetY() << " " << gpt.GetY() << " " << pt.GetZ() << " " << gpt.GetZ() << "\n";
			//std::cout << (x == gx) << " " <<  (y == gy) << " " <<  (z = gz) << "\n";
			if(x == gx && y == gy && z == gz){
				pt.SetClassification(cls);
				//std::cout << ((int) pt.GetClassification().GetClass()) << "\n";
				++ground;
				found = true;
			}
			//std::cin.get();
			outwtr.WritePoint(pt);
			if(found) {
				found = false;
				break;
			}
		}
	}

}