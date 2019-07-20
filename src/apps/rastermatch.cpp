/*
 * rastermatch.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: rob
 */
#include <grid.hpp>
#include <vector>
#include <sstream>
#include <thread>
#include <iostream>


using namespace geo::raster;

double memAvg(std::vector<double>& buf, int size, double nd) {
	double sum = 0, div = 0;
	double rad = (double) size / 2.0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			double v = buf[r * size + c];
			double d = std::sqrt(std::pow(r - size / 2.0, 2.0) + std::pow(c - size / 2.0, 2.0));
			if(d <= rad && v != nd) {
				double w = 1.0 - d / rad;
				sum += v * w;
				div += w;
			}
		}
	}
	return div > 0 ? sum / div : nd;
}

double memMed(std::vector<double>& buf, double nd) {
	std::vector<double> lst;
	for(double d : buf) {
		if(d != nd)
			lst.push_back(d);
	}
	if(lst.empty()) {
		return nd;
	} else {
		std::sort(lst.begin(), lst.end());
		if(lst.size() % 2 == 0) {
			size_t i = lst.size() / 2;
			return (lst[i] + lst[i + 1]) / 2.0;
		} else {
			return lst[lst.size() / 2];
		}
	}
}

void doInterp(MemRaster& tmem, MemRaster& amem, MemRaster& difmem, int size) {

	const GridProps& tprops = tmem.props();
	const GridProps& aprops = amem.props();

	double nd = tprops.nodata();
	if(size % 2 == 0) ++size;

	std::vector<double> buf(size * size);

	difmem.fillFloat(nd, 1);

	double tavg, aavg, tv;
	for(int trow = 0; trow < tprops.rows(); ++trow) {
		for(int tcol = 0; tcol < tprops.cols(); ++tcol) {

			if((tv = tmem.getFloat(tcol, trow, 1)) == nd)
				continue;

			double x = tprops.toCentroidX(tcol);
			double y = tprops.toCentroidY(trow);

			int acol = aprops.toCol(x);
			int arow = aprops.toRow(y);

			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int cc = acol - size / 2 + c;
					int rr = arow - size / 2 + r;
					if(aprops.hasCell(cc, rr)) {
						buf[r * size + c] = amem.getFloat(cc, rr, 1);
					} else {
						buf[r * size + c] = nd;
					}
				}
			}
			aavg = memAvg(buf, size, nd);

			if(aavg == nd) {
				difmem.setFloat(tcol, trow, tv, 1);
				continue;
			}

			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int cc = tcol - size / 2 + c;
					int rr = trow - size / 2 + r;
					if(tprops.hasCell(cc, rr)) {
						buf[r * size + c] = tmem.getFloat(cc, rr, 1);
					} else {
						buf[r * size + c] = nd;
					}
				}
			}
			tavg = memAvg(buf, size, nd);

			difmem.setFloat(tcol, trow, tv + (aavg - tavg), 1);
		}
	}

	for(int trow = 0; trow < tprops.rows(); ++trow) {
		for(int tcol = 0; tcol < tprops.cols(); ++tcol)
			tmem.setFloat(tcol, trow, difmem.getFloat(tcol, trow, 1), 1);
	}
}

double meanDif(std::vector<double>& tvec, double tn, std::vector<double>& avec, double an, int size) {
	double av, tv, sum = 0;
	int count = 0;
	for(int i = 0; i < size * size; ++i) {
		tv = tvec[i];
		av = avec[i];
		if(isnan(tv) || isnan(av))
			continue;
		if(tv == tn || tv == 0) //  This shouldn't happen.
			continue;
		++count;
		if(av == an)
			continue;
		sum += (av - tv);
	}
	return count > 0 ? sum / count : 0;
}

void doInterp2(MemRaster& tmem, MemRaster& amem, MemRaster& dmem, int size) {

	const GridProps& tprops = tmem.props();
	const GridProps& aprops = amem.props();

	double tn = tprops.nodata();
	double an = aprops.nodata();

	if(size % 2 != 0) ++size;
	int side = 2 * size + 1;

	dmem.fillFloat(tn, 1);

	std::vector<double> avec(side * side);
	std::vector<double> tvec(side * side);

	int tcols = tprops.cols();
	int trows = tprops.rows();
	double invalid = std::nan("");

	int minp = -1;
	double tv;
	for(int tcol = 0; tcol < tcols; ++tcol) {

		int p = (int) (((float) tcol / tcols) * 100.0);
		if(p > minp) {
			if(p % 25 == 0)
				std::cout << " " << p << "% ";
			if(p % 10 == 0)
				std::cout << ".";
			minp = p;
			std::cout << std::flush;
		}

		int lastRow = 0;

		for(int trow = 0; trow < trows; ++trow) {

			if((tv = tmem.getFloat(tcol, trow, 1)) != tn) {

				if((trow - lastRow) != 1) {
					tmem.writeToVector(tvec, tcol - size, trow - size, side, side, 1, invalid);
					amem.writeToVector(avec, tcol - size, trow - size, side, side, 1, invalid);
				} else {
					int cols = side;
					int rc = tcol - size, rr = trow + size; // Only interested in the new bottom row. All others stay the same.
					int vc = 0, vr = (trow - 1) % side;
					for(int i = 0; i < side; ++i) {
						avec[vr * side + i] = invalid;
						tvec[vr * side + i] = invalid;
					}
					if(rr < trows) {
						if(rc < 0) {
							cols += rc;
							vc -= rc;
							rc = 0;
						}
						if(rc + cols > tcols)
							cols = tcols - rc;
						std::memcpy(avec.data() + vr * side + vc, ((double*) amem.grid()) + rr * tcols + rc, cols * sizeof(double));
						std::memcpy(tvec.data() + vr * side + vc, ((double*) tmem.grid()) + rr * tcols + rc, cols * sizeof(double));
					}
				}
				lastRow = trow;

				double dif = meanDif(tvec, tn, avec, an, side);

				dmem.setFloat(tcol, trow, tv + dif, 1);
			} else {
				dmem.setFloat(tcol, trow, tn, 1);
			}
		}
	}
	std::cout << "100%\n";

	dmem.writeTo(tmem);
}

int main(int argc, char** argv) {

	std::vector<std::string> anchors;
	std::vector<int> abands;
	std::string target;
	int tband = 1;
	std::string adjusted;
	std::string mask;
	int mband = 1;
	std::vector<int> sizes;
	bool test = true;
	bool mapped = false;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-a") {
			anchors.push_back(argv[++i]);
			abands.push_back(atoi(argv[++i]));
			continue;
		} else if(arg == "-t") {
			target = argv[++i];
			tband = atoi(argv[++i]);
			continue;
		} else if(arg == "-o") {
			adjusted = argv[++i];
			continue;
		} else if(arg == "-m") {
			mask = argv[++i];
			mband = atoi(argv[++i]);
			continue;
		} else if(arg == "-d") {
			std::string val = argv[++i];
			if(val.find(',') != std::string::npos) {
				std::stringstream ss(val);
				std::string part;
				while(std::getline(ss, part, ','))
					sizes.push_back(atoi(part.c_str()));
				std::sort(sizes.begin(), sizes.end());
			} else {
				sizes.push_back(atoi(val.c_str()));
			}
			continue;
		}
	}

	std::sort(sizes.begin(), sizes.end());
	std::reverse(sizes.begin(), sizes.end());

	Bounds bounds;
	GridProps tprops;
	MemRaster tmem;
	{
		Raster trast(target);
		tprops = trast.props();
		tprops.setWritable(true);
		tprops.setBands(1);
		tmem.init(tprops, mapped);
		trast.writeTo(tmem, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);
	}

	GridProps mprops;
	MemRaster mmem;
	if(!mask.empty()) {
		Raster mrast(mask);
		mprops = mrast.props();
		mmem.init(mprops, mapped);
		mrast.writeTo(mmem, mprops.cols(), mprops.rows(), 0, 0, 0, 0, mband, 1);
	}

	GridProps aprops(tprops);
	MemRaster amem;
	{
		amem.init(tprops, mapped);
		amem.fillFloat(tprops.nodata(), 1);

		for(size_t i = 0; i < anchors.size(); ++i) {
			MemRaster cmem;
			GridProps cprops;
			{
				Raster crast(anchors[i]);
				cprops = crast.props();
				cprops.setBands(1);
				cprops.setWritable(true);
				cmem.init(cprops, mapped);
				crast.writeTo(cmem, cprops.cols(), cprops.rows(), 0, 0, 0, 0, abands[i], 1);
			}

			double cv, cn = cprops.nodata();
			//double tn = tprops.nodata();
			bool hasMask = !mask.empty();
			for(int crow = 0; crow < cprops.rows(); ++crow) {
				for(int ccol = 0; ccol < cprops.cols(); ++ccol) {
					double x = cprops.toCentroidX(ccol);
					double y = cprops.toCentroidY(crow);
					int mcol = mprops.toCol(x);
					int mrow = mprops.toRow(y);
					int tcol = tprops.toCol(x);
					int trow = tprops.toRow(y);
					if(tprops.hasCell(tcol, trow)) {
							if((cv = cmem.getFloat(ccol, crow, 1)) != cn
									&& (!hasMask || mmem.getInt(mcol, mrow, 1) == 1)) {
							amem.setFloat(tcol, trow, cv, 1);
						}
					}
				}
			}
		}
	}

	{
		MemRaster dmem(tprops, mapped);
		for(const int& d : sizes) {

			doInterp2(tmem, amem, dmem, d);

			if(test) {
				std::stringstream tmp;
				tmp << "/tmp/match_" <<  d << ".tif";
				std::cerr << tmp.str() << "\n";
				Raster rtmp(tmp.str(), tprops);
				tmem.writeTo(rtmp);
			}
		}
	}

	Raster adjrast(adjusted, tprops);
	tmem.writeTo(adjrast);

}


