/*
 * rastermatch.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: rob
 */
#include <vector>
#include <sstream>

#include "raster.hpp"

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

double meanDif(std::vector<double>& tmem, std::vector<double>& amem, int tcol, int trow, int tcols, int trows, int size, double tn, double an) {
	double av, tv;
	double sum = 0;
	int count = 0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			int cc = tcol + c - size / 2;
			int rr = trow + r - size / 2;
			if(cc < 0 || rr < 0 || rr >= trows || cc >= tcols)
				continue;
			tv = tmem[rr * tcols + cc];
			av = amem[rr * tcols + cc];
			if(tv != tn && av != an)
				sum += (av - tv);
		}
	}
	count = size * size;
	return count > 0 ? sum / count : 0;
}

void doInterp2(std::vector<double>& tmem, const GridProps& tprops, std::vector<double>& amem, const GridProps& aprops, std::vector<double>& difmem, int size) {

	double tn = tprops.nodata();
	double an = aprops.nodata();
	int tcols = tprops.cols();
	int trows = tprops.rows();

	if(size % 2 == 0) ++size;

	std::fill(difmem.begin(), difmem.end(), tn);

	double tv;
	for(int trow = 0; trow < tprops.rows(); ++trow) {
		if(trow % 100 == 0)
			std::cout << "Row: " << trow << " of " << tprops.rows() << "\n";
		for(int tcol = 0; tcol < tprops.cols(); ++tcol) {

			if((tv = tmem[trow * tcols + tcol]) == tn)
				continue;

			double dif = meanDif(tmem, amem, tcol, trow, tcols, trows, size, tn, an);

			difmem[trow * tcols + tcol] = tv + dif;
		}
	}

	tmem.swap(difmem);
}

int main(int argc, char** argv) {

	std::vector<std::string> anchors;
	std::vector<int> abands;
	std::string target;
	int tband = 1;
	std::string adjusted;
	std::string adjustment;
	std::string mask;
	int mband = 1;
	std::vector<int> sizes;
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
		} else if(arg == "-f") {
			adjustment = argv[++i];
			continue;
		} else if(arg == "-m") {
			mask = argv[++i];
			mband = atoi(argv[++i]);
			continue;
		} else if(arg == "-d") {
			std::stringstream ss(argv[++i]);
			std::string part;
			while(std::getline(ss, part, ','))
				sizes.push_back(atoi(part.c_str()));
			continue;
		} else if(arg == "-p") {
			mapped = true;
			continue;
		}
	}

	Bounds bounds;
	GridProps tprops;
	MemRaster tmem;
	{
		Raster trast(target);
		tprops = trast.props();
		tprops.setWritable(true);
		tprops.setBands(1);
		tmem.init(tprops, mapped);
		tmem.fillFloat(tprops.nodata(), 1);
		trast.writeTo(tmem, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);
		bounds.extend(tprops.bounds());
	}

	GridProps mprops;
	MemRaster mmem;
	if(!mask.empty()) {
		Raster mrast(mask);
		mprops = mrast.props();
		mmem.init(mprops, mapped);
		mmem.fillFloat(mprops.nodata(), 1);
		mrast.writeTo(mmem, mprops.cols(), mprops.rows(), 0, 0, 0, 0, mband, 1);
	}

	GridProps aprops(tprops);
	MemRaster amem;
	{
		bool first = true;
		for(size_t i = 0; i < anchors.size(); ++i) {
			MemRaster cmem;
			GridProps cprops;
			{
				Raster crast(anchors[i]);
				cprops = crast.props();
				cprops.setBands(1);
				cprops.setWritable(true);
				cmem.init(cprops, mapped);
				cmem.fillFloat(cprops.nodata(), 1);
				if(first) {
					aprops.setNoData(cprops.nodata());
					first = false;
				}
				crast.writeTo(cmem, cprops.cols(), cprops.rows(), 0, 0, 0, 0, abands[i], 1);
			}

			amem.init(aprops, mapped);
			amem.fillFloat(aprops.nodata(), 1);

			double cv, cn = cprops.nodata();
			bool hasMask = !mask.empty();
			for(int crow = 0; crow < cprops.rows(); ++crow) {
				for(int ccol = 0; ccol < cprops.cols(); ++ccol) {
					double x = cprops.toCentroidX(ccol);
					double y = cprops.toCentroidY(crow);
					int tcol = tprops.toCol(x);
					int trow = tprops.toRow(y);
					int mcol = mprops.toCol(x);
					int mrow = mprops.toRow(y);
					if(tprops.hasCell(tcol, trow)) {
							if((cv = cmem.getFloat(ccol, crow, 1)) != cn
									&& (!hasMask || mmem.getInt(mcol, mrow, 1) == 1)) {
							amem.setFloat(tcol, trow, cv, 1);
						} else {
							amem.setFloat(tcol, trow, cn, 1);
						}
					}
				}
			}
		}
	}

	{
		std::vector<double> tvec(tprops.cols() * tprops.rows());
		{
			std::vector<double> dvec(tvec.size());
			std::vector<double> avec(dvec.size());
			tmem.readIntoVector(tvec);
			amem.readIntoVector(avec);
			for(int d : sizes) {
				std::cout << "Interp at " << d << "\n";
				doInterp2(tvec, tprops, avec, aprops, dvec, d);
			}
		}
		//MemRaster difmem(tprops, mapped);
		//difmem.writeFromVector(dvec);
		tmem.writeFromVector(tvec);
	}

	Raster adjrast(adjusted, tprops);
	tmem.writeTo(adjrast);

}


