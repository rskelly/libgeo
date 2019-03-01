/*
 * rastermatch.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: rob
 */
#include <vector>
#include <sstream>
#include <thread>

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

double meanDif(std::vector<double>& tvec, double tn, std::vector<double>& avec, double an, int size) {
	double av, tv, sum = 0;
	int count = 0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			tv = tvec[r * size + c];
			av = avec[r * size + c];
			if(isnan(tv) || isnan(av))
				continue;
			++count;
			if(tv == tn || av == an)
				continue;
			sum += (av - tv);
		}
	}
	return count > 0 ? sum / count : 0;
}

void doInterp2(MemRaster* tmem, MemRaster* amem, MemRaster* dmem, int size) {

	const GridProps& tprops = tmem->props();
	const GridProps& aprops = amem->props();

	double tn = tprops.nodata();
	double an = aprops.nodata();

	if(size % 2 == 0) ++size;

	dmem->fillFloat(tn, 1);

	std::vector<double> avec(size * size);
	std::vector<double> tvec(size * size);

	double tv;
	for(int trow = 0; trow < tprops.rows(); ++trow) {
		if(trow % 100 == 0)
			std::cout << "Row: " << trow << " of " << tprops.rows() << "\n";
		for(int tcol = 0; tcol < tprops.cols(); ++tcol) {

			if(trow == 321 && tcol == 100)
				std::cerr << "c\n";

			tmem->writeToVector(tvec, tcol - size / 2, trow - size / 2, size, size, 1);
			amem->writeToVector(avec, tcol - size / 2, trow - size / 2, size, size, 1);

			double dif = meanDif(tvec, tn, avec, an, size);

			if((tv = tmem->getFloat(tcol, trow, 1)) != tn) {
				dmem->setFloat(tcol, trow, tv + dif, 1);
			} else {
				dmem->setFloat(tcol, trow, tn, 1);
			}
		}
	}

	dmem->writeTo(*tmem);
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
	bool test = true;
	bool mapped = false;
	int threads = 4;

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
			for(int& s : sizes) {
				if(s % 2 == 0) s++;
			}
			continue;
		}
	}

	Bounds bounds;
	GridProps tprops;
	int maxSize = sizes[sizes.size() - 1];
	int bufSize = (int) std::ceil((float) maxSize / 2.0);
	int rowHeight = 0;
	int bufHeight = 0;
	std::vector<MemRaster> tmems;
	{
		Raster trast(target);
		tprops = trast.props();
		tprops.setWritable(true);
		tprops.setBands(1);
		const Bounds& bounds = tprops.bounds();
		rowHeight = (int) std::ceil((float) tprops.rows() / threads);
		bufHeight = rowHeight + bufSize * 2;
		for(int i = 0; i < threads; ++i) {
			int row = i * rowHeight - bufSize;
			double y0 = tprops.toCentroidY(i * rowHeight) - tprops.resolutionY() * bufSize;
			double y1 = tprops.toCentroidY((i + 1) * rowHeight) + tprops.resolutionY() * bufSize;
			Bounds bounds0;
			bounds0.extend(bounds.minx(), y0);
			bounds0.extend(bounds.maxx(), y1);
			GridProps tprops0(tprops);
			tprops0.setBounds(bounds0);
			tmems.emplace_back(tprops0, mapped);
			trast.writeTo(tmems[i], tprops0.cols(), bufHeight, 0, row, 0, 0, tband, 1);
		}
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
	std::vector<MemRaster> amems;
	{
		for(int i = 0; i < threads; ++i) {
			amems.emplace_back(tmems[i].props(), mapped);
			amems[i].fillFloat(amems[i].props().nodata(), 1);
		}

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
			bool hasMask = !mask.empty();
			for(int crow = 0; crow < cprops.rows(); ++crow) {
				for(int ccol = 0; ccol < cprops.cols(); ++ccol) {
					double x = cprops.toCentroidX(ccol);
					double y = cprops.toCentroidY(crow);
					int mcol = mprops.toCol(x);
					int mrow = mprops.toRow(y);
					for(int i = 0; i < threads; ++i) {
						const GridProps& tprops0 = tmems[i].props();
						int tcol = tprops0.toCol(x);
						int trow = tprops0.toRow(y);
						if(tprops0.hasCell(tcol, trow)) {
								if((cv = cmem.getFloat(ccol, crow, 1)) != cn
										&& (!hasMask || mmem.getInt(mcol, mrow, 1) == 1)) {
								amems[i].setFloat(tcol, trow, cv, 1);
							} else {
								amems[i].setFloat(tcol, trow, 0, 1);
							}
						}
					}
				}
			}
		}
	}

	{
		std::vector<std::thread> threadv;
		std::vector<MemRaster> dmems;
		for(int i = 0; i < threads; ++i)
			dmems.emplace_back(tmems[i].props(), mapped);
		for(const int& d : sizes) {
			for(int i = 0; i < threads; ++i)
				threadv.emplace_back(doInterp2, &(tmems[i]), &(amems[i]), &(dmems[i]), d);
			for(int i = 0; i < threads; ++i) {
				if(threadv[i].joinable())
					threadv[i].join();
			}
			if(test) {
				std::stringstream tmp;
				tmp << "/tmp/match_" <<  d << ".tif";
				std::cerr << tmp.str() << "\n";
				Raster rtmp(tmp.str(), tprops);
				for(int i = 0; i < threads; ++i) {
					const GridProps& props = tmems[i].props();
					tmems[i].writeTo(rtmp, props.cols(), rowHeight, 0, bufSize, 0, i * rowHeight, 1, 1);
				}
			}
		}
	}

	Raster adjrast(adjusted, tprops);
	for(int i = 0; i < threads; ++i) {
		const GridProps& props = tmems[i].props();
		tmems[i].writeTo(adjrast, props.cols(), rowHeight, 0, bufSize, 0, i * rowHeight, 1, 1);
	}

}


