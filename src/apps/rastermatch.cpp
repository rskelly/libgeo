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

#include "grid.hpp"

using namespace geo::grid;

float memAvg(std::vector<float>& buf, int size, float nd) {
	float sum = 0, div = 0;
	float rad = (float) size / 2.0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			float v = buf[r * size + c];
			float d = std::sqrt(std::pow(r - size / 2.0, 2.0) + std::pow(c - size / 2.0, 2.0));
			if(d <= rad && v != nd) {
				float w = 1.0 - d / rad;
				sum += v * w;
				div += w;
			}
		}
	}
	return div > 0 ? sum / div : nd;
}

float memMed(std::vector<float>& buf, float nd) {
	std::vector<float> lst;
	for(float d : buf) {
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

void doInterp(Band<float>& tmem, Band<float>& amem, Band<float>& difmem, int size) {

	const GridProps& tprops = tmem.props();
	const GridProps& aprops = amem.props();

	float nd = tprops.nodata();
	if(size % 2 == 0) ++size;

	std::vector<float> buf(size * size);

	difmem.fill(nd);

	float tavg, aavg, tv;
	for(int trow = 0; trow < tprops.rows(); ++trow) {
		for(int tcol = 0; tcol < tprops.cols(); ++tcol) {

			if((tv = tmem.get(tcol, trow)) == nd)
				continue;

			float x = tprops.toX(tcol);
			float y = tprops.toY(trow);

			int acol = aprops.toCol(x);
			int arow = aprops.toRow(y);

			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int cc = acol - size / 2 + c;
					int rr = arow - size / 2 + r;
					if(aprops.hasCell(cc, rr)) {
						buf[r * size + c] = amem.get(cc, rr);
					} else {
						buf[r * size + c] = nd;
					}
				}
			}
			aavg = memAvg(buf, size, nd);

			if(aavg == nd) {
				difmem.set(tcol, trow, tv);
				continue;
			}

			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int cc = tcol - size / 2 + c;
					int rr = trow - size / 2 + r;
					if(tprops.hasCell(cc, rr)) {
						buf[r * size + c] = tmem.get(cc, rr);
					} else {
						buf[r * size + c] = nd;
					}
				}
			}
			tavg = memAvg(buf, size, nd);

			difmem.set(tcol, trow, tv + (aavg - tavg));
		}
	}

	for(int trow = 0; trow < tprops.rows(); ++trow) {
		for(int tcol = 0; tcol < tprops.cols(); ++tcol)
			tmem.set(tcol, trow, difmem.get(tcol, trow));
	}
}

float meanDif(std::vector<float>& tvec, float tn, std::vector<float>& avec, float an, int size) {
	float av, tv, sum = 0;
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

void doInterp2(Band<float>& tmem, Band<float>& amem, Band<float>& dmem, int size) {

	const GridProps& tprops = tmem.props();
	const GridProps& aprops = amem.props();

	float tn = tprops.nodata();
	float an = aprops.nodata();

	if(size % 2 != 0) ++size;
	int side = 2 * size + 1;

	dmem.fill(tn);

	std::vector<float> avec(side * side);
	std::vector<float> tvec(side * side);

	int tcols = tprops.cols();
	int trows = tprops.rows();
	float invalid = std::nan("");

	int minp = -1;
	float tv;
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

			if((tv = tmem.get(tcol, trow)) != tn) {

				if((trow - lastRow) != 1) {
					tmem.writeToVector(tvec, tcol - size, trow - size, side, side, invalid);
					amem.writeToVector(avec, tcol - size, trow - size, side, side, invalid);
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
						std::memcpy(avec.data() + vr * side + vc, amem.grid() + rr * tcols + rc, cols * sizeof(float));
						std::memcpy(tvec.data() + vr * side + vc, tmem.grid() + rr * tcols + rc, cols * sizeof(float));
					}
				}
				lastRow = trow;

				float dif = meanDif(tvec, tn, avec, an, side);

				dmem.set(tcol, trow, tv + dif);
			} else {
				dmem.set(tcol, trow, tn);
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
	Band<float> tmem;
	{
		Band<float> trast(target, tband - 1, false, true);
		tprops = trast.props();
		tprops.setWritable(true);
		tprops.setBands(1);
		tmem.init(tprops, mapped);
		trast.writeTo(tmem);
	}

	GridProps mprops;
	Band<float> mmem;
	if(!mask.empty()) {
		Band<float> mrast(mask, mband - 1, false, true);
		mprops = mrast.props();
		mprops.setBands(1);
		mmem.init(mprops, mapped);
		mrast.writeTo(mmem);
	}

	GridProps aprops(tprops);
	Band<float> amem;
	{
		amem.init(tprops, mapped);
		amem.fill(tprops.nodata());

		for(size_t i = 0; i < anchors.size(); ++i) {
			Band<float> cmem;
			GridProps cprops;
			{
				Band<float> crast(anchors[i], abands[i] - 1, false, true);
				cprops = crast.props();
				cprops.setBands(1);
				cprops.setWritable(true);
				cmem.init(cprops, mapped);
				crast.writeTo(cmem);
			}

			float cv, cn = cprops.nodata();
			//float tn = tprops.nodata();
			bool hasMask = !mask.empty();
			for(int crow = 0; crow < cprops.rows(); ++crow) {
				for(int ccol = 0; ccol < cprops.cols(); ++ccol) {
					float x = cprops.toX(ccol);
					float y = cprops.toY(crow);
					int mcol = mprops.toCol(x);
					int mrow = mprops.toRow(y);
					int tcol = tprops.toCol(x);
					int trow = tprops.toRow(y);
					if(tprops.hasCell(tcol, trow)) {
							if((cv = cmem.get(ccol, crow)) != cn
									&& (!hasMask || mmem.get<int>(mcol, mrow) == 1)) {
							amem.set(tcol, trow, cv);
						}
					}
				}
			}
		}
	}

	{
		Band<float> dmem(tprops, mapped);
		for(const int& d : sizes) {

			doInterp2(tmem, amem, dmem, d);

			if(test) {
				std::stringstream tmp;
				tmp << "/tmp/match_" <<  d << ".tif";
				std::cerr << tmp.str() << "\n";
				Band<float> rtmp(tmp.str(), tprops);
				tmem.writeTo(rtmp);
			}
		}
	}

	Band<float> adjrast(adjusted, tprops);
	tmem.writeTo(adjrast);

}


