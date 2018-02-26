/*
 * rast_polygonize_tacet.cpp
 *
 *  Created on: Feb 23, 2018
 *      Author: rob
 */


bool computeDirection(int tl, int tr, int bl, int br, int target, int& c, int& r, int& dir) {
	int q = (tl == target && tr == target ? 1 : 0) | (bl == target && br == target ? 2 : 0)
			| (tl == target && bl == target ? 4 : 0) | (tr == target && br == target ? 8 : 0);

	int p = (tl == target ? 1 : 0) | (tr == target ? 2 : 0) | (bl == target ? 4 : 0) | (br == target ? 8 : 0);

	// dir:
	// 0 - none
	// 1 - right
	// 2 - up
	// 4 - left
	// 8 - down

	switch(q) {
	case 1:
	case 9:
		switch(dir) {
		case 2: ++c; dir = 1; return true;
		case 8: --c; dir = 4; return true;
		}
		break;
	case 2:
	case 6:
		switch(dir) {
		case 1: ++r; dir = 8; return true;
		case 4: --r; dir = 2; return true;
		}
		break;
	case 4:
	case 5: ++r; dir = 8; return true;
	case 8:
	case 10: --r; dir = 2; return true;
	case 0:
		switch(p) {
		case 1: --c; dir = 4; return true;
		case 2: --r; dir = 2; return true;
		case 4: ++r; dir = 8; return true;
		case 8: ++c; dir = 1; return true;
		case 6:
		case 9: ++c; dir = 1; return true; // Default to right when crossed.
		}
	}

	/*
	std::cerr << "tl: " << tl << ", tr: " << tr << "\n;";
	std::cerr << "bl: " << bl << ", br: " << br << "\n;";
	std::cerr << "target: " << target << "\n";
	std::cerr << "c : " << c << ", r: " << r << "\n";
	g_runerr("Failed to determine a direction.");
	*/
	return false;
}

void rightCoord(int& c, int& r, int dir) {
	switch(dir) {
	case 1: break;
	case 2: --r; break;
	case 4: --r; --c; break;
	case 8: --c; break;
	}
}

/*
template <class T>
T getXatD(T x, int d) {
	switch(d % 8) {
	case 3:
	case 4:
	case 5: return x + 1;
	case 0:
	case 1:
	case 7: return x - 1;
	default: return x;
	}
}

template <class T>
T getYatD(T y, int d) {
	switch(d % 8) {
	case 1:
	case 2:
	case 3: return y + 1;
	case 5:
	case 6:
	case 7: return y - 1;
	default: return y;
	}
}

uint64_t _tPolyId = 0;
uint64_t _tArcId = 0;

class TPoly {
public:
	uint64_t id;
	std::vector<TArc*> leftArcs;
	std::vector<TArc*> rightArcs;

	TPoly() : id(++_tPolyId) {}

	~TPoly() {
		for(TArc* a : leftArcs)
			delete a;
		for(TArc* a : rightArcs)
			delete a;
	}
}

class TArc {
public:
	uint16_t id;
	int state; // 0 - ignore; 1 - virtual; 2 solid
	std::vector<int> vertices; // Alternating x, y.

	TArc() : id(++_tArcId) {}

	void vertex(int x, int y) {
		vertices.push_back(x);
		vertices.push_back(y);
	}
}

class TArm {
public:
	int pixValue;
	int col;
	TPoly* insidePoly;
	TPoly* abovePoly;
	TPoly* leftPoly;
	TArc vertArm;
	TArc horizArm;

	TArm(int v, int c) : pixValue(v), col(c),
		insidePoly(nullptr), abovePoly(nullptr), leftPoly(nullptr) {}

	~TArm() {
		if(insidePoly) delete insidePoly;
		if(abovePoly) delete abovePoly;
		if(leftPoly) delete leftPoly;
	}
}

void Raster::polygonizeTacet(const std::string& filename, const std::string& layerName,
		const std::string& driver, uint16_t srid, uint16_t band, uint16_t threads,
		bool removeHoles, bool removeDangles, Status *status, bool *cancel) {

	std::cerr << "starting polygonize\n";

	if(!cancel)
		cancel = &s_cancel;

	const GridProps& gp = props();
	int cols = gp.cols();
	int rows = gp.rows();
	int nodata = gp.nodata();

	GridProps rp(props);
	rp.setRows(2);
	MemRaster rowBuf(rp, false);

	//std::cerr << "initializing writer\n";

	//PolyContext ctx(this, status, 0, cancel, 0, band, removeHoles, removeDangles);
	//ctx.initOutput(driver, filename, layerName, srid);

	//std::thread transT(&PolyContext::polyQueueTransfer, &ctx);
	//std::thread writeT(&PolyContext::polyWriteQueue, &ctx);

	// Organized by ID; each edge is stored as a pair of coordinates,
	// in i and i+1.
	//std::unordered_map<int, std::vector<double> > edges;

	std::vector<TArm*> arms(cols);

	writeto(rowBuf, 0, 0, cols, 2);

	int v, v0, lastV = nodata;

	// Build the arms using the first row of the buffer. From here
	// on, we use the second row and compare it to the first.
	for(int c = 0; c < cols; ++c) {
		if(c == 0 || (v = getInt(c, 0)) != lastV) {
			TArm* arm = arms[c] = new TArm(v, c);
			arm->vertArm.state = 2;
			arm->vertArm.vertex(c, 0);
			arm->vertArm.vertex(c, 1);
			arm->horizArm.state = c < cols - 1 ? 2 : 0;
			arm->horizArm.vertex(c, 0);
			arm->horizArm.vertex(c + 1, 0);
		}
		lastV = v;
	}

	// For each seed, build a path.
	for(int r = 0; r < rows - 1; ++r) {

		lastV = nodata;
		writeTo(rowBuf, 0, r, cols, 2);

		for(int c = 0; c < cols; ++c) {

			v = getInt(c, 1);
			v0 = getInt(c, 0);

			if(c == 0 || (v = getInt(c, 0)) != lastV) {

				TArm* arm;

				if(arms[c]) {
					if(arms[c].pixValue != v) {
						arm = arms[c] = new TArm(v, c);
					} else {

					}
				}

				arm->vertArm.state = 2;
				arm->vertArm.vertex(c, 0);
				arm->vertArm.vertex(c, 1);
				arm->horizArm.state = c < cols - 1 ? 2 : 0;
				arm->horizArm.vertex(c, 0);
				arm->horizArm.vertex(c + 1, 0);
			}
			lastV = v;



			int v = getInt()
			if(c == 0 || )

			for(int d = 0; d < 8; d += 2) {

				int cc = getXatD(c, d);
				int rr = getYatD(r, d);

				if(!gp.hasCell(cc, rr) || rast.getInt(cc, rr) != id) {
					 ++edges;
					 cc = c * 2 + 1;
					 rr = r * 2 + 1;
					 erast.setInt(cc, rr, id);
					 erast.setInt(getXatD(cc, d), getYatD(rr, d), -1);
					 erast.setInt(getXatD(cc, d - 1), getYatD(rr, d - 1), -1);
				}

			}
		}
	}


	ctx.readFinish();
	ctx.notifyTransfer();
	transT.join();

	ctx.transferFinish();
	ctx.notifyWrite();
	writeT.join();

	ctx.commitOutput();

}
*/



