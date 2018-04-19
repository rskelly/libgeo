#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <liblas/liblas.hpp>
#include <proj_api.h>
#include <gdal_priv.h>
#include <gdal.h>
#include <cpl_conv.h> // for CPLMalloc()

#include "geod/NAD83VG.hpp"

#define PI 3.14159265358979323846
#define SEC2RAD 4.84813681 / 1000000000.0

#define LAS2CSRS_DATA "LAS2CSRS_DATA"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;


// Convert from arcsec to radians.
double _sec2rad(double x) {
	return x * SEC2RAD;
}

// Convert from radians to degrees.
double _deg(double x) {
	return x * 180.0 / PI;
}

double _rad(double x) {
	return x * PI / 180.0;
}

double _sq(double x) {
	return x*x;
}

/**
 * Convert projected distance in mm to latlon in radians.
 * dx, dy  -- Distance in m.
 * lat     -- The latitude at which distances are computed
 * h       -- Ellipsoidal height.
 * a       -- Semi-major axis
 * e2      -- Eccentricity^2
 * count   -- Number of points.
 * dlat    -- The distance in rad (out).
 * dlon    -- The distance in rad (out).
 */
void _shift2latlon(std::vector<double>& dx, std::vector<double>& dy,
		std::vector<double>& lat, std::vector<double>& h, double a, double e2,
		std::vector<double>& dlat, std::vector<double>& dlon) {
	double r, m, n;
	for(int i = 0; i < dx.size(); ++i) {
		m = a * (1.0 - e2) / std::pow((1.0 - e2 * std::pow(std::sin(lat[i]), 2.0)), 3.0 / 2.0); 	// Meridional radius of curvature.
		n = a / std::pow((1.0 - e2 * std::pow(sin(lat[i]), 2.0)), 1.0 / 2.0); 				// Parallel radius of curvature.
		r = n * std::cos(lat[i]); 												// Radius of parallel.
		dlon[i] = dx[i] / (r + h[i]);
		dlat[i] = dy[i] / (m + h[i]);
	}
}

/**
 * Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
 */
class Transformer {
private:
	geo::geod::NAD83VG m_shift;

	// Source reference frame.
	std::string ffrom;
	// From and to epochs.
	double efrom;
	double eto;
	// Projection objects.
	projPJ p1;
	projPJ p2;
	projPJ p3;
	projPJ p4;
	// SRIDs of the from and to CRSes.
	int fsrid;
	int tsrid;

	// Transform parameters; loaded from the itrf file.
	double tx, ty, tz, dtx, dty, dtz;		// Shifts, rates.
	double rx, ry, rz, drx, dry, drz;		// Rotations, rates.
	double epoch;							// ITRF Transform epoch.
	double dt; 								// Time delta
	double d, dd; 							// Scale, scale rate.

	/**
		Transform the coordinate using the procedure listed in Craymer (2006).

		x, y, z -- 	The coordinate arrays.
		count 	-- 	The number of points.
	*/
	void epochTransform(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double dt) {
		double a0 = tx + dtx * dt;
		double a1 = ty + dty * dt;
		double a2 = tz + dtz * dt;
		double bsx = 1.0 + d + dd * dt;
		double b01 = -_sec2rad(rz + drz * dt);
		double b02 = _sec2rad(ry + dry * dt);
		double b10 = _sec2rad(rz + drz * dt);
		double b12 = -_sec2rad(rx + drx * dt);
		double b20 = -_sec2rad(ry + dry * dt);
		double b21 = _sec2rad(rx + drx * dt);
		for(int i = 0; i < x.size(); ++i) {
			x[i] = a0 + bsx * x[i] + b01 * y[i] + b02 * z[i];
			y[i] = a1 + b10 * x[i] + bsx * y[i] + b12 * z[i];
			z[i] = a2 + b20 * x[i] + b21 * y[i] + bsx * z[i];
		}
	}

	/**
	 * Initialize the proj4 projection objects.
	 */
	void initProjections() {
		// Initialize projections.
		if(!fsrid || !tsrid)
			throw "SRIDs are not set.";
		char str[128];
		sprintf(str, "+init=EPSG:%u", fsrid);
		p1 = pj_init_plus(str);
		if(!p1) {
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		sprintf(str, "+init=EPSG:%u", tsrid);
		p3 = pj_init_plus(str);
		if(!p3){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p2 = pj_init_plus("+init=EPSG:4978");
		if(!p2){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
		p4 = pj_init_plus("+init=EPSG:4326");
		if(!p4){
			const char *err = pj_strerrno(pj_errno);
			throw err;
		}
	}

	/**
	 * Load the transformation database.
	 */
	void loadHelmert() {

		char path[PATH_MAX];
		sprintf(path, "%s%s", std::getenv(LAS2CSRS_DATA), "/itrf.csv");

		char ffrom[64], fto[64];
		bool found = false;
		float epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd;
		FILE *f = fopen(path, "r");
		if(f == NULL)
			throw "ITRF database file not found.";
		while(fscanf(f, " %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
						ffrom, fto, &epoch, &tx, &ty, &tz, &rx, &ry, &rz, &d,
						&dtx, &dty, &dtz, &drx, &dry, &drz, &dd) > 0) {
			std::string _ffrom(ffrom);
			if(_ffrom == this->ffrom) {
				this->epoch = epoch;
				this->d = d / 1000000000.0; 		// Listed as ppb.
				this->dd = dd / 1000000000.0;
				this->tx = tx;
				this->ty = ty;
				this->tz = tz;
				this->rx = rx;
				this->ry = ry;
				this->rz = rz;
				this->dtx = dtx;
				this->dty = dty;
				this->dtz = dtz;
				this->drx = drx;
				this->dry = dry;
				this->drz = drz;
				found = true;
				break;
			}
		}
		fclose(f);

		if(!found)
			throw "Failed to find a transformation matching the parameters.";
	}
public:

	/**
	 * Prepares the Helmert transformation parameters for the given transformation.
	 * These will be used in later method calls.
	 * The constructor will load an process the grid shift file and the transformation database.
	 * ffrom -- The name of the reference frame, e.g. 'itrf90'
	 * efrom -- The epoch of data collection (decimal years), e.g. 1994.2
	 * eto   -- The target epoch. One might select 1997.0 for BC or 2002.0 for Albera (etc.)
	 * fsrid -- The SRID (EPSG code) of the source.
	 * tsrid -- The SRID of the destination. This is the code for the UTM zone under
	 *          NAD83(CSRS), e.g. 2956 for UTM12N.
	 */
	Transformer(std::string &ffrom, float efrom, float eto, int fsrid, int tsrid) {
		this->efrom = efrom;
		this->eto = eto;
		this->fsrid = fsrid;
		this->tsrid = tsrid;
		this->ffrom.assign(ffrom);

		initProjections();
		loadHelmert();
	}

	/**
	 * Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	 * x, y, z -- Coordinate arrays.
	 * count   -- The number of coordinates.
	 * bounds  -- The bounds of the transformed coordinate list.
	 */
	void transformPoints(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double bounds[6]) {

		// Project to Cartesian 3D. (c)
		pj_transform(p1, p2, x.size(), 1, (double*) x.data(), (double*) y.data(), (double*) z.data());

		// Transform to csrs using Helmert (etc.) params. (c)
		epochTransform(x, y, z, this->efrom - 1997.0);

		// Only use the grid shift if the epoch changes.
		if(efrom != eto) {

			// Copy the coordinate arrays for transformation.
			std::vector<double> x0(x.size());
			std::vector<double> y0(y.size());
			std::vector<double> z0(z.size());
			std::copy(x.begin(), x.end(), std::back_inserter(x0));
			std::copy(y.begin(), y.end(), std::back_inserter(y0));
			std::copy(z.begin(), z.end(), std::back_inserter(z0));

			// Initalize shift arrays.
			std::vector<double> dx(x.size());
			std::vector<double> dy(y.size());
			std::vector<double> dz(z.size());
			
			// Transform to latlon. (b)
			pj_transform(p2, p4, x.size(), 1, (double*) x0.data(), (double*) y0.data(), (double*) z0.data());
			
			for(size_t i = 0; i < x.size(); ++i)
				m_shift.getVelocity(x[i], y[i], dx[i], dy[i], dz[i]);
			
			// Transform mm shifts to latlon
			std::vector<double> dlat(x.size());
			std::vector<double> dlon(y.size());

			// Get projection's spheroid props.
			double a, e2;
			pj_get_spheroid_defn(p4, &a, &e2);

			// Change the shift distance in projected coords to latlon.
			// This avoids the scale distortion associated with projection.
			_shift2latlon(dx, dy, y0, z0, a, e2, dlat, dlon);

			double dt = this->eto - this->efrom;
			// Apply shifts to latlon coords.
			for(size_t i = 0; i < x0.size(); ++i) {
				x0[i] += dlon[i] * dt;
				y0[i] += dlat[i] * dt;
				z0[i] += dz[i] * dt;
			}
			
			// Transform latlon to target proj
			pj_transform(p4, p3, x0.size(), 1, (double*) x0.data(), (double*) y0.data(), (double*) z0.data());
			
			// Assign the shifted coords to the output arrays.
			std::copy(x0.begin(), x0.end(), std::back_inserter(x));
			std::copy(y0.begin(), y0.end(), std::back_inserter(y));
			std::copy(z0.begin(), z0.end(), std::back_inserter(z));

		} else {

			// Reproject to the dest coordinates
			pj_transform(p2, p3, x.size(), 1, (double*) x.data(), (double*) y.data(), (double*) z.data());

		}

		// Expand the bounds for the new header.
		for(size_t i = 0; i < x.size(); ++i) {
			if(x[i] < bounds[0]) bounds[0] = x[i];
			if(x[i] > bounds[1]) bounds[1] = x[i];
			if(y[i] < bounds[2]) bounds[2] = y[i];
			if(y[i] > bounds[3]) bounds[3] = y[i];
			if(z[i] < bounds[4]) bounds[4] = z[i];
			if(z[i] > bounds[5]) bounds[5] = z[i];
		}			
	}

	/**
	 * Transform the given LAS file(s) from the given reference frame to the one
	 * configured in this Transformer. 
	 * srcfile  -- The folder containing las files, or the path to a single las file.
	 * dstdir   -- The destination folder.
	 */
	int transformLas(std::string &srcfile, std::string &dstdir, bool overwrite) {
		
		const fs::path src(srcfile);
		const fs::path dst(dstdir);

		// Check the files for sanity.
		if(src == dst)
			throw std::string("Destination and source are the same: ") + dst.string() + std::string("==") + src.string();
		//if(!fs::is_directory(dst))
		//	throw std::string("Destination is a file: ") + dst.string();
		if(!fs::exists(dst)) {
			if(!fs::create_directory(dst))
				throw std::string("Failed to create output directory: ") + dst.string();
		}

		// Get the list of las files.
		// TODO: Ignore completed files.
		std::list<fs::path> files;
		if(fs::is_regular_file(src)) {
			if(overwrite || !fs::exists(dst / src.leaf()))
				files.push_back(src);
		} else {
			std::string ext(".las");
			fs::directory_iterator end;
			fs::directory_iterator di(src);
			for(; di != end; ++di) {
				std::string p(di->path().string());
				alg::to_lower(p);
				if((overwrite || !fs::exists(dst / di->path().leaf())) && alg::ends_with(p, ext))
					files.push_back(di->path());
			}
		}
		if(files.size() == 0)
			throw "No matching files found.";

		std::cout << "Processing " << files.size() << " files." << std::endl;

		// Start
		las::WriterFactory wf;
		las::ReaderFactory rf;

		// The overall bounds: min x, max x, min y, max y, min z, max z
		double bounds[] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };

		for(std::list<fs::path>::iterator it = files.begin(); it != files.end(); ++it) {

			const char * filename = it->c_str();

			std::cout << "Processing file " << filename << std::endl;

			// Open the source file.
			std::ifstream in(filename, std::ios::in | std::ios::binary);
			las::Reader r = rf.CreateWithStream(in);
			las::Header h = r.GetHeader();

			// Open the destination file.
			fs::path outfile(dst / it->leaf());
			las::Header dsth(h);

			int count = h.GetPointRecordsCount();

			std::vector<double> x(count);
			std::vector<double> y(count);
			std::vector<double> z(count);

			// Iterate over the points.
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				x[i] = pt.GetX();
				y[i] = pt.GetY();
				z[i] = pt.GetZ();
			}

			// Transform the points in-place.
			transformPoints(x, y, z, bounds);

			// Set bounds.
			dsth.SetMin(bounds[0], bounds[2], bounds[4]);
			dsth.SetMax(bounds[1], bounds[3], bounds[5]);
			std::ofstream out(outfile.c_str(), std::ios::out | std::ios::binary);
			las::Writer w(out, dsth);

			// Iterate over the points.
			r.Reset();
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				pt.SetX(x[i]);
				pt.SetY(y[i]);
				pt.SetZ(z[i]);
				// Write to the output
				w.WritePoint(pt);
			}

			// Close files.
			in.close();
			out.close();

		}

		return 0;

	}

	~Transformer() {
		pj_free(p1);
		pj_free(p2);
		pj_free(p3);
		pj_free(p4);
	}


};

	
void usage() {
	std::cout << "Usage: las2csrs [options] <src file or dir> <dst dir> <src ref frame> <src epoch> <dst epoch> <srd srid> <dst srid>" << std::endl;
	std::cout << " -o     Overwrite existing files. Defaults to false." << std::endl;
}

int test(int argc, char **argv) {

	if(argc < 10) {
		std::cerr << "Too few arguments." << std::endl;
		return 1;
	}

	int i = 1;
	std::string ffrom(argv[++i]);
	double efrom = atof(argv[++i]);
	double eto = atof(argv[++i]);
	int fsrid = atoi(argv[++i]);
	int tsrid = atoi(argv[++i]);
	double x = atof(argv[++i]);
	double y = atof(argv[++i]);
	double z = atof(argv[++i]);

	Transformer trans(ffrom, efrom, eto, fsrid, tsrid);

	std::vector<double> xx; xx.push_back(x);
	std::vector<double> yy; yy.push_back(y);
	std::vector<double> zz; zz.push_back(z);
	double bounds[6];

	trans.transformPoints(xx, yy, zz, bounds);

	std::cout << std::setprecision(12) << xx[0] << " " << yy[0] << " " << zz[0] << std::endl;

	return 0;
}

int main(int argc, char **argv) {

	std::string comp("test");
	if(argc >= 2 && comp.compare(argv[1]) == 0) {
		return test(argc, argv);
	}

	if(argc < 8) {
		usage();
		return 1;
	}

	bool overwrite = false;
	try {

		int i = 0;
		for(; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg.c_str()[0] != '-') {
				break;
			} else if(arg == "-o") {
				overwrite = true;
			}
		}

		std::string srcfile(argv[++i]);
		std::string dstdir(argv[++i]);
		std::string ffrom(argv[++i]);
		double efrom = atof(argv[++i]);
		double eto = atof(argv[++i]);
		int fsrid = atoi(argv[++i]);
		int tsrid = atoi(argv[++i]);

		if(!std::getenv(LAS2CSRS_DATA))
			setenv(LAS2CSRS_DATA, "..", 1);

		Transformer trans(ffrom, efrom, eto, fsrid, tsrid);

		trans.transformLas(srcfile, dstdir, overwrite);

	} catch(char const *err) {
		std::cerr << err << std::endl;
		return 1;
	} catch(const std::string &err) {
		std::cerr << err << std::endl;
	}

	return 0;

}

