/**
 */

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <liblas/liblas.hpp>

#include "util.hpp"

#define RASTER 1
#define VECTOR 2
#define POINT_CLOUD 3

std::vector<std::string> charToVector(char** arr) {
	std::vector<std::string> vec;
	char** tmp = arr;
	while(tmp != NULL) {
		vec.push_back(std::string(*tmp));
		++tmp;
	}
	CSLDestroy(arr);
	return vec;
}

class Dataset {
public:
	int type;
	double bounds[4];
	std::string projection;
	std::vector<std::string> files;

	Dataset(int type) : 
		type(type),
		bounds{99999999.,99999999.,-99999999.,-99999999.} {
	}

};

class RasterBand {
public:
	std::string dataType;
	int cols;
	int rows;
	double nodata;
	bool hasNodata;
	double min;
	double max;
	double mean;
	double stdDev;
	std::unordered_map<std::string, std::string> metadata;

	RasterBand(GDALRasterBand* band) {
		
		dataType = band->GetRasterDataType();
		cols = band->GetXSize();
		rows = band->GetYSize();
		nodata = band->GetNoDataValue((int*) &hasNodata);
		band->GetStatistics(0, 1, &min, &max, &mean, &stdDev);

		if(!hasNodata)
			nodata = 0;
	}
};

class Raster : public Dataset {
public:
	int cols;
	int rows;
	double resX;
	double resY;
	std::vector<RasterBand> bands;

	Raster(GDALDataset* ds) : Dataset(RASTER) {
		
		projection = std::string(ds->GetProjectionRef());
		files = charToVector(ds->GetFileList());

		cols = ds->GetRasterXSize();
		rows = ds->GetRasterYSize();

		double trans[6];
		ds->GetGeoTransform(trans);

		resX = trans[1];
		resY = trans[5];

		bounds[0] = trans[0];
		bounds[1] = trans[3];
		bounds[2] = bounds[0] + cols * resX;
		bounds[3] = bounds[1] + rows * resY;
		if(bounds[0] > bounds[2]) {
			double tmp = bounds[0];
			bounds[0] = bounds[2];
			bounds[2] = tmp;
		}
		if(bounds[1] > bounds[3]) {
			double tmp = bounds[1];
			bounds[1] = bounds[3];
			bounds[3] = tmp;
		}

		for(int i = 0; i < ds->GetRasterCount(); ++i)
			bands.push_back(RasterBand(ds->GetRasterBand(i)));
	}
};

class VectorLayer {
public:
	int geomType;
	std::string name;
	double bounds[4];

	VectorLayer(OGRLayer* layer) :
		geomType(layer->GetGeomType()),
		name(layer->GetName()) {

		OGREnvelope* env = nullptr;
		if(OGRERR_FAILURE != layer->GetExtent(env, true) && env) {
			bounds[0] = env->MinX;
			bounds[1] = env->MinY;
			bounds[2] = env->MaxX;
			bounds[3] = env->MaxY;
		}
	}
};

class Vector : public Dataset {
public:
	std::vector<VectorLayer> layers;

	Vector(GDALDataset* ds) : Dataset(VECTOR) {

		projection = std::string(ds->GetProjectionRef());
		files = charToVector(ds->GetFileList());

		for(int i = 0; i < ds->GetLayerCount(); ++i) {
			layers.push_back(VectorLayer(ds->GetLayer(i)));
			VectorLayer& lyr = layers[i];
			if(lyr.bounds[0] < bounds[0]) bounds[0] = lyr.bounds[0];
			if(lyr.bounds[1] < bounds[1]) bounds[1] = lyr.bounds[1];
			if(lyr.bounds[2] > bounds[2]) bounds[2] = lyr.bounds[2];
			if(lyr.bounds[3] > bounds[3]) bounds[3] = lyr.bounds[3];
		}
	}
};


class PointCloud : public Dataset {
public:

	PointCloud(liblas::Reader* reader) : Dataset(POINT_CLOUD) {
		const liblas::Header& hdr = reader->GetHeader();
		bounds[0] = hdr.GetMinX();
		bounds[1] = hdr.GetMinY();
		bounds[2] = hdr.GetMaxX();
		bounds[3] = hdr.GetMaxY();
	}

};

bool tryGDAL(const std::string& filename, std::unique_ptr<Dataset>& ds) {

	GDALDataset* gds;
	GDALDriver* gdrv;

	try {
		
		if(!(gds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly)))
			g_runerr("Failed to open GDAL dataset.");

		if(!(gdrv = gds->GetDriver()))
			g_runerr("Could not retrieve driver.");

		if(gds->GetRasterCount()) {
			// This is a raster.
			
			ds.reset(new Raster(gds));

		} else if(gds->GetLayerCount()) {
			// This is a vector.

			ds.reset(new Vector(gds));

		} else {
			g_runerr("Dataset is not a vector or raster. Something is amiss (missing driver?)");
		}

		GDALClose(gds);
		return true;

	} catch(std::exception& ex) {
		g_debug("It's not a valid GDAL dataset: " << ex.what());
		GDALClose(gds);
		return false;
	}

}

bool tryLAS(const std::string& filename, std::unique_ptr<Dataset>& ds) {

	try {
		liblas::ReaderFactory rfact;
		std::ifstream str(filename, std::ios::in | std::ios::binary);
		liblas::Reader reader = rfact.CreateWithStream(str);
		
		ds.reset(new PointCloud(&reader));
		ds->files.push_back(filename);

		return true;

	} catch(const std::exception& ex) {
		g_debug("It's not a valid LAS dataset: " << ex.what());
		return false;
	}

}

int handleFile(const std::string& filename) {

	std::unique_ptr<Dataset> ds;

	if(!tryGDAL(filename, ds)) {
		if(!tryLAS(filename, ds)) {
			g_runerr("Unknown file type: " << filename);
		}
	}

	return 0;

}

void usage() {
	std::cerr << "Usage: info <filename>\n";
}

int main(int argc, char** argv) {
	
	if(argc < 2) {
		usage();
		return 1;
	}

	GDALAllRegister();

	std::string filename(argv[1]);

	try {
		return handleFile(filename);
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
		usage();
		return 1;
	}

	return 0;
}