/**
 */

#include <string>
#include <vector>

#include <gdal.h>
#include <ogrsf_frmts.h>

#include <liblas/liblas.hpp>

#include <json/value.h>
#include <json/writer.h>

#include "util.hpp"

#define RASTER 1
#define VECTOR 2
#define POINT_CLOUD 3

Json::Value asGeom(double* bounds) {
	Json::Value node(Json::objectValue);
	Json::Value geom(Json::objectValue);
	Json::Value coords(Json::arrayValue);
	{
		Json::Value coord(Json::arrayValue);
		coord.append(bounds[0]);
		coord.append(bounds[1]);
		coords.append(coord);
	}
	{
		Json::Value coord(Json::arrayValue);
		coord.append(bounds[2]);
		coord.append(bounds[1]);
		coords.append(coord);
	}
	{
		Json::Value coord(Json::arrayValue);
		coord.append(bounds[2]);
		coord.append(bounds[3]);
		coords.append(coord);
	}
	{
		Json::Value coord(Json::arrayValue);
		coord.append(bounds[0]);
		coord.append(bounds[3]);
		coords.append(coord);
	}
	{
		Json::Value coord(Json::arrayValue);
		coord.append(bounds[0]);
		coord.append(bounds[1]);
		coords.append(coord);
	}
	geom["type"] = "polygon";
	geom["coordinates"] = coords;
	node["geometry"] = geom;
	return node;

}
std::vector<std::string> charToVector(char** arr) {
	std::vector<std::string> vec;
	char** tmp = arr;
	while(*tmp != NULL) {
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

	virtual Json::Value asJSON() {
		Json::Value node(Json::objectValue);
		node["type"] = type;
		node["projection"] = projection;
		node["bounds"] = asGeom(bounds);
		Json::Value jfiles(Json::arrayValue);
		for(const std::string& file : files)
			jfiles.append(file);
		node["files"] = jfiles;
		return node;
	}

	virtual ~Dataset() {}

};

class RasterBand {
public:
	int dataType;
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

	Json::Value asJSON() {
		Json::Value node(Json::objectValue);
		node["dataType"] = dataType;
		node["cols"] = cols;
		node["rows"] = rows;
		node["hasNodata"] = hasNodata;
		if(hasNodata)
			node["nodata"] = nodata;
		Json::Value stats(Json::objectValue);
		stats["min"] = min;
		stats["max"] = max;
		stats["mean"] = mean;
		stats["stddev"] = stdDev;
		node["stats"] = stats;
		return node;
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
			bands.push_back(RasterBand(ds->GetRasterBand(i + 1)));
	}

	Json::Value asJSON() {
		Json::Value node = Dataset::asJSON();
		node["cols"] = cols;
		node["rows"] = rows;
		node["resolutionX"] = resX;
		node["resolutionY"] = resY;
		node["bands"] = Json::Value(Json::arrayValue);
		for(RasterBand& band : bands)
			node["bands"].append(band.asJSON());
		return node;
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

		OGREnvelope env;
		if(OGRERR_FAILURE != layer->GetExtent(&env, true)) {
			bounds[0] = env.MinX;
			bounds[1] = env.MinY;
			bounds[2] = env.MaxX;
			bounds[3] = env.MaxY;
		}
	}

	Json::Value asJSON() {
		Json::Value node(Json::objectValue);
		node["geomType"] = geomType;
		node["name"] = name;
		node["bounds"] = asGeom(bounds);
		return node;
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

	Json::Value asJSON() {
		Json::Value node = Dataset::asJSON();
		node["layers"] = Json::Value(Json::arrayValue);
		for(VectorLayer& layer : layers)
			node["layers"].append(layer.asJSON());
		return node;
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
		
		if(!(gds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL)))
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
			// Zippo.
		}
	}

	Json::Value node(Json::objectValue);

	if(ds.get()) {
		node["dataset"] = ds->asJSON();
	} else {
		node = Json::Value(Json::objectValue);
		node["error"] = "Failed to understand file.";
	}

	Json::StreamWriterBuilder wb;
	std::cout << Json::writeString(wb, node);

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
