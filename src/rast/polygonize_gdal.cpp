/*
 * polygonize_gdal.cpp
 *
 *  Created on: Aug 8, 2018
 *      Author: rob
 */

#include <sstream>
#include <iostream>

#include <gdal.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "raster.hpp"
#include "db.hpp"

using namespace geo::raster;
using namespace geo::db;

int progress(double prg, const char* message, void* data) {
	static_cast<Status*>(data)->update(prg, std::string(message));
	return 1;
}

void Grid::polygonize(const std::string& filename, const std::string& layerName,
	const std::string& driver, int srid, int band, bool removeHoles, bool removeDangles,
	const std::string& mask, int maskBand, int threads,
	bool& cancel, geo::util::Status& status) {

	// If not sqlite, use a temporary filename, otherwise use the given name.
	std::string tmpFile;
	if(driver == "SQLite") {
		tmpFile = filename;
	} else {
		tmpFile = Util::tmpFile();
	}

	// Register drivers.
	GDALAllRegister();
	OGRRegisterAll();

	// Set up a memory raster to hold the data.
	const GridProps& prop = props();

	// Load the raster driver.
	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* rdrv = dm->GetDriverByName("GTiff");
	if(!rdrv)
		g_runerr("Failed to get raster driver for polygonization.");

	std::string tmpRast = Util::tmpFile();
	GDALDataset* rast = rdrv->Create(tmpRast.c_str(), prop.cols(), prop.rows(), 1, GDT_Int32, nullptr);
	if(!rast)
		g_runerr("Failed to create raster dataset for polygonization.");
	GDALRasterBand* bnd = rast->GetRasterBand(1);
	if(!bnd) {
		GDALClose(rast);
		g_runerr("Failed to create raster band for polygonization.");
	}
	double trans[6];
	props().trans(trans);
	rast->SetGeoTransform(trans);
	bnd->SetNoDataValue(0);

	// Copy the grid data to a gdal raster.
	{
		GridProps mp(prop);
		mp.setSize(prop.cols(), 1);
		mp.setWritable(true);
		MemRaster buf(mp);
		for(int r = 0; r < prop.rows(); ++r) {
			writeTo(buf, prop.cols(), 1, 0, r, 0, 0, band, 1);
			if(CE_None != bnd->RasterIO(GF_Write, 0, r, prop.cols(), 1, buf.grid(), prop.cols(), 1, GDT_Int32, 0, 0, 0)) {
				GDALClose(rast);
				g_runerr("Failed to write to temporary raster.");
			}
		}
	}

	// Get the vector driver.
	GDALDriver* sdrv = dm->GetDriverByName("SQLite");
	if(!sdrv) {
		GDALClose(rast);
		g_runerr("Failed to get vector driver for polygonization.");
	}

	// Configure the output database.
	char** options = nullptr;
	if(driver == "SQLite")
		options = CSLSetNameValue(options, "SPATIALITE", "YES");
	GDALDataset* poly = sdrv->Create(tmpFile.c_str(), 0, 0, 0, GDT_Unknown, options);
	CSLDestroy(options);
	if(!poly) {
		GDALClose(rast);
		g_runerr("Failed to create dataset for polygonization.");
	}

	// Configure the output layer.
	OGRSpatialReference sr;
	OGRSpatialReference* srp = nullptr;
	if(srid > 0) {
		sr.importFromEPSG(srid);
		srp = &sr;
	}
	options = nullptr;
	if(driver == "SQLite")
		options = CSLSetNameValue(options, "FORMAT", "SPATIALITE");
	OGRLayer* lyr = poly->CreateLayer(layerName.c_str(), srp, wkbPolygon, options);
	CSLDestroy(options);
	if(!lyr) {
		GDALClose(rast);
		GDALClose(poly);
		g_runerr("Failed to create layer for polygonization.");
	}

	// There's only one field -- an ID.
	OGRFieldDefn field( "id", OFTInteger);
	lyr->CreateField(&field);

	// Perform polygonization.
	if(CE_None != GDALPolygonize(bnd, nullptr, lyr, 0, nullptr, progress, &status)) {
		GDALClose(rast);
		GDALClose(poly);
		g_runerr("Failed to polygonize raster.");
	}

	// Attempt to delete invalid polys. Should be a no-op but you never know.
	poly->StartTransaction();
	poly->ExecuteSQL("DELETE FROM crowns WHERE id < 1", nullptr, nullptr, nullptr);
	poly->CommitTransaction();

	GDALClose(poly);
	GDALClose(rast);

	// If necessary move or convert the file.
	if(tmpFile != filename) {
		if(driver == "SQLite") {
			Util::rename(tmpFile, filename);
		} else {
			DB db(tmpFile, layerName);
			db.convert(filename, driver);
		}
	}
}


