/*
 * grid.cpp
 *
 *  Created on: Jun 16, 2019
 *      Author: rob
 */

#include <string>
#include <algorithm>

#include "grid.hpp"

using namespace geo::grid;

Interleave geo::grid::interleaveFromString(const std::string& str) {
	std::string lower;
	std::transform(str.begin(), str.end(), lower.begin(), ::tolower);
	if(lower == "band") {
		return Interleave::BIL;
	} else {
		return Interleave::BIP;
	}
}

size_t geo::grid::detail::minValue(std::unordered_map<size_t, double>& m) {
	double min = std::numeric_limits<double>::max();
	size_t key = 0;
	for(const auto& it : m) {
		if(it.second < min) {
			min = it.second;
			key = it.first;
		}
	}
	return key;
}

void geo::grid::detail::polyWriteToFile(std::unordered_map<int, std::vector<GEOSGeometry*> >* geoms,
		std::unordered_set<int>* finalIds,
		const std::string* idField, OGRLayer* layer, GEOSContextHandle_t* gctx,
		bool removeHoles, bool removeDangles, bool* running, geo::Monitor* monitor,
		std::atomic<size_t>* poly_fid, std::mutex* poly_gmtx, std::mutex* poly_fmtx,
		std::mutex* poly_omtx, std::condition_variable* poly_cv) {


	std::vector<GEOSGeometry*> polys;
	GEOSGeometry* geom = nullptr;
	int id;

	while(!monitor->canceled() && (*running || !finalIds->empty())) {

		// Get an ID and the list of polys from the queue.
		{
			std::unique_lock<std::mutex> lk(*poly_fmtx);
			// Wait for a notification if the queue is empty.
			while(!monitor->canceled() && *running && finalIds->empty())
				poly_cv->wait(lk);
			// If the wakeup is spurious, skip.
			if(finalIds->empty())
				continue;
			// Get the ID.
			id = *(finalIds->begin());
			finalIds->erase(id);
		}
		// Get the list of polys and remove the list from the map.
		{
			std::lock_guard<std::mutex> lk(*poly_gmtx);
			// If the geoms list is still here, grab it.
			if(geoms->find(id) == geoms->end())
				continue;
			polys = geoms->at(id);
			geoms->erase(id);
		}

		if(polys.empty())
			continue;

		if(monitor->canceled()) {
			for(GEOSGeometry* p : polys)
				GEOSGeom_destroy(p);
			continue;
		}

		// Union the polys.
		geom = polys.front();
		for(size_t i = 1; i < polys.size(); ++i) {
			GEOSGeometry* tmp = GEOSUnion_r(*gctx, geom, polys[i]); //
			GEOSGeom_destroy(geom);
			GEOSGeom_destroy(polys[i]);
			geom = tmp;
		}
		polys.clear();

		if(monitor->canceled()) {
			GEOSGeom_destroy(geom);
			continue;
		}

		// If we're removing dangles, throw away all but the
		// largest single polygon. If it was originally a polygon, there are no dangles.
		size_t numGeoms;
		if(!monitor->canceled() && removeDangles && (numGeoms = GEOSGetNumGeometries(geom)) > 1) {
			size_t idx = 0;
			double a, area = 0;
			for(size_t i = 0; i < numGeoms; ++i) {
				const GEOSGeometry* p = GEOSGetGeometryN(geom, i);
				GEOSArea(p, &a);
				if(a > area) {
					area = a;
					idx = i;
				}
			}
			GEOSGeometry *g = GEOSGeom_clone(GEOSGetGeometryN(geom, idx)); // Force copy.
			GEOSGeom_destroy(geom);
			geom = g;
		}

		// If we're removing holes, extract the exterior rings of all constituent polygons.
		if(!monitor->canceled() && removeHoles) {
			std::vector<GEOSGeometry*> geoms0;
			for(int i = 0; i < GEOSGetNumGeometries(geom); ++i) {
				const GEOSGeometry* p = GEOSGetGeometryN(geom, i);
				const GEOSGeometry* l = GEOSGetExteriorRing(p);
				const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq(l);
				GEOSGeometry* r = GEOSGeom_createLinearRing(GEOSCoordSeq_clone(seq));
				GEOSGeometry* npoly = GEOSGeom_createPolygon(r, 0, 0);
				geoms0.push_back(npoly);
			}
			GEOSGeometry* g = GEOSGeom_createCollection(GEOSGeomTypes::GEOS_MULTIPOLYGON, geoms0.data(), geoms0.size()); // Do not copy -- take ownership.
			GEOSGeom_destroy(geom);
			geom = g;
		}

		// If the result is not a multi, make it one.
		if(!monitor->canceled() && GEOSGeomTypeId(geom) != GEOSGeomTypes::GEOS_MULTIPOLYGON) {
			std::vector<GEOSGeometry*> gs;
			gs.push_back(geom);
			geom = GEOSGeom_createCollection(GEOSGeomTypes::GEOS_MULTIPOLYGON, gs.data(), gs.size());
		}

		if(!geom)
			g_runerr("Null geometry.");

		if(monitor->canceled()) {
			GEOSGeom_destroy(geom);
			continue;
		}

		// Create and write the OGR geometry.
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(*gctx, (GEOSGeom) geom);
		GEOSGeom_destroy(geom);

		// Create and configure the feature. Feature owns the OGRGeometry.
		OGRFeatureDefn* fdef = layer->GetLayerDefn();
		OGRFeature* feat = OGRFeature::CreateFeature(fdef); // Creates on the OGR heap.
		feat->SetGeometryDirectly(ogeom);
		feat->SetField(idField->c_str(), (GIntBig) id);
		feat->SetFID(++*poly_fid);

		// Write to the output file.
		int err;
		{
			std::lock_guard<std::mutex> lk(*poly_omtx);
			err = layer->CreateFeature(feat);
		}

		OGRFeature::DestroyFeature(feat);

		if(OGRERR_NONE != err)
			g_runerr("Failed to add geometry.");
	}
}

/**
 * \brief Produce a rectangular polygon from the four corners.
 *
 * \param x0 The top left corner x-coordinate.
 * \param y0 The top left corner y-coordinate.
 * \param x1 The top left corner x-coordinate.
 * \param y1 The top left corner y-coordinate.
 * \param dims The number of dimensions.
 */
GEOSGeometry* geo::grid::detail::polyMakeGeom(double x0, double y0, double x1, double y1, int dims) {

	// Build the geometry.
	// TODO: Necessary to give z coord here?
	GEOSCoordSequence* seq = GEOSCoordSeq_create(5, dims);
	GEOSCoordSeq_setX(seq, 0, x0);
	GEOSCoordSeq_setY(seq, 0, y0);
	if(dims == 3)
		GEOSCoordSeq_setZ(seq, 0, 0);
	GEOSCoordSeq_setX(seq, 1, x0);
	GEOSCoordSeq_setY(seq, 1, y1);
	if(dims == 3)
		GEOSCoordSeq_setZ(seq, 0, 0);
	GEOSCoordSeq_setX(seq, 2, x1);
	GEOSCoordSeq_setY(seq, 2, y1);
	if(dims == 3)
		GEOSCoordSeq_setZ(seq, 0, 0);
	GEOSCoordSeq_setX(seq, 3, x1);
	GEOSCoordSeq_setY(seq, 3, y0);
	if(dims == 3)
		GEOSCoordSeq_setZ(seq, 0, 0);
	GEOSCoordSeq_setX(seq, 4, x0);
	GEOSCoordSeq_setY(seq, 4, y0);
	if(dims == 3)
		GEOSCoordSeq_setZ(seq, 0, 0);
	GEOSGeometry* ring = GEOSGeom_createLinearRing(seq);
	GEOSGeometry* poly = GEOSGeom_createPolygon(ring, 0, 0);
	return poly;
}

/**
 * \brief Make or open an OGR database to write the polygons to.
 *
 * \param filename The output filename.
 * \param driver The output driver.
 * \param layerName The layer name.
 * \param idField The field for the geometry ID.
 * \param sr The spatial reference object.
 * \param gType The geometry type.
 * \param ds The dataset.
 * \param layer the layer.
 */
void geo::grid::detail::polyMakeDataset(const std::string& filename, const std::string& driver, const std::string& layerName,
		const std::string& idField,
		OGRSpatialReference* sr, OGRwkbGeometryType gType,
		GDALDataset** ds, OGRLayer** layer) {

	*ds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, nullptr, nullptr, nullptr);

	std::string drvl = lowercase(driver);

	if(!*ds) {
		// Get the vector driver.
		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		if(!drv)
			g_runerr("Failed to find driver for " << driver << ".");

		// Create an output dataset for the polygons.
		char** dopts = NULL;
		if(drvl == "sqlite")
			dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

		*ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
		CSLDestroy(dopts);
		if(!*ds)
			g_runerr("Failed to create dataset " << filename << ".");

	}

	// Create the layer.
	char** lopts = NULL;
	if(drvl == "sqlite") {
		lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
	} else if(drvl == "esri shapefile") {
		lopts = CSLSetNameValue(lopts, "2GB_LIMIT", "YES");
	}
	*layer = (*ds)->CreateLayer(layerName.c_str(), sr, gType, lopts);
	CSLDestroy(lopts);

	if(!*layer)
		g_runerr("Failed to create layer " << layerName << ".");

	// There's only one field -- an ID.
	OGRFieldDefn field(idField.c_str(), OFTInteger);
	(*layer)->CreateField(&field, TRUE);

}

/**
 * \brief Return the number of bytes required to store each type.
 *
 * \param type The libgeo data type.
 * \return The number of bytes required to store each type.
 */
int geo::grid::detail::getTypeSize(DataType type) {
	switch(type) {
	case DataType::Byte: return sizeof(uint8_t);
	case DataType::Float32: return sizeof(float);
	case DataType::Float64: return sizeof(double);
	case DataType::Int16: return sizeof(int16_t);
	case DataType::Int32: return sizeof(int32_t);
	case DataType::UInt16: return sizeof(uint16_t);
	case DataType::UInt32: return sizeof(uint32_t);
	default:
		g_runerr("No size for type: " << (int) type);
	}
}

/**
 * \brief Convert the libgeo data type to a GDAL data type.
 *
 * \param The libgeo data type.
 * \return type The GDAL data type.
 */
GDALDataType geo::grid::detail::dataType2GDT(DataType type) {
	switch(type) {
	case DataType::Byte:  	return GDT_Byte;
	case DataType::UInt16: 	return GDT_UInt16;
	case DataType::UInt32:	return GDT_UInt32;
	case DataType::Int16:	return GDT_Int16;
	case DataType::Int32:	return GDT_Int32;
	case DataType::Float64:	return GDT_Float64;
	case DataType::Float32:	return GDT_Float32;
	case DataType::None:
	default:
		break;
	}
	return GDT_Unknown;
}

/**
 * \brief Convert the GDAL data type to a libgeo data type.
 *
 * \param type The GDAL data type.
 * \return The libgeo data type.
 */
DataType geo::grid::detail::gdt2DataType(GDALDataType type) {
	switch(type) {
	case GDT_Byte:	  	return DataType::Byte;
	case GDT_UInt16: 	return DataType::UInt16;
	case GDT_UInt32:	return DataType::UInt32;
	case GDT_Int16:		return DataType::Int16;
	case GDT_Int32:		return DataType::Int32;
	case GDT_Float64:	return DataType::Float64;
	case GDT_Float32:	return DataType::Float32;
	case GDT_Unknown:
	case GDT_CInt16:
	case GDT_CInt32:
	case GDT_CFloat32:
	case GDT_CFloat64:
	case GDT_TypeCount:
	default:
		break;
	}
	return DataType::None;
}


void geo::grid::detail::fixWriteBounds(int& cols, int& rows, int& srcCol, int& srcRow,
		int& dstCol, int& dstRow, int rcols, int rrows, int gcols, int grows) {

	if(cols <= 0)
		cols = gcols;
	if(rows <= 0)
		rows = grows;

	if(srcCol < 0) {
		dstCol -= srcCol;
		cols += srcCol;
		srcCol = 0;
	}
	if(srcRow < 0) {
		dstRow -= srcRow;
		rows += srcRow;
		srcRow = 0;
	}
	if(dstCol < 0) {
		cols -= dstCol;
		dstCol = 0;
	}
	if(dstRow < 0) {
		rows -= dstRow;
		dstRow = 0;
	}

	if(srcCol + cols >= rcols)
		cols = rcols - srcCol;
	if(srcRow + rows >= rrows)
		rows = rrows - srcRow;

	if(dstCol + cols > gcols)
		cols = gcols - dstCol;
	if(dstRow + rows > grows)
		rows = grows - dstRow;

}

bool geo::grid::detail::fixCoords(int& srcCol, int& srcRow, int& dstCol, int& dstRow,
		int& cols, int& rows, int srcCols, int srcRows, int dstCols, int dstRows) {

	if(cols <= 0) cols = srcCols;
	if(rows <= 0) rows = srcRows;

	if(srcCol >= srcCols || srcRow >= srcRows || srcCol + cols < 0 || srcRow + rows < 0) {
		g_warn("Col/row out of range." << srcCol << ", " << srcRow << ", " << srcCols << ", " << srcRows << ", " << cols << ", " << rows);
		return false;
	}
	if(srcCol < 0) {
		cols += srcCol;
		srcCol = 0;
	}
	if(srcRow < 0) {
		rows += srcRow;
		srcRow = 0;
	}
	if(srcCol + cols > srcCols) {
		cols = srcCols - srcCol;
	}
	if(srcRow + rows > srcRows) {
		rows = srcRows - srcRow;
	}
	if(dstCol < 0) {
		cols += dstCol;
		dstCol = 0;
	}
	if(dstRow < 0) {
		rows += dstRow;
		dstRow = 0;
	}
	if(dstCol + cols > dstCols) {
		cols = dstCols - dstCol;
	}
	if(dstRow + rows > dstRows) {
		rows = dstRows - dstRow;
	}
	if(cols <= 0 || rows <= 0) {
		g_warn("Zero copy area.")
		return false;
	}
	return true;
}


