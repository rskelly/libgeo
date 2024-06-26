/*
 * grid.cpp
 *
 *  Created on: Jun 16, 2019
 *      Author: rob
 */

#include <string>
#include <algorithm>
#include <limits>

#include "grid.hpp"

bool gdalInited = false;

void dijital::grid::init() {
	if(!gdalInited) {
		GDALAllRegister();
		gdalInited = true;
	}
}

Interleave dijital::grid::interleaveFromString(const std::string& str) {
	std::string lower;
	std::transform(str.begin(), str.end(), lower.begin(), ::tolower);
	if(lower == "band") {
		return Interleave::BIL;
	} else {
		return Interleave::BIP;
	}
}

size_t dijital::grid::detail::minValue(std::unordered_map<size_t, double>& m) {
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


void dijital::grid::detail::polyWriteToDB(PolygonContext* pc) {

	GEOSGeometry* geom = nullptr;
	int id;

	while(!Monitor::get().canceled() && (pc->writeRunning || !pc->geoms.empty())) {

		// Get an ID and the list of polys from the queue.
		{
			// Wait for if the queue is empty.
			std::unique_lock<std::mutex> lk(pc->gmtx);
			while (!Monitor::get().canceled() && pc->writeRunning && pc->geoms.empty())
				pc->gcv.wait(lk);
			// If the wakeup is spurious, skip.
			if(pc->geoms.empty())
				continue;
			// Get the ID and geom.
			id = pc->geoms.front().first;
			geom = pc->geoms.front().second;
			pc->geoms.pop_front();
		}

		if(Monitor::get().canceled()) {
			GEOSGeom_destroy_r(pc->gctx, geom);
			continue;
		}

		// Create and write the OGR geometry.
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(pc->gctx, (GEOSGeom) geom);

		// Create and configure the feature. Feature owns the OGRGeometry.
		OGRFeatureDefn* fdef = pc->layer->GetLayerDefn();
		OGRFeature* feat = OGRFeature::CreateFeature(fdef); // Creates on the OGR heap.
		if(pc->geomField.empty()) {
			feat->SetGeometryDirectly(ogeom);
		} else {
			feat->SetGeomFieldDirectly(feat->GetGeomFieldIndex(pc->geomField.c_str()), ogeom);
		}
		feat->SetFID(-1);
		feat->SetField(pc->idField.c_str(), (GIntBig) id);
		for(PolygonValue v : pc->fieldValues) {
			switch(v.type) {
			case OFTInteger:
				feat->SetField(v.name.c_str(), v.ivalue);
				break;
			case OFTReal:
				feat->SetField(v.name.c_str(), v.dvalue);
				break;
			case OFTWideString:
			case OFTString:
				feat->SetField(v.name.c_str(), v.svalue.c_str());
				break;
			default:
				g_runerr("Unknown field type: " << v.type);
			}
		}

		// Write to the output file.
		int err = OGRERR_NONE;
		{
			std::lock_guard<std::mutex> lk(pc->lmtx);
			err = pc->layer->CreateFeature(feat);
		}

		OGRFeature::DestroyFeature(feat);
		GEOSGeom_destroy_r(pc->gctx, geom);

		if(OGRERR_NONE != err)
			g_runerr("Failed to add geometry.");
	}
}

void dijital::grid::detail::printGEOSGeom(GEOSGeometry* geom, GEOSContextHandle_t gctx) {
	GEOSWKTWriter* wtr = GEOSWKTWriter_create_r(gctx);
	char* a = GEOSWKTWriter_write_r(gctx, wtr, geom);
	g_warn(a);
	free(a);
	GEOSWKTWriter_destroy_r(gctx, wtr);
}

GEOSGeometry* dijital::grid::detail::polyMakeGeom(GEOSContextHandle_t gctx, double x0, double y0, double x1, double y1, double eps, int dims) {

	// Build the geometry.
	GEOSCoordSequence* seq = GEOSCoordSeq_create_r(gctx, 5, dims);
	GEOSCoordSeq_setX_r(gctx, seq, 0, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 0, y0);
	GEOSCoordSeq_setX_r(gctx, seq, 1, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 1, y1);
	GEOSCoordSeq_setX_r(gctx, seq, 2, x1);
	GEOSCoordSeq_setY_r(gctx, seq, 2, y1);
	GEOSCoordSeq_setX_r(gctx, seq, 3, x1);
	GEOSCoordSeq_setY_r(gctx, seq, 3, y0);
	GEOSCoordSeq_setX_r(gctx, seq, 4, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 4, y0);
	GEOSGeometry* ring = GEOSGeom_createLinearRing_r(gctx, seq);
	GEOSGeometry* poly = GEOSGeom_createPolygon_r(gctx, ring, 0, 0);
	GEOSGeometry* prec = GEOSGeom_setPrecision_r(gctx, poly, eps, 0);
	GEOSGeom_destroy_r(gctx, poly);
	return prec;
}

void dijital::grid::detail::polyMakeDataset(const std::string& filename, const std::string& driver,
		const std::string& layerName, const std::string& idField,
		const std::vector<PolygonValue>& fields,
		OGRSpatialReference* sr, OGRwkbGeometryType gType,
		GDALDataset** ds, OGRLayer** layer) {

	rem(filename);

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

	// Make fields if necessary.
	for(const PolygonValue& v : fields) {
		if(v.name != idField) {
			OGRFieldDefn field(v.name.c_str(), v.type);
			(*layer)->CreateField(&field, TRUE);
		}
	}

}

int dijital::grid::detail::getTypeSize(DataType type) {
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

GDALDataType dijital::grid::detail::dataType2GDT(DataType type) {
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

DataType dijital::grid::detail::gdt2DataType(GDALDataType type) {
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

void dijital::grid::detail::fixWriteBounds(int& cols, int& rows, int& srcCol, int& srcRow,
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

bool dijital::grid::detail::fixCoords(int& srcCol, int& srcRow, int& dstCol, int& dstRow,
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

int dijital::grid::detail::gdalProgress(double dfComplete, const char* pszMessage, void*pProgressArg) {
	struct gdalprg* prg = (struct gdalprg*) pProgressArg;
	int p0 = (int) (dfComplete * 100);
	if(prg->p != p0 && p0 % 10 == 0) {
		Monitor::get().status((float) dfComplete, std::string(pszMessage));
		prg->p = p0;
	}
	return Monitor::get().canceled() ? 0 : 1;
};
