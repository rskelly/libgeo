#include <unordered_set>

#include "db.hpp"
#include "geo.hpp"
#include "util.hpp"

using namespace geo::db;
using namespace geo::util;

namespace geo {

	namespace db {

		namespace util {

			GeomType geomType(OGRwkbGeometryType type) {
				switch(type) {
				case wkbPoint: return GeomType::GTPoint;
				case wkbPoint25D: return GeomType::GTPointZ;
				case wkbLineString: return GeomType::GTLine;
				case wkbPolygon: return GeomType::GTPolygon;
				case wkbMultiPoint: return GeomType::GTMultiPoint;
				case wkbMultiLineString: return GeomType::GTMultiLine;
				case wkbMultiPolygon: return GeomType::GTMultiPolygon;
				case wkbMultiPolygon25D: return GeomType::GTMultiPolygonZ;
				default: return GeomType::GTUnknown;
				}
			}

			OGRwkbGeometryType geomType(GeomType type) {
				switch(type) {
				case GeomType::GTPoint: return wkbPoint;
				case GeomType::GTPointZ: return wkbPoint25D;
				case GeomType::GTLine: return wkbLineString;
				case GeomType::GTPolygon: return wkbPolygon;
				case GeomType::GTMultiPoint: return wkbMultiPoint;
				case GeomType::GTMultiLine: return wkbMultiLineString;
				case GeomType::GTMultiPolygon: return wkbMultiPolygon;
				case GeomType::GTMultiPolygonZ: return wkbMultiPolygon25D;
				default: return wkbUnknown;
				}
			}

			FieldType fieldType(OGRFieldType type) {
				switch(type) {
				case OFTInteger:
				case OFTInteger64: return FieldType::FTInt;
				case OFTString: return FieldType::FTString;
				case OFTReal: return FieldType::FTDouble;
				case OFTBinary: return FieldType::FTBlob;
				default: return FieldType::FTUnknown;
				}
			}

			OGRFieldType fieldType(FieldType type) {
				switch(type) {
				case FieldType::FTInt: return OFTInteger64;
				case FieldType::FTString: return OFTString;
				case FieldType::FTDouble: return OFTReal;
				case FieldType::FTBlob: return OFTBinary;
				case FieldType::FTUnknown:
				default:
					g_argerr("Unknown or unimplemented type: " << type);
				}
			}

			bool isRast(GDALDriver *drv) {
				const char* cc = drv->GetMetadataItem(GDAL_DCAP_RASTER);
				return cc != nullptr && std::strncmp(cc, "YES", 3) == 0;
			}

		} // util
	} // db
} // geo

using namespace geo::db::util;

DB::DB(const std::string& file, const std::string& layer, const std::string& driver,
		const std::unordered_map<std::string, FieldType>& fields,
		GeomType type, const std::string& projection, bool replace) :
    m_type(type),
	m_srid(0),
    m_projection(projection),
    m_file(file),
	m_layerName(layer),
	m_driver(driver),
	m_fieldTypes(fields),
	m_ds(nullptr),
	m_layer(nullptr),
	m_fdef(nullptr) {

	GDALAllRegister();

	// If the driver was not given, try to discover it from an existing file, otherwise fail.
	if (m_driver.empty()) {
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(m_file.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr));
		if (!ds)
			g_runerr("Driver not given and no existing file to read from.");
		const char* drv = ds->GetDriverName();
		GDALClose(ds);
		if(!drv)
			g_runerr("Driver not given and no existing file to read from.");
		m_driver = drv;
		g_warn("Using driver " << drv << " from existing file.");
	}

	// Check layer name.
	if (m_layerName.empty()) {
		m_layerName = "data";
		g_warn("Using layer name data as none given.");
	}

	// If replace and file exists, delete the existing file.
    if(replace && isfile(file))
        rem(file);

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(m_driver.c_str());
    if(!drv)
        g_runerr("Driver not found for " << m_file << " (" << m_driver << ")");

	// If the file is sqlite, use the spatialite driver.
	char** dopts = nullptr;
	if(m_driver == "SQLite")
		dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

	m_ds = drv->Create(m_file.c_str(), 0, 0, 0, GDT_Unknown, dopts);

	if(dopts)
		CSLDestroy(dopts);

	if(!m_ds)
        g_runerr("Failed to create data set for " << m_file);

	OGRSpatialReference* sr = nullptr;
	if(!m_projection.empty())
		sr = new OGRSpatialReference(m_projection.c_str());

	dopts = nullptr;
	if(m_driver == "SQLite") {
		dopts = CSLSetNameValue(dopts, "FORMAT", "SPATIALITE");
		dopts = CSLSetNameValue(dopts, "GEOMETRY_NAME", "geom");
	} else if(m_driver == "ESRI Shapefile") {
		dopts = CSLSetNameValue(dopts, "2GB_LIMIT", "YES");
	}

	m_layer = m_ds->CreateLayer(m_layerName.c_str(), sr, geomType(m_type), dopts);

	if(sr)
		sr->Release();

	if(dopts)
		CSLDestroy(dopts);

	if (!m_layer) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to create layer, " << m_layerName << ".");
	}

	for(const auto& it : m_fieldTypes) {
		OGRFieldDefn def(it.first.c_str(), fieldType(it.second));
		m_layer->CreateField(&def);
	}

	m_fdef = m_layer->GetLayerDefn();

	if (!m_fdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retrieve layer definition.");
	}

    OGRGeomFieldDefn* gdef = m_fdef->GetGeomFieldDefn(0);

	if (!gdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retrieve geometry field definition.");
	}

    m_geomName = std::string(gdef->GetNameRef());

}

DB::DB(const std::string& file, const std::string& layer, const std::string& driver,
		const std::unordered_map<std::string, FieldType>& fields,
		GeomType type, int srid, bool replace) :
    m_type(type),
    m_srid(srid),
    m_file(file),
	m_layerName(layer),
	m_driver(driver),
	m_fieldTypes(fields),
	m_ds(nullptr),
	m_layer(nullptr),
	m_fdef(nullptr) {

	GDALAllRegister();

	// If the driver was not given, try to discover it from an existing file, otherwise fail.
	if (m_driver.empty()) {
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(m_file.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr));
		if (!ds)
			g_runerr("Driver not given and no existing file to read from.");
		const char* drv = ds->GetDriverName();
		GDALClose(ds);
		if(!drv)
			g_runerr("Driver not given and no existing file to read from.");
		m_driver = drv;
		g_warn("Using driver " << drv << " from existing file.");
	}

	// Check layer name.
	if (m_layerName.empty()) {
		m_layerName = "data";
		g_warn("Using layer name data as none given.");
	}

	// If replace and file exists, delete the existing file.
    if(replace && isfile(file))
        rem(file);

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(m_driver.c_str());
    if(!drv)
        g_runerr("Driver not found for " << m_file << " (" << driver << ")");

	// If the file is sqlite, use the spatialite driver.
	char** dopts = nullptr;
	if(m_driver == "SQLite")
		dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

	m_ds = drv->Create(m_file.c_str(), 0, 0, 0, GDT_Unknown, dopts);

	if(dopts)
		CSLDestroy(dopts);

	if(!m_ds)
        g_runerr("Failed to create data set for " << m_file);

	OGRSpatialReference* sr = nullptr;
	if(m_srid) {
		sr = new OGRSpatialReference();
		sr->importFromEPSG(m_srid);
	}

	dopts = nullptr;
	if(m_driver == "SQLite") {
		dopts = CSLSetNameValue(dopts, "FORMAT", "SPATIALITE");
		dopts = CSLSetNameValue(dopts, "GEOMETRY_NAME", "geom");
	} else if(m_driver == "ESRI Shapefile") {
		dopts = CSLSetNameValue(dopts, "2GB_LIMIT", "YES");
	}

	m_layer = m_ds->CreateLayer(m_layerName.c_str(), sr, geomType(m_type), dopts);

	if(sr)
		sr->Release();

	if(dopts)
		CSLDestroy(dopts);

	if (!m_layer) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to create layer, " << m_layerName << ".");
	}

	for(const auto& it : m_fieldTypes) {
		OGRFieldDefn def(it.first.c_str(), fieldType(it.second));
		m_layer->CreateField(&def);
	}

	m_fdef = m_layer->GetLayerDefn();

	if (!m_fdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retrieve layer definition.");
	}

    OGRGeomFieldDefn* gdef = m_layer->GetLayerDefn()->GetGeomFieldDefn(0);

    if (!gdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retrieve geometry field definition.");
	}

    m_geomName = std::string(gdef->GetNameRef());
}

DB::DB(const std::string& file, const std::string& layer) :
    m_type(GeomType::GTUnknown),
    m_srid(0),
    m_file(file),
	m_layerName(layer),
	m_ds(nullptr),
	m_layer(nullptr),
	m_fdef(nullptr) {

	open();

}

void DB::open() {

    GDALAllRegister();

    m_ds = static_cast<GDALDataset*>(GDALOpenEx(m_file.c_str(), GDAL_OF_VECTOR|GDAL_OF_UPDATE, nullptr, nullptr, nullptr));
    if(!m_ds)
        g_runerr("Failed to open data set for " << m_file);

    if(!m_layerName.empty())
		m_layer = m_ds->GetLayerByName(m_layerName.c_str());
    if(!m_layer)
		m_layer = m_ds->GetLayer(0);
	if(!m_layer) {
		GDALClose(m_ds);
		m_ds = nullptr;
		if(m_layerName.empty()) {
			g_runerr("No layer, " << m_layerName << " was found on this data set, and no default was available.");
		} else {
			g_runerr("No layer was found on this data set.");
		}
	}

	m_type = geomType(m_layer->GetGeomType());

	OGRGeomFieldDefn* gdef = m_layer->GetLayerDefn()->GetGeomFieldDefn(0);
	if (!gdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retreive geometry field definition.");
	}

	m_geomName = std::string(gdef->GetNameRef());

	m_fdef = m_layer->GetLayerDefn();
	if (!m_fdef) {
		GDALClose(m_ds);
		m_ds = nullptr;
		g_runerr("Failed to retrieve layer definition.");
	}

	for(int i = 0; i < m_fdef->GetFieldCount(); ++i) {
		OGRFieldDefn* def = m_fdef->GetFieldDefn(i);
		m_fieldTypes[std::string(def->GetNameRef())] = fieldType(def->GetType());
	}
}

const std::string& DB::geomColumnName() const {
	return m_geomName;
}

void DB::flush() {
	if(OGRERR_NONE != m_layer->SyncToDisk())
		g_warn("Failed to sync to disk.");
}

void DB::close(bool remove) {
	if(m_ds) {
		flush();
		GDALClose(m_ds);
		m_ds = nullptr;
		if(remove)
			rem(m_file);
	}
}

DB::~DB() {
	close();
}

std::map<std::string, std::set<std::string> > DB::extensions() {
	GDALAllRegister();
	std::map<std::string, std::set<std::string> > extensions;
	GDALDriverManager* mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver* drv = mgr->GetDriver(i);
		if(!isRast(drv)) {
			const char* desc = drv->GetDescription();
			if(desc != nullptr) {
				const char *ext = drv->GetMetadataItem(GDAL_DMD_EXTENSION);
				if(ext != nullptr ) {
					std::list<std::string> lst;
					split(std::back_inserter(lst), std::string(ext));
					for(const std::string& item : lst)
						extensions[desc].insert("." + lowercase(item));
				}
			}
		}
	}
	return extensions;
}

std::map<std::string, std::string> DB::drivers() {
	std::vector<std::string> filter;
	return drivers(filter);
}

std::map<std::string, std::string> DB::drivers(const std::vector<std::string>& filter) {
	GDALAllRegister();
	std::map<std::string, std::string> drivers;
	GDALDriverManager *mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver *drv = mgr->GetDriver(i);
		if(!isRast(drv)) {
			const char* name = drv->GetMetadataItem(GDAL_DMD_LONGNAME);
			const char* desc = drv->GetDescription();
			if(name != nullptr && desc != nullptr) {
				bool found = true;
				if(!filter.empty()) {
					found = false;
					for(const std::string& f : filter) {
						if(f == desc) {
							found = true;
							break;
						}
					}
				}
				if(found)
					drivers[desc] = name;
			}
		}
	}
	return drivers;
}

std::string DB::getDriverForFilename(const std::string& filename) {
	std::string ext = extension(filename);
	std::map<std::string, std::set<std::string> > drivers = extensions();
	std::string result;
	for(const auto& it : drivers) {
		if(it.second.find(ext) != it.second.end())
			result = it.first;
	}
	return result;
}

void DB::clear() {
	g_runerr("Not implemented.");
}

void DB::setCacheSize(size_t) {
	g_runerr("Not implemented.");
}

/**
 * Join the names together as a comma-delimited string of quoted names.
 */
std::string nameJoin(std::unordered_set<std::string>& names) {
	std::stringstream ss;
	bool first = true;
	for(const std::string& n : names) {
		if(!first)
			ss << ",";
		ss << '"' << n << '"';
		first = false;
	}
	return ss.str();
}

void DB::convert(const std::string& filename, const std::string& driver) {

	// If replace and file exists, delete the existing file.
    if(isfile(filename))
        rem(filename);

    GDALAllRegister();

	{
		GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		if(!drv)
			g_runerr("Driver not found for " << filename << " (" << driver << ")");

		// If the file is sqlite, use the spatialite driver.
		char **dopts = nullptr;
		if(driver == "SQLite")
			dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

		GDALDataset* ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);

		if(dopts)
			CSLDestroy(dopts);

		if(!ds)
			g_runerr("Failed to create data set for " << filename);

		OGRLayer* layer = m_ds->GetLayerByName(m_layerName.c_str());
		if (!layer) {
			GDALClose(ds);
			g_runerr("Failed to get layer " << m_layerName << ".");
		}

		dopts = nullptr;
		if(driver == "SQLite") {
			dopts = CSLSetNameValue(dopts, "FORMAT", "SPATIALITE");
			dopts = CSLSetNameValue(dopts, "GEOMETRY_NAME", "geom");
		} else if(m_driver == "ESRI Shapefile") {
			dopts = CSLSetNameValue(dopts, "2GB_LIMIT", "YES");
		}

		OGRLayer* newLayer = ds->CopyLayer(layer, m_layerName.c_str(), dopts);

		if(dopts)
			CSLDestroy(dopts);

		if(!newLayer) {
			GDALClose(ds);
			g_runerr("Failed to copy layer to new database.");
		}

		ds->FlushCache();
		GDALClose(ds);
    }
}

void DB::dropGeomIndex(const std::string& table, const std::string& column) {
	if (!m_ds)
		g_runerr("No database open.");
	std::string _table;
	if(table.empty()) {
		_table = m_layerName;
	} else {
		_table = table;
	}
	std::string sql = "SELECT DisableSpatialIndex('" + _table + "', '" + column + "'); DropTable idx_" + _table + "_Geometry; VACUUM;";
	m_ds->ExecuteSQL(sql.c_str(), nullptr, nullptr);
}

void DB::createGeomIndex(const std::string& table, const std::string& column) {
	if (!m_ds)
		g_runerr("No database open.");
	std::string _table;
	if(table.empty()) {
		_table = m_layerName;
	} else {
		_table = table;
	}
	std::string sql = "SELECT CreateSpatialIndex('" + _table + "', '" + column + "');";
	m_ds->ExecuteSQL(sql.c_str(), nullptr, nullptr);
}

std::string typeStr(FieldType type) {
	switch(type) {
	case FieldType::FTDouble:
		return "DOUBLE PRECISION";
	case FieldType::FTInt:
		return "INTEGER";
	case FieldType::FTString:
		return "TEXT";
	default:
		g_runerr("Unknown or unimplemented field type: " << type);
	}
}

uint64_t DB::getGeomCount() const {
	return static_cast<uint64_t>(m_layer->GetFeatureCount(1));
}

void DB::execute(const std::string&) {
	g_runerr("Not implemented.");
}

void DB::begin() {
    if(CPLE_None != m_layer->StartTransaction())
        g_warn("Failed to start transaction.");
}

void DB::rollback() {
    if(CPLE_None != m_layer->RollbackTransaction())
        g_warn("Failed to roll back transaction.");
}

void DB::commit() {
    if(CPLE_None != m_layer->CommitTransaction())
        g_warn("Failed to commit transaction.");
    flush();
}

int DB::srid() const {
    return m_srid;
}

