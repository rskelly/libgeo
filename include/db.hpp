#ifndef __DB_HPP__
#define __DB_HPP__
#include <map>

#include "geo.hpp"
#include "util.hpp"

#include <ogr_core.h>
#include <ogr_feature.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <list>
#include <vector>
#include <set>
#include <unordered_map>

using namespace geo::util;

namespace geo {

    namespace db {

        enum GeomType {
            GTUnknown = 0,
            GTPoint = 1,
            GTLine = 2,
            GTPolygon = 3,
            GTMultiPoint = 4,
            GTMultiLine = 5,
            GTMultiPolygon = 6
        };

        enum FieldType {
            FTUnknown = 0,
            FTInt = 1,
            FTDouble = 2,
            FTString = 3,
            FTBlob = 4
        };


        // Provides some easy database read/write methods.
        class DB {
        protected:
            GeomType m_type;
            int m_srid;
            std::string m_file;
            std::string m_layerName;
            std::string m_driver;
            std::string m_geomName;
            std::unordered_map<std::string, FieldType> m_fieldTypes;
            GDALDataset *m_ds;
            OGRLayer *m_layer;
            OGRFeatureDefn *m_fdef;

        public:

            DB(const std::string &file, const std::string &layer, const std::string &driver,
            		const std::unordered_map<std::string, FieldType> &fields,
            		GeomType type, int srid = 0, bool replace = false);

            DB(const std::string &file, const std::string &layer);

            ~DB();

            // Returns a map with file extensions for keys, and a list of driver names
            // as values.
            static std::map<std::string, std::set<std::string> > extensions();

            // Returns a map with driver short names as keys and long names for values.
            static std::map<std::string, std::string> drivers();

            // Returns a vector driver for the given filename.
            static std::string getDriverForFilename(const std::string &filename);

            void clear();

            void setCacheSize(size_t size);

            void dropGeomIndex(const std::string& table = "", const std::string& column = "GEOMETRY");

            void createGeomIndex(const std::string& table = "", const std::string& column = "GEOMETRY");

            uint64_t getGeomCount() const;

            void execute(std::string &sql);

            void begin();

            void rollback();

            void commit();

            int srid() const;

        };

    } // db

} // geo

#endif

