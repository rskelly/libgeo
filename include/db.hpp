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


        /**
         * Provides some easy database read/write methods.
         */
        class DB {
        protected:
            GeomType m_type;											///<! The geometry column type.
            int m_srid;													///<! The horizontal spatial reference ID.
            std::string m_file;											///<! The output filename.
            std::string m_layerName;									///<! The name of the table or layer.
            std::string m_driver;										///<! The OGR driver to use. Availability may differ from system to system.
            std::string m_geomName;										///<! The name of the geometry field.
            std::unordered_map<std::string, FieldType> m_fieldTypes;	///<! A map containing field names and their types.
            GDALDataset *m_ds;											///<! Pointer to dataset.
            OGRLayer *m_layer;											///<! Pointer to layer.
            OGRFeatureDefn *m_fdef;										///<! Pointer to an OGR vector feature.

        public:

            /**
             * Construct a database object.
             *
             * @param file The filename for the database. Presumably a connection string could be used here.
             * @param layer The layer or table name.
             * @param driver The OGR driver.
             * @param fields A map of field names and their types.
             * @param type The geometry type.
             * @param srid The horizontal spatial reference ID.
             * @param replace If true, removes and replaces an existing database.
             */
            DB(const std::string &file, const std::string &layer, const std::string &driver,
            		const std::unordered_map<std::string, FieldType> &fields,
            		GeomType type, int srid = 0, bool replace = false);

            /**
             * Construct a database object.
             *
             * @param file The filename for the database. Presumably a connection string could be used here.
             * @param layer The layer or table name.
             */
            DB(const std::string &file, const std::string &layer);

            virtual ~DB();

            /**
             * Returns a map with file extensions for keys, and a list of driver names as values.
             *
             * @return A map with file extensions for keys, and a list of driver names as values.
             */
            static std::map<std::string, std::set<std::string> > extensions();

            /**
             * Returns a map with driver short names as keys and long names for values.
             *
             * @return A map of drivers with short name as key and long name as value.
             */
            static std::map<std::string, std::string> drivers();

            /**
             * Returns a map with driver short names as keys and long names for values. The filter
             * contains short names of drivers. Only those names will be included if the filter is
             * not empty.
             *
             * @param filter A vector of short driver names.
             * @return A map of drivers with short name as key and long name as value.
             */
            static std::map<std::string, std::string> drivers(const std::vector<std::string>& filter);

            // Returns a vector driver for the given filename.
            static std::string getDriverForFilename(const std::string &filename);

            /**
             * Clear the database -- delete all rows.
             */
            virtual void clear();

            /**
             * Re-save the database to the new filename using the given driver.
             *
             * @param filename 	The filename for the new database. Will be overwritten if it exists.
             * @param driver	The database driver.
             */
            virtual void convert(const std::string& filename, const std::string& driver);

            virtual void setCacheSize(size_t size);

            void dropGeomIndex(const std::string& table = "", const std::string& column = "GEOMETRY");

            void createGeomIndex(const std::string& table = "", const std::string& column = "GEOMETRY");

            uint64_t getGeomCount() const;

            void execute(std::string &sql);

            void begin();

            void rollback();

            void commit();

            int srid() const;

            void flush();

            /**
             * Close the database. Implies flush. Called by destructor.
             *
             * @param remove If true, the file is deleted if it exists.
             */
            void close(bool remove = false);
        };

    } // db

} // geo

#endif

