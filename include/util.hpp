#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <set>
#include <list>
#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/join.hpp>

#include "geotools.hpp"

#ifdef _MSC_VER
#include <float.h>
namespace std {

    inline bool isnan(double value) {
        return _isnan(value) == 1;
    }
    
}
#endif

namespace geotools {

    namespace util {

        // Provides access to an allocated buffer
        // which will be safely disposed of.
        class Buffer {
        public:
            void *buf;

            Buffer(uint64_t size) {
                buf = std::calloc(size, 1);
            }

            ~Buffer() {
                std::free(buf);
            }
        };

        // Provides methods for handling status callbacks.
        class Callbacks {
        public:
            virtual ~Callbacks() = 0;
            virtual void stepCallback(float status) const = 0;
            virtual void overallCallback(float status) const = 0;
            virtual void statusCallback(const std::string &msg) const = 0;
        };

		// Simple class for capturing status from utility functions.
		class Status {
		public:
			Callbacks *callbacks;
			float start, end;

			Status(Callbacks *callbacks, float start, float end);

			void update(float s);
		};

        class Point {
        public:
            double x, y, z;
            std::unordered_map<std::string, std::string> fields;
            Point();
            Point(double x, double y, double z = 0);
            Point(double x, double y, double z, const std::map<std::string, std::string> &fields);
        };

        class Bounds {
        private:
            double m_minx, m_miny, m_minz;
            double m_maxx, m_maxy, m_maxz;
        public:
            Bounds();

            Bounds(double minx, double miny, double maxx, double maxy);

            Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz);

            bool contains(double x, double y) const;

            bool contains(double x, double y, double z) const;

            bool contains(const geotools::util::Bounds &b, int dims = 2) const;

            bool intersects(const geotools::util::Bounds &b, int dims = 2) const;

            Bounds intersection(const Bounds &other) const;

            double minx() const;

            void minx(double minx);

            double miny() const;

            void miny(double miny);

            double minz() const;

            void minz(double minz);

            double maxx() const;

            void maxx(double maxx);

            double maxy() const;

            void maxy(double maxy);

            double maxz() const;

            void maxz(double maxz);

            double width() const;

            double height() const;

            double depth() const;

            int maxCol(double resolution) const;

            int maxRow(double resolution) const;

            int toCol(double x, double resolution) const;
            
            int toRow(double y, double resolution) const;
            
            double toX(int col, double resolution) const;

            double toY(int row, double resolution) const;

            void extend(const geotools::util::Bounds &b);

            void extendX(double x);

            void extendY(double y);

            void extendZ(double z);

            void extend(double x, double y);

            void extend(double x, double y, double z);

            void collapse(int dims = 2);

            double operator[](size_t pos) const;

            void snap(double resolution);

            std::string print() const;

            void print(std::ostream &str) const;
            
            std::string toString() const;
            
            void fromString(const std::string&);

            void align(double x, double y, double xres, double yres);
            
        };

        class Util;

        // Maintains a memory-mapped file, and gives access to the mapped data.
        class MappedFile {
            friend class Util;
        private:
            std::string m_filename;
            uint64_t m_size;
            boost::interprocess::file_mapping *m_mapping;
            boost::interprocess::mapped_region *m_region;
            bool m_remove;
        protected:
            MappedFile(const std::string &filename, uint64_t size, bool remove);
        public:
            void* data();
            uint64_t size();
            ~MappedFile();
        };

        // Provides utility methods for working with LiDAR data.
        class Util {
        public:

            // Computes the area of the triangle given by the coordinates.
            static double computeArea(double x1, double y1, double z1, 
                double x2, double y2, double z2, 
                double x3, double y3, double z3);

            static void parseRanges(std::set<double> &values, const char *str, double step = 1.0);

            static void parseRanges(std::set<int> &values, const char *str);

            // Split a comma-delimited string into a set of unique integers.
            static void intSplit(std::set<int> &values, const char *str);

            // Split a comma-delimited string into a set of unique integers.
            static void intSplit(std::list<int> &values, const char *val);

            // Split a comma-delimited string into a set of unique integers.
            static void intSplit(std::vector<int> &values, const char *str);

            static void intSplit(std::set<uint8_t> &values, const char *str);

            // Return true if the integer is in the set, or the set is empty.
            static bool inList(std::set<int> &values, int value);

            static bool inList(std::vector<int> &values, int value);

            static void splitString(const std::string &str, std::list<std::string> &lst);

            // TODO: Use back inserter.
            static void splitString(const std::string &str, std::vector<std::string> &lst);

            static std::string join(const std::vector<std::string> &lst, const std::string &delim);

            static std::string& lower(std::string &str);

            static std::string& upper(std::string &str);

            static std::string lower(const std::string &str);

            static std::string upper(const std::string &str);

            // Prints out a status message; a percentage representing current
            // of total steps.
            static void status(int current, int total);

            static void copyfile(std::string &srcfile, std::string &dstfile);

            // Load the samples from a csv file. The file must have x, y and z headers.
            static void loadXYZSamples(std::string &datafile, std::vector<std::tuple<double, double, double> > &samples);

            static void loadIDXYZSamples(std::string &datafile, std::vector<std::tuple<std::string, double, double, double> > &samples);
            
            static void status(int step, int of, const std::string &message = "", bool end = false);

            static const std::string tmpFile(const std::string &root);

            static const std::string tmpFile();

            // Returns true if the file exists.
            static bool exists(const std::string &name);

            // Returns true if the path exists, even if the file does not.
            static bool pathExists(const std::string &name);
            
            static bool rm(const std::string &name);

            static bool mkdir(const std::string &dir);

            static std::string extension(const std::string &filename);

            // Populates the vector with the files contained in dir. If ext is specified, filters
            // the files by that extension (case-insensitive). If dir is a file, it is added to the list.
            // Returns the number of files found.
            static size_t dirlist(const std::string &dir, std::vector<std::string> &files, 
                const std::string &ext = std::string());

            static std::unique_ptr<MappedFile> mapFile(const std::string &filename, 
                uint64_t size, bool remove = true);

        };

        class CRS {
        public:

        	std::string epsg2Proj4(int crs) const;
        	std::string epsg2WKT(int crs) const;
        };

    } // util

} // geotools

#endif
