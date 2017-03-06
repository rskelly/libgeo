#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include "geo.hpp"

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/join.hpp>

#ifdef _MSC_VER
#include <float.h>
namespace std {

	inline bool isnan(double value) {
		return _isnan(value) == 1;
	}

}
#endif

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

namespace geo {

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

            bool contains(const geo::util::Bounds &b, int dims = 2) const;

            bool intersects(const geo::util::Bounds &b, int dims = 2) const;

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

            void extend(const geo::util::Bounds &b);

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
            bool m_remove;
            boost::interprocess::file_mapping *m_mapping;
            boost::interprocess::mapped_region *m_region;
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

            // Parse a string of numbers, commas and dashes, representing numbers and
            // ranges of numbers, into a list of doubles. The step parameter allows
            // the creation of ranges with fractional intermediate values.
            // e.g., 1-4,7,9-11 gives {1, 2, 3, 4, 7, 9, 10, 11}
            template <class T>
            static void parseFloatRanges(T iter, const std::string &str, double step = 1.0) {
                std::stringstream ss;
                double first = 0, second = 0;
                bool range = false;
                int i = 0;
                char c = str[i++];
                while (true) {
                    if (c == '-') {
                        range = true;
                        first = atof(ss.str().c_str());
                        ss.str(std::string());
                    } else if (c == ',' || c == '\0') {
                        if (!range) {
                            *iter = atof(ss.str().c_str());
                            ++iter;
                            ss.str(std::string());
                        } else {
                            second = atof(ss.str().c_str());
                            ss.str(std::string());
                            step = g_abs(step);
                            if (first > second) {
                                double tmp = second;
                                second = first;
                                first = tmp;
                            }
                            for (double j = first; j <= second; j += step) {
                                *iter = j;
                                ++iter;
                            }
                            range = false;
                        }
                        if (c == '\0')
                            break;
                    } else {
                        ss << c;
                    }
                    c = str[i++];
                }
            }

            template <class T>
            static void parseIntRanges(T iter, const std::string &str) {
                std::stringstream ss;
                int first = 0, second = 0;
                bool range = false;
                int i = 0;
                char c = str[i++];
                while (true) {
                    if (c == '-') {
                        range = true;
                        first = atoi(ss.str().c_str());
                        ss.str(std::string());
                    } else if (c == ',' || c == '\0') {
                        if (!range) {
                            iter = atoi(ss.str().c_str());
                            ++iter;
                            ss.str(std::string());
                        } else {
                            second = atoi(ss.str().c_str());
                            ss.str(std::string());
                            if (first > second) {
                                int tmp = second;
                                second = first;
                                first = tmp;
                            }
                            for (int j = first; j <= second; ++j) {
                                *iter = j;
                                ++iter;
                            }
                            range = false;
                        }
                        if (c == '\0')
                            break;
                    } else {
                        ss << c;
                    }
                    c = str[i++];
                }
            }

            // Split a comma-delimited string into a list of unique integers.
            template <class T>
            static void intSplit(T iter, const std::string &str, const std::string &delim = ",") {
                std::stringstream ss(str);
                std::string item;
                while (std::getline(ss, item, *(delim.c_str()))) {
                    *iter = atoi(item.c_str());
                    ++iter;
                }
            }

            // Return true if the integer is in the list, or the list is empty.
            template <class T, class U>
            static bool inList(T begin, T end, U value) {
                while(begin != end) {
                    if(value == *begin)
                        return true;
                    ++begin;
                } 
                return false;
            }

            // Split a string with the given delimiter
            template <class T>
            static void splitString(T iter, const std::string &str, const std::string &delim = ",") {
                std::stringstream ss(str);
                std::string item;
                while (std::getline(ss, item, *(delim.c_str()))) {
                    *iter = item;
                    ++iter;
                }
            }

            // Join the string with the given delimiter
            template <class T>
            static std::string join(T begin, T end, const std::string &delim = ",") {
                std::vector<std::string> lst;
                while(begin != end) {
                    lst.push_back(*begin);
                    ++begin;
                }
                return boost::algorithm::join(lst, delim);
            }

            // Lowercase the string.
            static std::string& lower(std::string &str);

            // Uppercase the string.
            static std::string& upper(std::string &str);

            // Lowercase and return a copy of the string.
            static std::string lower(const std::string &str);

            // Uppercase and return a copy of the string.
            static std::string upper(const std::string &str);

            // Move the file.
            static void copyfile(std::string &srcfile, std::string &dstfile);

            // Create a temporary file at the given root folder. If no root is given,
            // a relative path is created.
            static const std::string tmpFile(const std::string &root = "");

            // Returns true if the file exists.
            static bool exists(const std::string &name);

            // Returns true if the parent folder exists.
            static bool pathExists(const std::string &name);
            
            // Remove a file.
            static bool rm(const std::string &name);

            // Make a directory.
            static bool mkdir(const std::string &dir);

            // Get the file extension.
            static std::string extension(const std::string &filename);

            // Populates the list with the files contained in dir. If ext is specified, filters
            // the files by that extension (case-insensitive). If dir is a file, it is added to the list.
            // Returns the number of files found.
            template <class T>
            static size_t dirlist(T iter, const std::string &dir, const std::string &ext = std::string()) {
                using namespace boost::filesystem;
                using namespace boost::algorithm;
                int i = 0;
                if (is_regular_file(dir)) {
                    *iter = dir;
                    ++iter;
                    ++i;
                } else {
                    directory_iterator end;
                    directory_iterator di(dir);
                    for (; di != end; ++di) {
                        if (!ext.empty()) {
                            std::string p(di->path().string());
                            to_lower(p);
                            if (ends_with(p, ext)) {
                                *iter = p;
                                ++iter;
                                ++i;
                            }
                        } else {
                            *iter = di->path().string();
                            ++iter;
                            ++i;
                        }
                    }
                }
                return i;
            }

            // Returns a memory-mapped file.
            static std::unique_ptr<MappedFile> mapFile(const std::string &filename, 
                uint64_t size, bool remove = true);

        };

        class CRS {
        public:

        	std::string epsg2Proj4(int crs) const;
        	std::string epsg2WKT(int crs) const;
        };

    } // util

} // geo

#endif
