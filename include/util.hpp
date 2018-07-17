#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <chrono>
#include <unordered_map>
#include <map>
#include <memory>
#include <random>

#ifdef _MSC_VER
#include <float.h>
namespace std {

	inline bool isnan(double value) {
		return _isnan(value) == 1;
	}

}
#endif

#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "geo.hpp"

namespace geo {

    namespace util {

    	// Provides a simple way to get a formatted
    	// time beginning at an arbitrary point.
    	class G_DLL_EXPORT Stopwatch {
    	private:
    		bool m_reset;
    		std::chrono::time_point<std::chrono::system_clock> m_start;
    		std::chrono::time_point<std::chrono::system_clock> m_stop;
    		bool m_running;

    	public:
    		// Create a stopwatch with time formatted as hh:mm:ss
    		Stopwatch();
    		// Start the stopwatch at its current time.
    		void start();
    		// Stop the stopwatch. Does not reset.
    		void stop();
    		// Reset the stopwatch to zero.
    		void reset();
    		// Return the formatted time.
    		std::string time();
    		// Return the number of milliseconds elapsed.
    		uint64_t millis();
    	};

        // Provides access to an allocated buffer
        // which will be safely disposed of.
        class G_DLL_EXPORT Buffer {
        public:
            void *buf;

            Buffer() : buf(nullptr) {}

            Buffer(size_t size) :
            	buf(nullptr) {
            	resize(size);
            }

            void resize(size_t size) {
            	if(buf) {
            		std::free(buf);
            		buf = nullptr;
            	}
            	if(size > 0) {
					buf = std::malloc(size);
					if (!buf)
						g_runerr("Failed to allocate buffer.");
            	}
            }

            ~Buffer() {
            	resize(0);
            }
        };

        // Provides methods for handling status callbacks.
        class G_DLL_EXPORT Callbacks {
        public:
            virtual ~Callbacks() {}
            virtual void stepCallback(float status) const = 0;
            virtual void overallCallback(float status) const = 0;
            virtual void statusCallback(const std::string& msg) const = 0;
        };

		// Simple class for capturing status from utility functions.
		class G_DLL_EXPORT Status {
		private:
			geo::util::Callbacks *m_callbacks;
			float m_start, m_end;
		public:
			/**
			 * Create a Status object to wrap the given callbacks. The Status object
			 * does not own the Callbacks.
			 *
			 * @param callbacks The Callbacks implementation instance.
			 * @param start The start value for the status counter, at zero.
			 * @param end The end value for the status counter at 100.
			 */
			Status(geo::util::Callbacks *callbacks, float start, float end);
			void update(float s, const std::string& msg = "");
			float start() const;
			float end() const;
			geo::util::Callbacks* callbacks() const;
		};

        class G_DLL_EXPORT Point {
        public:
            double x, y, z;
            std::unordered_map<std::string, std::string> fields;
            Point();
            Point(double x, double y, double z = 0);
            Point(double x, double y, double z, const std::map<std::string, std::string> &fields);
        };

        class G_DLL_EXPORT Bounds {
        private:
            double m_minx, m_miny;
            double m_maxx, m_maxy;
            double m_minz, m_maxz;
        public:
            Bounds();

            Bounds(double minx, double miny, double maxx, double maxy);

            Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz);

            void set(double minx, double miny, double maxx, double maxy, double minz = 0, double maxz = 0);

            void assign(const Bounds& bounds);

            bool contains(double x, double y) const;

            bool contains(double x, double y, double z) const;

            bool contains(const geo::util::Bounds &b, int dims = 2) const;

            bool intersects(const geo::util::Bounds &b, int dims = 2) const;

            Bounds intersection(const Bounds &other) const;

            // Enlarge dimensions so the bounds becomes a cube.
            // If the bounds is 2D, will become a square.
            void cube();

            double midx() const;

            double midy() const;

            double midz() const;

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

            double volume() const;

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

            void collapse(int dims = 3);

            double operator[](size_t pos) const;

            void snap(double resolution);

            std::string print() const;

            void print(std::ostream &str) const;
            
            std::string toString() const;
            
            void fromString(const std::string&);

            void align(double x, double y, double xres, double yres);
            
        };

        // Maintains a memory-mapped file, and gives access to the mapped data.
        class G_DLL_EXPORT MappedFile {
        private:
            uint64_t m_size;
			std::string m_name;
			std::unique_ptr<Buffer> m_data;
            boost::interprocess::mapped_region* m_region;
			boost::interprocess::shared_memory_object* m_shm;
			boost::interprocess::file_mapping* m_fm;
			bool m_fileBacked;
			bool m_deleteOnDestruct;
			bool m_fileExists; // True if the file was preexisting -- assumed to contain data.
			bool m_fileInited;

        public:

            /**
             * Create a memory-mapped object (possibly file-backed) with the given name and size (in bytes).
             * If the name is a filename and the file exists, the object assumes that this file contains
             * data and should be treated as an already-mapped file. If the name does not describe an extant file
             * a new file is created with that name. For anonymous memory segments, the name is not
             * treated as a file, and is just used to identify the segment.
             * @param name 				A name for the memory segment; a filename if the mapped memory is to be file-backed.
             * @param size 				The number of bytes to map.
             * @param fileBacked 		Set to true to cause the object to map memory to a file.
             * @param deleteOnDestruct 	If true, the mapped file will be deleted on destruction.
             */
			MappedFile(const std::string& name, uint64_t size, bool fileBacked = false, bool deleteOnDestruct = true);

            /**
             * Create a memory-mapped object (possibly file-backed) with an automatically-assigned name.
             * @param size 				The number of bytes to map.
             * @param fileBacked 		Set to true to cause the object to map memory to a file.
             * @param deleteOnDestruct 	If true, the mapped file will be deleted on destruction.
             */
			MappedFile(uint64_t size, bool fileBacked = false, bool deleteOnDestruct = true);

            /**
             * Create a memory-mapped object (possibly file-backed) with an automatically-assigned name.
             */
			MappedFile();

            /**
             * Initialize the object with the given size (in bytes).
             * @param name 				A name for the memory segment; a filename if the mapped memory is to be file-backed.
             * @param size 				The number of bytes to map.
             * @param fileBacked 		Set to true to cause the object to map memory to a file.
             * @param deleteOnDestruct 	If true, the mapped file will be deleted on destruction.
             */
			void init(size_t size, bool fileBacked, bool deleteOnDestruct = true);

            /**
             * Initialize the object with the given name and size (in bytes).
             * If the name is a filename and the file exists, the object assumes that this file contains
             * data and should be treated as an already-mapped file. If the name does not describe an extant file
             * a new file is created with that name. For anonymous memory segments, the name is not
             * treated as a file, and is just used to identify the segment.
             * @param name 				A name for the memory segment; a filename if the mapped memory is to be file-backed.
             * @param size 				The number of bytes to map.
             * @param fileBacked 		Set to true to cause the object to map memory to a file.
             * @param deleteOnDestruct 	If true, the mapped file will be deleted on destruction.
             */
			void init(const std::string& name, size_t size, bool fileBacked, bool deleteOnDestruct = true);

			/**
			 * Flush the memory.
			 */
			void flush();

			/**
			 * Return the name of the mapped memory segment. May be a filename for file-backed objects.
			 * @return The name of the mapped memory segment.
			 */
            const std::string& name() const;

            /**
             * Return the address of the mapped data.
             * @return A pointer into the mapped segment.
             */
            void* data();

            /**
             * Resize the mapped memory segment. Not currently shrinkable.
             * @param The size in bytes of the resized segment.
             */
            void reset(uint64_t size);

            // Return the page size.
            static size_t pageSize();
            static size_t fixSize(size_t size);

            // Return the size of mapped memory.
            uint64_t size() const;

            ~MappedFile();

            /**
             * Write the object to memory, resize if necessary.
             * @param item 		The item to write.
             * @param position	The offset to write to, relative to the start address of the segment.
             * @return True if the operation is successful.
             */
            template <class T>
            bool write(const T& item, uint64_t position) {
                if(size() < position + sizeof(T))
                    reset(position + sizeof(T));
                std::memcpy((char*) data() + position, &item, sizeof(T));
                return true;
            }

            /**
             * Write the object to memory, resize if necessary.
             * @param item 		The item to write.
             * @param position	The offset to write to, relative to the start address of the segment.
             * @return True if the operation is successful.
             */
            bool write(void* input, uint64_t position, uint64_t size);

            /**
             * Read the object from memory.
             * @param item 		The item to read.
             * @param position	The offset to from, relative to the start address of the segment.
             * @return True if the operation is successful.
             */
            template <class T>
            bool read(T& item, uint64_t position) {
                if(position + sizeof(T) > size())
                    return false;
                std::memcpy(&item, (char*) data() + position, sizeof(T));
                return true;
            }

            /**
             * Read the object from memory.
             * @param item 		The item to read.
             * @param position	The offset to from, relative to the start address of the segment.
             * @return True if the operation is successful.
             */
            bool read(void* output, uint64_t position, uint64_t size);

        };

        /**
         * A file-backed allocator for stl classes.
         */
        template<typename T>
        class G_DLL_EXPORT mmap_allocator {
        private:
			MappedFile m_file;
			std::string m_filename;
			bool m_inited;

        public:
			typedef T value_type;

			mmap_allocator(const std::string& filename) :
			  m_filename(filename),
			  m_inited(false) {
			}

			mmap_allocator() :
			  m_inited(false) {
			}

			T* allocate (size_t n) {
			  if(m_inited) {
				  m_file.reset(n);
			  } else {
				  m_inited = true;
				  if(m_filename.empty()) {
					  m_file.init(m_filename, n, true);
				  } else {
					  m_file.init(n, true);
				  }
			  }
			  return reinterpret_cast<T*>(m_file.data());
			}

			void deallocate (T* p, size_t n) {
			}
        };


		/**
		 * The radius is the number of pixels to search *not including* the center.
		 * Therefore even and odd inputs are allowed. Returns a vector of offsets
		 * from a centre position at 0, 0.
		 * T must be constructible using T(col, row).
		 */
		template <class T>
		std::vector<T> circularKernel(int outerRadius, int innerRadius = 0, bool includeCenter = false) {
			double outr = g_sq((double) outerRadius) + 1.0;
			double inr = g_sq((double) innerRadius);
			std::vector<T> offsets;
			for(int r = -outerRadius; r < outerRadius + 1; ++r) {
				for(int c = -outerRadius; c < outerRadius + 1; ++c) {
					double d0 = g_sq((double) c) + g_sq((double) r);
					if((d0 <= outr && d0 >= inr) || (includeCenter && r == 0 && c == 0)) {
						offsets.emplace_back(c, r);
						std::cerr << "x ";
					} else {
						std::cerr << "  ";
					}
				}
				std::cerr << "\n";
			}
			return std::move(offsets);
		}

		template <class T>
		std::vector<T> squareKernel(int size, bool includeCenter = false) {
			std::vector<T> offsets;
			for(int r = -size / 2; r < size / 2 + 1; ++r) {
				for(int c = -size / 2; c < size / 2 + 1; ++c) {
					if(includeCenter || !(r == 0 && c == 0))
						offsets.emplace_back(c, r);
				}
			}
			return std::move(offsets);
		}

		/**
         * Provides utility methods.
         */
        class G_DLL_EXPORT Util {
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
            static void parseFloatRanges(T iter, const std::string& str, double step = 1.0) {
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
            static void parseIntRanges(T iter, const std::string& str) {
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
            static void intSplit(T iter, const std::string& str, const std::string& delim = ",") {
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
            static void splitString(T iter, const std::string& str, const std::string& delim = ",") {
                std::stringstream ss(str);
                std::string item;
                while (std::getline(ss, item, *(delim.c_str()))) {
                    *iter = item;
                    ++iter;
                }
            }

            // Join the string with the given delimiter
            template <class T>
            static std::string join(T begin, T end, const std::string& delim = ",") {
                std::vector<std::string> lst;
                while(begin != end) {
                    lst.push_back(*begin);
                    ++begin;
                }
                return boost::algorithm::join(lst, delim);
            }

            // Return the system tempp directory.
            static std::string tmpDir();

            // Join the path in a system-appropriate way.
            static std::string pathJoin(const std::string& a, const std::string& b);

            // Lowercase the string.
            static std::string& lower(std::string& str);

            // Uppercase the string.
            static std::string& upper(std::string& str);

            // Lowercase and return a copy of the string.
            static std::string lower(const std::string& str);

            // Uppercase and return a copy of the string.
            static std::string upper(const std::string& str);

			// Generate an MD5 hash of the string.
			static std::string md5(const std::string& input);
			
			static std::string sha256(const std::string& input);

			static std::string sha256File(const std::string& file);

            // Move the file.
            static void copyfile(const std::string& srcfile, const std::string& dstfile);

            // Return the basename of the file.
            static std::string basename(const std::string& filename);

            static std::string filename(const std::string& filename);

            // Create a temporary file at the given root folder. If no root is given,
            // a relative path is created.
            static std::string tmpFile(const std::string& root = "");

            // Returns true if the file exists.
            static bool exists(const std::string& name);

            // Returns true if the parent folder exists.
            static bool pathExists(const std::string& name);

            // Return the number of bytes available on the device containing
            // the given path.
            static uint64_t diskSpace(const std::string& path);
            
            // Generate a random number in the range.
			static double random(double from, double to) {
				std::random_device r;
				std::uniform_real_distribution<double> unif(from, to);
				std::default_random_engine re(r());
				return unif(re);
			}

            // Remove a file.
            static bool rm(const std::string& name);

            // Return the file size in bytes.
            static uint64_t filesize(const std::string& name);

            // Make a directory.
            static bool mkdir(const std::string& dir);

            static bool isFile(const std::string& path);

            static bool isDir(const std::string& path);

            // Get the parent directory
            static std::string parent(const std::string& filename);

            // Get the file extension.
            static std::string extension(const std::string& filename);

            // Populates the list with the files contained in dir. If ext is specified, filters
            // the files by that extension (case-insensitive). If dir is a file, it is added to the list.
            // Returns the number of files found.
            template <class T>
            static size_t dirlist(T iter, const std::string& dir, const std::string& ext = std::string()) {
                using namespace boost::filesystem;
                using namespace boost::algorithm;
                using namespace boost::system;
                int i = 0;
                if (is_regular_file(dir)) {
                    *iter = dir;
                    ++iter;
                    ++i;
                } else {
                    error_code ec;
                    directory_iterator end;
                    directory_iterator di(dir, ec);
                    for (; di != end; ++di) {
                        if (!ext.empty()) {
                            std::string p(di->path().string());
                            std::string tmp = p;
                            to_lower(tmp);
                            if (ends_with(tmp, ext)) {
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
            static std::unique_ptr<MappedFile> mapFile(const std::string& filename,
                uint64_t size, bool remove = true);

        };

        class G_DLL_EXPORT CRS {
        public:
        	std::string epsg2Proj4(int crs) const;
        	std::string epsg2WKT(int crs) const;
        };

    } // util

} // geo

#endif
