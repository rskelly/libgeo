/*
 * contrem_util.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#ifndef INCLUDE_UTIL_HPP_
#define INCLUDE_UTIL_HPP_

#include <array>
#include <cstring>

#include <gdal_priv.h>

#include "geo.hpp"

namespace {

	/**
	 * Take the value and split its bits apart by 2.
	 * This is not fast, but it's adaptable.
	 *
	 * \param x The value to split.
	 * \return The split value.
	 */
	template <class T>
	inline uint64_t bitsplit2(T x) {
		// https://lemire.me/blog/2018/01/09/how-fast-can-you-bit-interleave-32-bit-integers-simd-edition/
	    uint64_t word = (uint32_t) x;
	    word = (word ^ (word << 16)) & 0x0000ffff0000ffff;
	    word = (word ^ (word << 8 )) & 0x00ff00ff00ff00ff;
	    word = (word ^ (word << 4 )) & 0x0f0f0f0f0f0f0f0f;
	    word = (word ^ (word << 2 )) & 0x3333333333333333;
	    word = (word ^ (word << 1 )) & 0x5555555555555555;
	    return word;
	}

} // anon


namespace geo {
namespace util {

/**
 * Input and output file types.
 */
enum class FileType {
	GTiff,
	ENVI,
	ROI,
	SHP,
	CSV,
	SQLITE,
	Unknown
};

/**
 * Enumeration containing common data types.
 */
/**
 * The allowable types for a raster.
 */
enum class DataType {
	Float64 = 7,
	Float32 = 6,
	UInt32 = 5,
	UInt16 = 4,
	Byte = 3,
	Int32 = 2,
	Int16 = 1,
	None = 0
};

/**
 * Interleave methods.
 */
enum class Interleave {
	BIL,
	BSQ,
	BIP
};

/**
 * Normalization methods.
 */
enum class NormMethod {
	ConvexHull,
	ConvexHullLongestSeg,
	Line,
	Unknown
};

constexpr std::array<FileType, 3> OUTPUT_TYPES = {FileType::GTiff, FileType::ENVI, FileType::CSV};			///<! Allowed output types for results.

constexpr std::array<NormMethod, 3> NORM_METHODS = {NormMethod::ConvexHull, NormMethod::ConvexHullLongestSeg, NormMethod::Line};

/**
 * A class that when instantiated creates a temporary file
 * whose lifecycle is automatically managed. When the class
 * destructs, the file is deleted. Maintains the filename
 * and file descriptor.
 */
class TmpFile {
public:
	std::string filename;	///<! The filename of the temporary file.
	int fd;					///<! The file descriptor.
	size_t size;			///<! The file size.

	/**
	 * Create a file with the given size. The contents of the file are not defined.
	 *
	 * \param size The file size.
	 */
	TmpFile(size_t size = 0);

	/**
	 * Resize the file.
	 *
	 * \param size The size.
	 */
	void resize(size_t size);

	/**
	 * Close the file.
	 */
	void close();

	/**
	 * Destroy. Closes and deletes the file.
	 */
	~TmpFile();
};


FileType getFileType(const std::string& filename);

std::string fileTypeAsString(FileType type);

FileType fileTypeFromString(const std::string& type);

NormMethod normMethodFromString(const std::string& method);

std::string normMethodAsString(NormMethod method);

bool isnonzero(const double& v);

/**
 * Return true if it's a dir and it exists.
 */
bool isdir(const std::string& path);

/**
 * Return true if it's a file and it exists.
 */
bool isfile(const std::string& path);

/**
 * Remove the directory or file.
 */
bool rem(const std::string& dir);

/**
 * Return the parent directory of the path.
 */
std::string parent(const std::string& path);

/**
 * Recursively make the directory.
 */
bool makedir(const std::string& filename);

/**
 * Join to paths together using the appropriate separator for the system.
 *
 * \param a The first path part.
 * \param b The second path part.
 * \return The joined path.
 */
std::string join(const std::string& a, const std::string& b);

/**
 * Join the items represented by the iterator using the given delimited.
 *
 * \param begin The start iterator.
 * \param end The end iterator.
 * \param delim The delimiter.
 * \return The joined string.
 */
template <class Iter>
std::string join(Iter begin, Iter end, const std::string& delim = ",") {
	std::stringstream ss;
	ss << *begin;
	++begin;
	while(begin != end) {
		ss << delim << *begin;
		++begin;
	}
	return ss.str();
}

std::string basename(const std::string& path);

std::string extension(const std::string& path);

/**
 * Remove non-alphanumeric characters and replace with underscores.
 */
std::string sanitize(const std::string& str);

/**
 * Return the byte size of the given GDAL type.
 *
 * \param The GDAL type.
 * \return The type size.
 */
int gdalTypeSize(GDALDataType type);

/**
 * Convert the given char buffer to a buffer of the templated type.
 *
 * \param[in] The source buffer.
 * \param[out] The output buffer.
 */
template <class T>
void inline convertBuffer(std::vector<char>& raw, std::vector<T>& buf) {
	buf.reserve(raw.size() / sizeof(T));
	std::memcpy(buf.data(), raw.data(), raw.size());
}

/**
 * Split a string with the given delimiter
 *
 * \param An iterator to send chunks to.
 * \param str The source string.
 * \param delim The delimiter.
 */
template <class T>
void split(T iter, const std::string& str, const std::string& delim = ",") {
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, *(delim.c_str()))) {
        *iter = item;
        ++iter;
    }
}

/**
 * Return a lower case version of a string.
 *
 * \param str A string.
 */
std::string lowercase(const std::string& str);

/**
 * Convert the given raw char buffer to the typed output buffer
 * using the GDALDataType as a guide.
 *
 * \param type the GDALDataType.
 * \param[in] rawBuf The raw character buffer.
 * \param[out] buf The output buffer.
 */
template <class T>
void convertBuffer(GDALDataType type, std::vector<char>& rawBuf, std::vector<T>& buf) {

	switch(type) {
	case GDT_Float32:
	{
		std::vector<float> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Int32:
	{
		std::vector<int32_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_UInt32:
	{
		std::vector<uint32_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Int16:
	{
		std::vector<int16_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_UInt16:
	{
		std::vector<uint16_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Float64:
	{
		std::vector<double> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Byte:
		std::copy(rawBuf.begin(), rawBuf.end(), buf.begin());
		break;
	default:
		throw std::runtime_error("Unknown GDAL data type: " + std::to_string((int) type));
	}
}

/***
 * Calculate the 2-dimensional Morton coordinate of the point.
 * The scale is a multiplier that should be 10^n. This converts
 * floating-point coordinates to integers with reasonable precision.
 * Integers are 32 bits. The morton code is 64 bits.
 *
 * \param x The x coordinate scaled to the full value of an unsigned int.
 * \param y The y coordinatescaled to the full value of an unsigned int.
 * \return The 1-dimensional Morton or z-order coordinate.
 */
uint64_t morton(uint32_t x, uint32_t y);

double random(double min, double max);

/**
 * Return the clock time in microseconds.
 *
 * \return The clock time in microseconds.
 */
uint64_t microtime();

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



} // util
} // geo

#endif /* INCLUDE_UTIL_HPP_ */
