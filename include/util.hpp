/*
 * contrem_util.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#ifndef INCLUDE_UTIL_HPP_
#define INCLUDE_UTIL_HPP_

#include "geo.hpp"

#include <array>
#include <cstring>
#include <chrono>

#include <gdal_priv.h>

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
G_DLL_EXPORT enum class FileType {
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
G_DLL_EXPORT enum class DataType {
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
G_DLL_EXPORT enum class Interleave {
	BIL,
	BSQ,
	BIP
};

/**
 * Normalization methods.
 */
G_DLL_EXPORT enum class NormMethod {
	ConvexHull,
	ConvexHullLongestSeg,
	Line,
	Unknown
};

G_DLL_EXPORT constexpr std::array<FileType, 3> OUTPUT_TYPES = {FileType::GTiff, FileType::ENVI, FileType::CSV};			///<! Allowed output types for results.

G_DLL_EXPORT constexpr std::array<NormMethod, 3> NORM_METHODS = {NormMethod::ConvexHull, NormMethod::ConvexHullLongestSeg, NormMethod::Line};

/**
 * A class that when instantiated creates a temporary file
 * whose lifecycle is automatically managed. When the class
 * destructs, the file is deleted. Maintains the filename
 * and file descriptor.
 */
class G_DLL_EXPORT TmpFile {
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


G_DLL_EXPORT FileType getFileType(const std::string& filename);

G_DLL_EXPORT std::string fileTypeAsString(FileType type);

G_DLL_EXPORT FileType fileTypeFromString(const std::string& type);

G_DLL_EXPORT NormMethod normMethodFromString(const std::string& method);

G_DLL_EXPORT std::string normMethodAsString(NormMethod method);

G_DLL_EXPORT bool isnonzero(const double& v);

/**
 * Return true if it's a dir and it exists.
 */
G_DLL_EXPORT bool isdir(const std::string& path);

/**
 * Return true if it's a file and it exists.
 */
G_DLL_EXPORT bool isfile(const std::string& path);

/**
 * Remove the directory or file.
 */
G_DLL_EXPORT bool rem(const std::string& dir);

/**
 * Attempt to return the system temp dir. Falls back to ".".
 */
G_DLL_EXPORT std::string gettmpdir();

/**
 * Return the processid.
 */
G_DLL_EXPORT int pid();

/**
 * Create a temporary directory and return the path.
 *
 * Attempts to put the directory in dir, otherwise in the system temp dir.
 *
 * The process ID is added to the prefix to create the dir name. Subsequent
 * calls from the same process with the same prefix return the same path.
 */
G_DLL_EXPORT std::string tmpdir(const std::string& prefix, const std::string& dir = "");

/**
 * Create a temporary file name. A random string is added to the end of the prefix.
 *
 * If dir is empty, creates a file in the system temp dir. Otherwise attempts to build the
 * directory if required and puts the file there.
 */
G_DLL_EXPORT std::string tmpfile(const std::string& prefix, const std::string& dir = "");

/**
 * \brief Return the size of the file in bytes.
 *
 * \param f The file path.
 * \return The size of the file in bytes.
 */
G_DLL_EXPORT size_t filesize(const std::string& f);

/**
 * \brief Renames or moves the file, across device boundaries if necessary.
 *
 * Uses rename for files on the same device. Throws an exception if the destination
 * is an extant directory; overwrites an extant file.
 *
 * \param from The existing file.
 * \param to The destination file.
 * \return True on success.
 */
G_DLL_EXPORT bool rename(const std::string& from, const std::string& to);

/**
 * Return the parent directory of the path.
 */
G_DLL_EXPORT std::string parent(const std::string& path);

/**
 * Recursively make the directory.
 */
G_DLL_EXPORT bool makedir(const std::string& filename);

/**
 * Join to paths together using the appropriate separator for the system.
 *
 * \param a The first path part.
 * \param b The second path part.
 * \return The joined path.
 */
G_DLL_EXPORT std::string join(const std::string& a, const std::string& b);

/**
 * Join the items represented by the iterator using the given delimited.
 *
 * \param begin The start iterator.
 * \param end The end iterator.
 * \param delim The delimiter.
 * \return The joined string.
 */
template <class Iter>
G_DLL_EXPORT std::string join(Iter begin, Iter end, const std::string& delim = ",") {
	std::stringstream ss;
	ss << *begin;
	++begin;
	while(begin != end) {
		ss << delim << *begin;
		++begin;
	}
	return ss.str();
}

/**
 * \brief Return the base name of the path. If it's a file with the extension, the extension (including separator) is removed.
 *
 * \param path A path.
 * \return The base name.
 */
G_DLL_EXPORT std::string basename(const std::string& path);

/**
 * \brief Return the file extension including separator (.) if there is one, else an empty string.
 *
 * \param A file path.
 * \return The file extension.
 */
G_DLL_EXPORT std::string extension(const std::string& path);

/**
 * Remove non-alphanumeric characters and replace with underscores.
 */
G_DLL_EXPORT std::string sanitize(const std::string& str);

/**
 * Return the byte size of the given GDAL type.
 *
 * \param The GDAL type.
 * \return The type size.
 */
G_DLL_EXPORT int gdalTypeSize(GDALDataType type);

/**
 * \brief Get the well-known text representation of the projection from the SRID.
 *
 * \param srid The SRID.
 * \return The projection string.
 */
G_DLL_EXPORT std::string projectionFromSRID(int srid);


/**
 * Convert the given char buffer to a buffer of the templated type.
 *
 * \param[in] The source buffer.
 * \param[out] The output buffer.
 */
template <class T>
G_DLL_EXPORT void inline convertBuffer(std::vector<char>& raw, std::vector<T>& buf) {
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
G_DLL_EXPORT void split(T iter, const std::string& str, const std::string& delim = ",") {
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
G_DLL_EXPORT std::string lowercase(const std::string& str);

/**
 * Convert the given raw char buffer to the typed output buffer
 * using the GDALDataType as a guide.
 *
 * \param type the GDALDataType.
 * \param[in] rawBuf The raw character buffer.
 * \param[out] buf The output buffer.
 */
template <class T>
G_DLL_EXPORT void convertBuffer(GDALDataType type, std::vector<char>& rawBuf, std::vector<T>& buf) {

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
G_DLL_EXPORT uint64_t morton(uint32_t x, uint32_t y);

G_DLL_EXPORT double random(double min, double max);

/**
 * Return the clock time in microseconds.
 *
 * \return The clock time in microseconds.
 */
G_DLL_EXPORT uint64_t microtime();

/**
 * \brief Get the distance from the point to the line.
 *
 * \param x The point x.
 * \param y The point y.
 * \param x1 The line end x.
 * \param y1 The line end y.
 * \param x2 The line end x.
 * \param y2 The line end y.
 * \return The distance.
 */
template <class T>
G_DLL_EXPORT T linedist(T x, T y, T x1, T y1, T x2, T y2) {
	return std::abs((y2 - y1) * x - (y2 - y1) * y + x2 * y1 - y2 * x1) /
		std::sqrt(geo::sq(y2 - y1) + geo::sq(x2 - x1));
}

template <class T> 
class G_DLL_EXPORT Bounds {
private:
	T m_minx, m_miny;
	T m_maxx, m_maxy;
	T m_minz, m_maxz;
public:
	Bounds() : Bounds(maxvalue<T>(), maxvalue<T>(),
		minvalue<T>(), minvalue<T>(),
		maxvalue<T>(), minvalue<T>()) {}

	Bounds(T minx, T miny, T maxx, T maxy) :
		Bounds(minx, miny, maxx, maxy, maxvalue<T>(), minvalue<T>()) {}

	Bounds(T minx, T miny, T maxx, T maxy, T minz, T maxz) :
		m_minx(minx), m_miny(miny), m_maxx(maxx), m_maxy(maxy), m_minz(minz), m_maxz(maxz) {}

	void set(T minx, T miny, T maxx, T maxy, T minz = 0, T maxz = 0) {
		m_minx = geo::min(minx, maxx);
		m_miny = geo::min(miny, maxy);
		m_minz = geo::min(minz, maxz);
		m_maxx = geo::max(minx, maxx);
		m_maxy = geo::max(miny, maxy);
		m_maxz = geo::max(minz, maxz);
	}

	void assign(const Bounds<T>& bounds) {
		set(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy(), bounds.minz(), bounds.maxz());
	}

	bool contains(T x, T y) const {
		return x >= m_minx && x < m_maxx && y >= m_miny && y < m_maxy;
	}

	bool contains(T x, T y, T z) const {
		return contains(x, y) && z >= m_minz && z < m_maxz;
	}

	bool contains(const geo::util::Bounds<T>& b, int dims = 2) const {
		if (dims == 3) {
			return contains(b.minx(), b.miny(), b.minz())
				&& contains(b.maxx(), b.maxy(), b.maxz());
		}
		else {
			return contains(b.minx(), b.miny()) && contains(b.maxx(), b.maxy());
		}
	}

	bool intersects(const geo::util::Bounds<T>& b, int dims = 2) const {
		if (dims == 3) {
			return !(b.maxx() < minx() || b.maxy() < miny() || b.minx() > maxx()
				|| b.miny() > maxy() || b.minz() > maxz() || b.maxz() < minz());
		}
		else {
			return !(b.maxx() < minx() || b.maxy() < miny() || b.minx() > maxx()
				|| b.miny() > maxy());
		}
	}
	
	Bounds intersection(const Bounds<T>& other) const {
		return Bounds(geo::max(minx(), other.minx()), geo::max(miny(), other.miny()),
			geo::min(maxx(), other.maxx()), geo::min(maxy(), other.maxy()));
	}



	void cube() {
		T max = geo::max(width(), height());
		if (m_maxz != maxvalue<T>()) {
			max = geo::max(max, depth()) / 2.0;
			set(midx() - max, midy() - max, midx() + max, midy() + max, midz() - max, midz() + max);
		}
		else {
			max /= 2.0;
			set(midx() - max, midy() - max, midx() + max, midy() + max);
		}
	}

	T midx() const {
		return (m_maxx + m_minx) / 2.0;
	}

	T midy() const {
		return (m_maxy + m_miny) / 2.0;
	}

	T midz() const {
		return (m_maxz + m_minz) / 2.0;
	}


	T minx() const {
		return m_minx;
	}

	void minx(T minx) {
		m_minx = minx;
	}

	T miny() const {
		return m_miny;
	}

	void miny(T miny) {
		m_miny = miny;
	}

	T minz() const {
		return m_minz;
	}

	void minz(T minz) {
		m_minz = minz;
	}

	T maxx() const {
		return m_maxx;
	}

	void maxx(T maxx) {
		m_maxx = maxx;
	}

	T maxy() const {
		return m_maxy;
	}

	void maxy(T maxy) {
		m_maxy = maxy;
	}

	T maxz() const {
		return m_maxz;
	}

	void maxz(T maxz) {
		m_maxz = maxz;
	}

	T width() const {
		return maxx() - minx();
	}

	T height() const {
		return maxy() - miny();
	}

	T depth() const {
		return maxz() - minz();
	}

	T volume() const {
		return width() * height() * depth();
	}

	int maxCol(T resolution) const {
		return (int)geo::abs(width() / resolution);
	}

	int maxRow(T resolution) const {
		return (int)geo::abs(height() / resolution);
	}

	int toCol(T x, T resolution) const {
		if (width() == 0.0)
			return 0;
		if (resolution > 0) {
			return (int)((x - m_minx) / width() * (width() / resolution));
		}
		else {
			return (int)((x - m_maxx) / width() * (width() / resolution));
		}
	}

	int toRow(T y, T resolution) const {
		if (height() == 0.0)
			return 0;
		if (resolution > 0) {
			return (int)((y - m_miny) / height() * (height() / resolution));
		}
		else {
			return (int)((y - m_maxy) / height() * (height() / resolution));
		}
	}

	T toX(int col, T resolution) const {
		if (resolution > 0) {
			return m_minx + resolution * col;
		}
		else {
			return m_maxx + resolution * col;
		}
	}

	T toY(int row, T resolution) const {
		if (resolution > 0) {
			return m_miny + resolution * row;
		}
		else {
			return m_maxy + resolution * row;
		}
	}

	void extend(const Bounds &b) {
		m_minx = geo::min(b.minx(), m_minx);
		m_maxx = geo::max(b.maxx(), m_maxx);
		m_miny = geo::min(b.miny(), m_miny);
		m_maxy = geo::max(b.maxy(), m_maxy);
		m_minz = geo::min(b.minz(), m_minz);
		m_maxz = geo::max(b.maxz(), m_maxz);
	}

	void extendX(T x) {
		m_minx = geo::min(x, m_minx);
		m_maxx = geo::max(x, m_maxx);
	}

	void extendY(T y) {
		m_miny = geo::min(y, m_miny);
		m_maxy = geo::max(y, m_maxy);
	}

	void extendZ(T z) {
		m_minz = geo::min(z, m_minz);
		m_maxz = geo::max(z, m_maxz);
	}

	void extend(T x, T y) {
		extendX(x);
		extendY(y);
	}

	void extend(T x, T y, T z) {
		extend(x, y);
		extendZ(z);
	}

	T operator[](size_t pos) const {
		switch (pos) {
		case 0:
			return m_minx;
		case 1:
			return m_miny;
		case 2:
			return m_maxx;
		case 3:
			return m_maxy;
		case 4:
			return m_minz;
		case 5:
			return m_maxz;
		default:
			g_argerr("Illegal position: " << pos);
		}
	}

	void snap(T resolution) {
		minx(std::floor(minx() / resolution) * resolution);
		miny(std::floor(miny() / resolution) * resolution);
		maxx(std::floor(maxx() / resolution) * resolution + resolution);
		maxy(std::floor(maxy() / resolution) * resolution + resolution);
	}

	void collapse(int dims) {
		minx(maxvalue<T>());
		miny(maxvalue<T>());
		maxx(minvalue<T>());
		maxy(minvalue<T>());
		if (dims == 3) {
			minz(maxvalue<T>());
			maxz(minvalue<T>());
		}
	}

	std::string print() const {
		std::stringstream s;
		print(s);
		return s.str();
	}

	void print(std::ostream &str) const {
		str << "[Bounds: " << minx() << ", " << miny() << ", " << minz() << "; "
			<< maxx() << ", " << maxy() << ", " << maxz() << "]";
	}

	void fromString(const std::string &str) {
		std::vector<std::string> parts;
		split(std::back_inserter(parts), str);
		if (parts.size() < 4)
			g_runerr("Bounds string must be 4 or 6 comma-separated values.");
		m_minx = atof(parts[0].c_str());
		m_miny = atof(parts[1].c_str());
		m_maxx = atof(parts[2].c_str());
		m_maxy = atof(parts[3].c_str());
		if (parts.size() >= 6) {
			m_maxz = atof(parts[4].c_str());
			m_maxz = atof(parts[5].c_str());
		}
	}

	std::string toString() const {
		std::stringstream ss;
		ss << minx() << "," << miny() << "," << maxx() << "," << maxy() << ","
			<< minz() << "," << maxz();
		return ss.str();
	}

    void align(T x, T y, T xres, T yres) {
		xres = geo::abs(xres);
		yres = geo::abs(yres);
		while (x < m_minx)
			x += xres;
		while (x > m_minx)
			x -= xres;
		m_minx = x;
		while (x < m_maxx)
			x += xres;
		m_maxx = x;
		while (y < m_miny)
			y += yres;
		while (y > m_miny)
			y -= yres;
		m_miny = y;
		while (y < m_maxy)
			y += yres;
		m_maxy = y;
	}


};




class G_DLL_EXPORT BivariateSpline {
private:
	std::vector<double> m_tx;
	std::vector<double> m_ty;
	std::vector<double> m_c;
	int m_nx;
	int m_ny;

public:

	/**
	 * \brief Initialize the spline.
	 *
	 * The x, y, z and weights lists are the same length.
	 *
	 * Recommended values for s depend on the weights w(i). if these are
	 * taken as 1/d(i) with d(i) an estimate of the standard deviation of
	 * z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
	 * sqrt(2*m)). if nothing is known about the statistical error in z(i)
	 * each w(i) can be set equal to one and s determined by trial and
	 * error, taking account of the comments above. the best is then to
	 * start with a very large value of s ( to determine the least-squares
	 * polynomial and the corresponding upper bound fp0 for s) and then to
	 * progressively decrease the value of s ( say by a factor 10 in the
	 * beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
	 * approximation shows more detail) to obtain closer fits.
	 *
	 * \param[inout] smooth The smoothness parameter. Recommended to be between m - sqrt(2 * m) and m + sqrt(2 * m) if the weights are 1/std(z). On exit has the smoothed value that was used.
	 * \param x The list of x-coordinates.
	 * \param y The list of x-coordinates.
	 * \param z The list of x-coordinates.
	 * \param[inout] weights The list of weights. If empty, estimated as 1/stddev(z). If this is the case, smooth should be chosen accordingly.
	 * \param x0 The minimum x coordinate of the computation area.
	 * \param y0 The minimum y coordinate of the computation area.
	 * \param x1 The maximum x coordinate of the computation area.
	 * \param y1 The maximum y coordinate of the computation area.
	 */
	int init(double& smooth, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
			std::vector<double>& weights,
			double x0, double y0, double x1, double y1);

	/**
	 * \brief Evaluate.
	 */
	int evaluate(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& z);

	int evaluate(const double* x, int nx, const double* y, int ny, double* z, int nz);

	double stddev(const std::vector<double>& v) const;

};

class G_DLL_EXPORT SmoothingSpline {
private:
	std::vector<double> m_tx;
	std::vector<double> m_c;
	int m_nx;

public:

	/**
	 * \brief Initialize the spline.
	 *
	 * The x, y, z and weights lists are the same length.
	 *
	 * Recommended values for s depend on the weights w(i). if these are
	 * taken as 1/d(i) with d(i) an estimate of the standard deviation of
	 * z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
	 * sqrt(2*m)). if nothing is known about the statistical error in z(i)
	 * each w(i) can be set equal to one and s determined by trial and
	 * error, taking account of the comments above. the best is then to
	 * start with a very large value of s ( to determine the least-squares
	 * polynomial and the corresponding upper bound fp0 for s) and then to
	 * progressively decrease the value of s ( say by a factor 10 in the
	 * beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
	 * approximation shows more detail) to obtain closer fits.
	 *
	 * \param[inout] smooth The smoothness parameter. Recommended to be between m - sqrt(2 * m) and m + sqrt(2 * m) if the weights are 1/std(z). On exit has the smoothed value that was used.
	 * \param x The list of x-coordinates.
	 * \param y The list of x-coordinates.
	 * \param[inout] weights The list of weights. If empty, estimated as 1/stddev(z). If this is the case, smooth should be chosen accordingly.
	 * \param x0 The minimum x coordinate of the computation area.
	 * \param x1 The maximum x coordinate of the computation area.
	 */
	int init(double& smooth, const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& weights,
			double x0, double x1);

	/**
	 * \brief Evaluate.
	 */
	int evaluate(const std::vector<double>& x, std::vector<double>& y);

	double stddev(const std::vector<double>& v) const;

	const std::vector<double>& knots() const;

	const std::vector<double>& coefficients() const;

};

namespace csv {

	enum G_DLL_EXPORT CSVType {
		Int,
		Double,
		String
	};

	class G_DLL_EXPORT CSVValue {
	public:
		std::string s;
		CSVType t;
		union {
			int i;
			double d;
		};
		int asInt() const;
		double asDouble() const;
		const std::string& asString() const;
	};

	class G_DLL_EXPORT CSVColumn {
	public:
		std::string name;
		CSVType type;
		std::vector<CSVValue> values;
	};

	class G_DLL_EXPORT CSV {
	private:
		std::vector<CSVColumn> m_values;
		std::string m_file;

		/**
		 * \brief Return true if the string contains a double.
		 *
		 * \param s A string, possibly representing a double.
		 * \return True if the string contains a double.
		 */
		bool isdouble(const std::string& s);

		/**
		 * \brief Return true if the string contains an int.
		 *
		 * \param s A string, possibly representing an int.
		 * \return True if the string contains an int.
		 */
		bool isint(const std::string& s);

	public:

		/**
		 * \brief Construct and possibly load a CSV file.
		 *
		 * \param file The filename.
		 * \param header True if the first line contains column names.
		 */
		CSV(const std::string& file = "", bool header = false);

		/**
		 * \brief Load a CSV file.
		 *
		 * \param file The filename.
		 * \param header True if the first line contains column names.
		 */
		void load(const std::string& file, bool header);

		/**
		 * \brief Return the list of column names.
		 *
		 * \return The list of column names.
		 */
		std::vector<std::string> columnNames() const ;

		/**
		 * \brief Return the row.
		 *
		 * \param idx The row index.
		 * \return The row of values.
		 */
		std::vector<CSVValue> row(size_t idx) const;

		/**
		 * \brief Return the column.
		 *
		 * \param name The column name.
		 * \return The column.
		 */
		std::vector<CSVValue> column(const std::string& name) const;

		/**
		 * \brief Return the column type.
		 *
		 * \param name The column name.
		 * \return The column type.
		 */
		CSVType columnType(const std::string& name) const;

		/**
		 * \brief Return the column.
		 *
		 * \param idx The column index.
		 * \return The column.
		 */
		std::vector<CSVValue> column(size_t i) const;

		/**
		 * \brief Return the column type.
		 *
		 * \param name The column index.
		 * \return The column type.
		 */
		CSVType columnType(size_t i) const;

	};

} // csv

/**
 * \brief A quick and dirty way to save a raster.
 *
 * The grid is a list of values in row-major order.
 *
 * \param file The path to the new files.
 * \param grid The grid of values to write to a raster.
 * \param cols The column width of teh raster.
 * \param rows The row height of the raster.
 * \param minx The x-coordinate at the top-left corner of the raster.
 * \param miny The y-coordinate at the top-left corner of the raster.
 * \param xres The resolution of the raster in x.
 * \param yres The resolution of the raster in y.
 * \param proj The projection.
 */
G_DLL_EXPORT void saveGrid(const std::string& file, const std::vector<double> grid,
		int cols, int rows, double minx, double miny,
		double xres, double yres, const std::string& proj);

/**
 * \brief Does what it says on the tin: keeps time.
 */
class G_DLL_EXPORT Stopwatch {
private:
	std::chrono::steady_clock::time_point m_begin;

public:
	void reset();
	void start();
	std::string time();

};

} // util
} // geo

#endif /* INCLUDE_UTIL_HPP_ */
