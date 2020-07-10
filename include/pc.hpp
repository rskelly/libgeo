/*
 * pc_io.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_IO_HPP_
#define INCLUDE_PC_IO_HPP_

#include <string>
#include <vector>

#include "geo.hpp"
#include "util.hpp"

namespace geo {
namespace pc {

	/**
	 * A point class used to retrieve information from Reader
	 * and provide information to Writer classes.
	 */
	class Point {
	protected:
		double m_x; ///<! The x coordinate.
		double m_y; ///<! The y coordinate.
		double m_z; ///<! The z coordinate.
		int m_cls;	///<! The point classification.

	public:
		/**
		 * \brief Get the x coordinate.
		 * \return The x coordinate.
		 */
		double x() const;
		/**
		 * \brief Get the y coordinate.
		 * \return The y coordinate.
		 */
		double y() const;
		/**
		 * \brief Get the z coordinate.
		 * \return The z coordinate.
		 */
		double z() const;

		/**
		 * \brief Get the point ASPRS classification.
		 * \return The point ASPRS classification.
		 */
		int cls() const;
	};

	/**
	 * A virtual class for handling individual point cloud files.
	 */
	class PointFile {
	private:
		size_t m_size;				///<! The point cloud size.
		geo::util::Bounds m_bounds;	///<! The point cloud bounding box.

	protected:

		/**
		 * \brief Scan the file for bounds and count if the header isn't trusted.
		 */
		virtual void scan() = 0;

	public:

		/**
		 * \brief Return the number of points in the file.
		 * \return The number of points in the file.
		 */
		virtual size_t size() const = 0;

		/**
		 * \brief Return the bounds of the point cloud.
		 * \return The bounds of the point cloud.
		 */
		virtual const geo::util::Bounds& bounds() const = 0;

		/**
		 * \brief Return the WKT projection string.
		 * \return The WKT projection string.
		 */
		virtual std::string projection() const;

		/**
		 * \brief Populate the given point reference with the next available point data.
		 * \param A reference to a point to populate.
		 * \return True if a point was found, false otherwise.
		 */
		virtual bool next(geo::pc::Point& pt) = 0;

		/**
		 * \brief Reset the object so it can read points again.
		 */
		virtual void reset();

		virtual ~PointFile() {}
	};

	/**
	 * A point cloud reader. Provides per-point and per-cloud
	 * access to underlying data.
	 */
	class Reader {
	private:
		size_t m_size;						///<! The number of points in the point cloud.
		size_t m_fileIdx;					///<! The index of the file currently being read.
		bool m_trustHeader;					///<! Set to true if the file headers are to be trusted.
		geo::util::Bounds m_bounds;			///<! The bounding box of the point cloud.
		std::string m_projection;			///<! The WKT projection.
		std::vector<std::string> m_files;	///<! The list of point cloud files.

		std::unique_ptr<PointFile> m_currentFile;	///<! Used to store the current open point cloud file.

		Reader();

		void init(const std::vector<std::string>& files, bool trustHeader);

	public:

		/**
		 * \brief Construct a point cloud reader.
		 *
		 * \param file The path to the point file.
		 * \param trustHeader If true, the bounds and count from
		 *                    the header are used. If false, the
		 *                    point cloud is scanned.
		 */
		Reader(const std::string& file, bool trustHeader);

		/**
		 * \brief Construct a point cloud reader.
		 *
		 * \param files A list of paths to point files.
		 * \param trustHeader If true, the bounds and count from
		 *                    the header are used. If false, the
		 *                    point cloud is scanned.
		 */
		Reader(const std::vector<std::string>& files, bool trustHeader);

		/**
		 * \brief Return the number of points in the point cloud.
		 *
		 * \return The number of points in the point cloud.
		 */
		size_t size() const;

		/**
		 * \brief Return the bounding box of the point cloud.
		 *
		 * \return The bounding box of the point cloud.
		 */
		const geo::util::Bounds& bounds() const;

		/**
		 * \brief Return the projection string.
		 *
		 * \return The projection string.
		 */
		const std::string& projection() const;

		/**
		 * \brief Read the next point in the point cloud.
		 *
		 * If no point is available, returns false, otherwise true.
		 *
		 * \param pt A reference to a point to populate.
		 * \return True if a point has been read successfully.
		 */
		bool next(geo::pc::Point& pt);

		/**
		 * \brief Reset the object so it can read points again.
		 */
		void reset();

	};

	class Writer {

	};

} // pc
} // geo


#endif /* INCLUDE_PC_IO_HPP_ */
