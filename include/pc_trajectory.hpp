/*
 * pc_trajectory.hpp
 *
 *  Created on: Feb 21, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_TRAJECTORY_HPP_
#define INCLUDE_PC_TRAJECTORY_HPP_

#include <string>
#include <fstream>

namespace geo {
namespace pc {

typedef struct {
	double time; //	seconds	double
	double latitude; //	radians	double
	double longitude; //	radians	double
	double altitude; //	meters	double
	double xVelocity; //	meters/second	double
	double yVelocity;	//meters/second	double
	double zVelocity; //	meters/second	double
	double roll; //	radians	double
	double pitch; //	radians	double
	double heading; //	radians	double
	double wanderAngle;//	radians	double
	double xAcceleration;//	meters/second2	double
	double yAcceleration;//	meters/second2	double
	double zAcceleration; //	meters/second2	double
	double xAngularRate; //	radians/second	double
	double yAngularRate; //	radians/second	double
	double zAngularRate; //	radians/second	double
} SBETRecord;

class SBETReader {
private:
	double m_lastTime;
	std::string m_filename;
	std::ifstream m_stream;

public:

	/**
	 * Construct an SBETReader using the given filename.
	 * @param filename A filename of an sbet trajectory file.
	 */
	SBETReader(const std::string& filename);

	/**
	 * Read the next record into the SBETRecord instance.
	 * @param  record An SBETRecord instance.
	 * @return True, if there was a record to be read. If false,
	 *         the contents of the record are useless.
	 */
	bool nextRecord(SBETRecord& record);

	/**
	 * Reset the reader to before the first record.
	 */
	void reset();

	/**
	 * Return the time of the last record read.
	 * @return The time of the last record read.
	 */
	double lastTime() const;

	/**
	 * Read until the first record after the given time. Use
	 * this method to locate the record after a given time, then
	 * use SBETReader::nextRecord thereafter.
	 * @param  time The time to search for.
	 * @param  record An SBETRecord to read into.
	 * @return True if there was a record after the given time.
	 */
	bool readToTime(double time, SBETRecord& record);

};

typedef struct {
	double time; //	seconds	double
	double latitude; //	radians	double
	double longitude; //	radians	double
	double altitude; //	meters	double
	double xVelocity; //	meters/second	double
	double yVelocity;	//meters/second	double
	double zVelocity; //	meters/second	double
	double roll; //	radians	double
	double pitch; //	radians	double
	double heading; //	radians	double
} SMRSMGRecord;

class SMRSMGReader {
private:
	double m_lastTime;
	std::string m_filename;
	std::ifstream m_stream;

public:

	/**
	 * Construct an SMRSMGReader using the given filename.
	 * @param filename A filename of an sbet trajectory file.
	 */
	SMRSMGReader(const std::string& filename);

	/**
	 * Read the next record into the SMRSMGRecord instance.
	 * @param  record An SMRSMGRecord instance.
	 * @return True, if there was a record to be read. If false,
	 *         the contents of the record are useless.
	 */
	bool nextRecord(SMRSMGRecord& record);

	/**
	 * Reset the reader to before the first record.
	 */
	void reset();

	/**
	 * Return the time of the last record read.
	 * @return The time of the last record read.
	 */
	double lastTime() const;

	/**
	 * Read until the first record after the given time. Use
	 * this method to locate the record after a given time, then
	 * use SMRSMGReader::nextRecord thereafter.
	 * @param  time The time to search for.
	 * @param  record An SMRSMGRecord to read into.
	 * @return True if there was a record after the given time.
	 */
	bool readToTime(double time, SMRSMGRecord& record);

};

} // pc
} // geo


#endif /* INCLUDE_PC_TRAJECTORY_HPP_ */
