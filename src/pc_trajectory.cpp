/*
 * pc_trajectory.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: rob
 */

#include <string>
#include <fstream>

#include "pc_trajectory.hpp"

using namespace geo::pc;

SBETReader::SBETReader(const std::string& filename) :
	m_lastTime(0),
	m_filename(filename) {

	m_stream.open(m_filename, std::ios::binary|std::ios::in);
}

void SBETReader::reset() {
	m_stream.seekg(0);
	m_lastTime = 0;
}

bool SBETReader::nextRecord(SBETRecord& record) {
	if(!m_stream.eof()) {
		m_stream.read((char*) &record, sizeof(SBETRecord));
		m_lastTime = record.time;
		return true;
	}
	return false;
}

double SBETReader::lastTime() const {
	return m_lastTime;
}

bool SBETReader::readToTime(double time, SBETRecord& record) {
	reset();
	while(nextRecord(record)) {
		if(record.time >= time)
			return true;
	}
	return false;
}


SMRSMGReader::SMRSMGReader(const std::string& filename) :
	m_lastTime(0),
	m_filename(filename) {

	m_stream.open(m_filename, std::ios::binary|std::ios::in);
}

void SMRSMGReader::reset() {
	m_stream.seekg(0);
	m_lastTime = 0;
}

bool SMRSMGReader::nextRecord(SMRSMGRecord& record) {
	if(!m_stream.eof()) {
		m_stream.read((char*) &record, sizeof(SMRSMGRecord));
		m_lastTime = record.time;
		return true;
	}
	return false;
}

double SMRSMGReader::lastTime() const {
	return m_lastTime;
}

bool SMRSMGReader::readToTime(double time, SMRSMGRecord& record) {
	reset();
	while(nextRecord(record)) {
		if(record.time >= time)
			return true;
	}
	return false;
}
