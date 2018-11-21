// http://asprs.org/a/society/committees/standards/LAS_1_4_r13.pdf

#include <fstream>
#include <string>
#include <vector>

class Header {
public:
	unsigned short fileSourceId;
	unsigned short globalEncoding;
	unsigned int guid1;
	unsigned short guid2;
	unsigned short guid3;
	unsigned char guid[4];
	unsigned char versionMajor;
	unsigned char versionMinor;
	char systemIdentifier[32];
	char generatingSoftware[32];
	unsigned short fileCreationDoY;
	unsigned short fileCreationYear;
	unsigned short headerSize;
	unsigned int offsetToPointData;
	unsigned int numVLRs;
	unsigned char pointRecordFormat;
	unsigned short pointRecordLength;
	unsigned int legNumPointRecords;
	unsigned int legNumPointsByReturn[5];
	double xScale;
	double yScale;
	double zScale;
	double xOffset;
	double yOffset;
	double zOffset;
	double maxX;
	double minX;
	double maxY;
	double minY;
	double maxZ;
	double minZ;
	unsigned long  waveDataPacketRecordStart;
	unsigned long  extVLRStart;
	unsigned int numExtVLR;
	unsigned long  numPointRecords;
	unsigned long  numPointsByReturn[15];

	void read(std::ifstream& str) {
		
		std::vector<char> buf(4);
		if(str.read(buf.data(), 4)) {
			if(std::string(buf.data()) != "LASF")
				throw std::runtime_error("Not a LAS file.");
		}

		buf.reserve(22);
		if(str.read(buf.data(), 22)) {
			int offset = 0;
			char* data = (char*) buf.data();
			fileSourceId = (unsigned short) *(data + offset); offset += 2;
			globalEncoding = (unsigned short) *(data + offset); offset += 2;
			guid1 = (unsigned long) *(data + offset); offset += 4;
			guid2 = (unsigned short) *(data + offset); offset += 2;
			guid3 = (unsigned short) *(data + offset); offset += 2;
			std::memcpy(guid4, data + offset, 4); offset += 4;
			versionMajor = (unsigned char) *(data + offset); offset++;
			versionMinor = (unsigned char) *(data + offset); offset++;
			std::memcpy(systemIdentifier, data + offset, 32); offset += 32;
			std::memcpy(generatingSoftware, data + offset, 32); offset += 32;
			fileCreationDoY = (unsigned short) *(data + offset); offset += 2;
			fileCreationYear = (unsigned short) *(data + offset); offset += 2;
			headerSize = (unsigned short) *(data + offset); offset += 2;
		}
	}

}
