#include "geo.hpp"

geo::Monitor defaultMonitor;		///<! A default monitor for when none is provided.

geo::Monitor* geo::getDefaultMonitor() {
	return &defaultMonitor;
}

int g__loglevel = 0;

void geo::loglevel(int x) {
	g__loglevel = x;
}

int geo::loglevel() {
	return g__loglevel;
}