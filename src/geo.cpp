#include "geo.hpp"

int g__loglevel = 0;

geo::Monitor defaultMonitor;		///<! A default monitor for when none is provided.

geo::Monitor* geo::getDefaultMonitor() {
	return &defaultMonitor;
}
